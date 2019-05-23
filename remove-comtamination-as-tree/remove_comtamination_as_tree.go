package main

import (
	"encoding/csv"
	"github.com/favorov/nwwreject"
	"time"
)



//we use https://gist.github.com/quwubin/fdf9a9b40f4c4fbbeb02
import (
	"flag"
	"math"
	"bytes"
	"fmt"
	"log"
	"os"
	"strings"
  "bufio"
  "io"
)


type fasta struct {
  id string
  desc string
  seq string
}

func build_fasta(header string, seq bytes.Buffer) (record fasta) {
  fields := strings.SplitN(header, " ", 2)

  if len(fields) > 1 {
    record.id = fields[0]
    record.desc = fields[1]
  }else{
    record.id = fields[0]
    record.desc = ""
  }

  record.seq = seq.String()

  return record
}

func parse(fastaFh io.Reader) chan fasta {

  outputChannel := make(chan fasta)

  scanner := bufio.NewScanner(fastaFh)
  // scanner.Split(bufio.ScanLines)
  header := ""
  var seq bytes.Buffer

  go func() {
    // Loop over the letters in inputString
    for scanner.Scan() {
      line := strings.TrimSpace(scanner.Text())
      if len(line) == 0 {
        continue
      }

      // line := scanner.Text()

      if line[0] == '>' {
        // If we stored a previous identifier, get the DNA string and map to the
        // identifier and clear the string
        if header != "" {
          // outputChannel <- build_fasta(header, seq.String())
          outputChannel <- build_fasta(header, seq)
          // fmt.Println(record.id, len(record.seq))
          header = ""
          seq.Reset()
        }

        // Standard FASTA identifiers look like: ">id desc"
        header = line[1:]
      } else {
        // Append here since multi-line DNA strings are possible
        seq.WriteString(line)
      }

    }

    outputChannel <- build_fasta(header, seq)

    // Close the output channel, so anything that loops over it
    // will know that it is finished.
    close(outputChannel)
  }()

  
  return outputChannel
}

//we use https://gist.github.com/quwubin/fdf9a9b40f4c4fbbeb02
//finish of the gist



func read_all_targets_from_fasta(fn string) map[string]string {
	// Open the fastq file specified on the command line
	// for reading:
	fh, err := os.Open(fn)
	// Check for open errors and abort:
	if err != nil {
		log.Fatalf("Error opening file: %v", err)
	}
	defer fh.Close()
	alltargtets:= make(map[string]string)
	for record := range parse(fh) {
		alltargtets[record.seq]=record.id
	}
	return alltargtets
}

func read_all_contaminators(fn string) []string {
	f, err := os.Open(fn)
	if err != nil {
		log.Fatalf("Error opening file: %v", err)
	}
	defer f.Close()

	contaminators:=make([]string,0); 
	rdr := csv.NewReader(bufio.NewReader(f))
	rdr.Comma = '\t'
	rdr.Read() //skip header line
	for {
		record, err := rdr.Read()
		if err != nil {
			if err == io.EOF {
				break
			}
			log.Fatal(err)
		}
		contaminators=append(contaminators,record[0])	
	}
	return contaminators
}

type relation struct {
  id1 string
  id2 string
}


//return is did we add new or just add a 
func add_to_contaminatoin_targets (contaminatoin_targets map[string]string, contamination_links map[relation]bool, seq, id, parent_id string) bool {
	cid, exist	:= contaminatoin_targets[seq]
	if exist {
		if(id!=cid){
			log.Fatalf("The target %v has different fasta ids: %s and %s. I am lost.", seq, id, cid)
		}
		rel:=relation{id,parent_id}
		_,linkexist := contamination_links[rel]
		if (linkexist) {return false} //we did nothing, link exist
		_,linkexist = contamination_links[relation{parent_id,id}]
		if (linkexist) {return false} //we did nothing, rev link exist
		contamination_links[rel]=true
		return true
	} else { //add
		contaminatoin_targets[seq]=id
		contamination_links[relation{id,parent_id}]=true
		return true
	}
}

var infile = flag.String("in", "", "file with contaminators")
var threshold = flag.Int("th", math.MaxInt32, "threshold for distance")

func main() {
	mismatch_cost:=1	
	indel_cost:=1
	
	flag.Parse()
	if *infile == "" || *threshold ==  math.MaxInt32 {
		log.Fatal("Please provide file for contaminators and threshold. See nwwreject --help")
	}

	distance_threshold:=*threshold

	time_start := time.Now()
	targets:=read_all_targets_from_fasta("unique_targets_nuc.fa")
	contaminators:=read_all_contaminators(*infile)
	contaminatoin_targets:=make(map[string]string)
	contaminatoin_links:=make(map[relation]bool)
	//key is seq desc.id is fasta id, desc.patrent.ids is why did it
	for nc,contaminator:= range contaminators{
		for target,targetname :=range targets {
			_,ok:=nwwreject.Distance(contaminator,target,mismatch_cost,indel_cost,distance_threshold)
			if ok {
				add_to_contaminatoin_targets(contaminatoin_targets,contaminatoin_links,target,targetname,fmt.Sprint("initial-target ",nc))
			}
		}
	}
	log.Println("contaminatoin_targets: ",len(contaminatoin_targets),"relations: ",len(contaminatoin_links))
	pass:=1
	exhausted:=false
	for !exhausted {
		exhausted=true 
		for seq,cid := range contaminatoin_targets{
			for target,targetname :=range targets {
				if targetname==cid {continue}
				_,ok:=nwwreject.Distance(seq,target,mismatch_cost,indel_cost,distance_threshold)
				if ok {
					if add_to_contaminatoin_targets(contaminatoin_targets,contaminatoin_links,target,targetname,cid) {
						exhausted=false
					}
				}
			}
		}
		if len(contaminatoin_targets)>=len(targets) {exhausted=true}
		log.Println("Pass ",pass, " finished...")
		log.Println("contaminatoin_targets: ",len(contaminatoin_targets),"relations: ",len(contaminatoin_links))
		pass++
	}
	//fmt.Println(contaminatoin_targets)
	//fmt.Println(contaminatoin_links)
	for target, id := range contaminatoin_targets {
		fmt.Println(target,"\t",id)
	}
	log.Println("contaminatoin_targets: ",len(contaminatoin_targets),"relations: ",len(contaminatoin_links)," targets: ", len(targets))
	log.Println(time.Since(time_start))
}
