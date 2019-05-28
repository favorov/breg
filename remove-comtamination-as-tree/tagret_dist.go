package main

import (
	"github.com/favorov/nwwreject"
	//"./nwwreject"
	"time"
)
import "runtime/pprof"


//we use https://gist.github.com/quwubin/fdf9a9b40f4c4fbbeb02
import (
	"flag"
	"math"
	"bytes"
	//"fmt"
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

func parseFastaRecords(fastaFh io.Reader) chan fasta {
//renamed parse to parseFastaRecords

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


func read_all_targets_from_fasta(fn string) (fasta map[string]string, maxlen int) {
	// Open the fastq file specified on the command line
	// for reading:
	fh, err := os.Open(fn)
	// Check for open errors and abort:
	if err != nil {
		log.Fatalf("Error opening file: %v", err)
	}
	defer fh.Close()
	fasta= make(map[string]string)
	maxlen=0
	for record := range parseFastaRecords(fh) {
		_,was:=fasta[record.id]
		if was { log.Fatalf("Nonunique Fasta record. %s",record.id)}
		fasta[record.id]=record.seq
		if maxlen<len(record.seq) {maxlen=len(record.seq)}
	}
	return
}

var threshold = flag.Int("th", math.MaxInt32, "threshold for distance")

func main() {
//profile
	if true {
		f, err := os.Create("pofile.target_dist.prof")
		if err != nil {
				log.Fatal("could not create CPU profile: ", err)
		}
		defer f.Close()
		if err := pprof.StartCPUProfile(f); err != nil {
				log.Fatal("could not start CPU profile: ", err)
		}
		defer pprof.StopCPUProfile()
	}
//remainder
	//log.Println("Using nwwreject version ",nwwreject.Version)
	mismatch_cost:=1
	indel_cost:=1

	flag.Parse()
	if  *threshold ==  math.MaxInt32 {
		log.Fatal("Please provide file for contaminators and threshold. See nwwreject --help")
	}

	distance_threshold:=*threshold

	time_start := time.Now()

	target_map,max_target_len:=read_all_targets_from_fasta("unique_targets_nuc.fa")

	nwwreject.Init_distance_matrix(max_target_len)

	count:=0
	succount:=0
	for id_1,t1 := range target_map{
		for id_2,t2 :=range target_map {
			if id_1==id_2 {continue}
			//dist,ok:=nwwreject.Distance(t1,t2,mismatch_cost,indel_cost,distance_threshold)
			_,ok:=nwwreject.Distance(t1,t2,mismatch_cost,indel_cost,distance_threshold)
			if(ok) {succount++}
		}
		count++
		log.Println(count," of ",len(target_map),"  ",succount," of ",count*len(target_map))
		if count*100>len(target_map) {break}
	}
	log.Println(time.Since(time_start))
}
