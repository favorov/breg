package main

import (
	"encoding/csv"
	"time"
)



//we use https://gist.github.com/quwubin/fdf9a9b40f4c4fbbeb02
import (
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


//being wise, we wold make a package with the string_nonstrict_match dunction, but we are not
//and, it is a bit different
//shift the start, test the match
func string_nonstrict_match(s1 string, s2 string, max_mismatch_no int, max_termdel_per_end int) bool {
	l1 := len(s1)
	l2 := len(s2)
	//too shifted
	if l2-l1 > 2*max_termdel_per_end || l1-l2 > 2*max_termdel_per_end {
		return false
	}
	//shift s2 to left (e.g. start not from start)
	for shift := 0; shift <= max_termdel_per_end; shift++ {
		shiftedstring := s2[shift:]
		if string_nonstrict_match_from_start(s1, shiftedstring, max_mismatch_no, max_termdel_per_end) {
			return true
		}
	}
	//shift s1 to left (e.g. start not from start)
	for shift := 1; shift <= max_termdel_per_end; shift++ {
		shiftedstring := s1[shift:]
		if string_nonstrict_match_from_start(shiftedstring, s2, max_mismatch_no, max_termdel_per_end) {
			return true
		}
	}
	return false
}

func string_nonstrict_match_from_start(s1 string, s2 string, maxmismatch, maxlendiff int) bool {
	//the strings have common start, and they can differ in lenght no more than maxlendiff
	if len(s1) > len(s2) {
		s3 := s1
		s1 = s2
		s2 = s3
	}
	//s2 now longer or eq
	if len(s2)-len(s1) > maxlendiff {
		return false
	}
	mismatches := 0
	for i := 0; i < len(s1); i++ {
		if s1[i] != s2[i] {
			mismatches++
			if mismatches > maxmismatch {
				return false
			}
		}
	}
	return true
}


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

type cont_desc struct {
  id string
  parent_ids []string
}

//return is did we add new or just add a 
func add_to_contaminatoin_targets (contaminatoin_targets map[string]cont_desc, seq, id, parent_id string) bool {
	cdesc, exist	:= contaminatoin_targets[seq]
	if exist {
		if(id!=cdesc.id){
			log.Fatalf("The target %v has different fasta ids: %s and %s. I am lost.", seq, id, cdesc.id)
		}
		cdesc.parent_ids=append(cdesc.parent_ids,parent_id)
		return false
	} else {
		var desc cont_desc
		desc.id=id
		desc.parent_ids=append(make([]string,0),parent_id)
		contaminatoin_targets[seq]=desc
		return true
	}
}

func main() {
	mutations_per_step:=2
	term_indel_per_step_per_end:=1

	time_start := time.Now()
	targets:=read_all_targets_from_fasta("unique_targets_nuc.fa")
	contaminators:=read_all_contaminators("contaminatoin_targets.txt")
	contaminatoin_targets:=make(map[string]cont_desc)
	//key is seq desc.id is fasta id, desc.patrent.ids is why did it
	added:=0
	for nc,contaminator:= range contaminators{
		for target,targetname :=range targets {
			if string_nonstrict_match(contaminator,target,mutations_per_step,term_indel_per_step_per_end) {
				if add_to_contaminatoin_targets(contaminatoin_targets,target,targetname,fmt.Sprint("initial-target ",nc)) {added++}
			}
		}
	}
	log.Println("contaminators: ",len(contaminators)," added: ",added," contaminatoin_targets: ",len(contaminatoin_targets))
	pass:=1
	exhausted:=false
	for !exhausted {
		exhausted=true 
		for seq,cont_desc := range contaminatoin_targets{
			for target,targetname :=range targets {
				if string_nonstrict_match(seq,target,mutations_per_step,term_indel_per_step_per_end) {
					if add_to_contaminatoin_targets(contaminatoin_targets,target,targetname,cont_desc.id) {added++}
					exhausted=false //we have what to add
				}
			}
		}
		if len(contaminatoin_targets)>=len(targets) {exhausted=true}
		log.Println("Pass ",pass, " finished...")
		log.Println(" added: ",added," contaminatoin_targets: ",len(contaminatoin_targets))
		pass++
	}
	fmt.Println(contaminatoin_targets)
	log.Println(time.Since(time_start))
}
