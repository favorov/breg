package main

import (
	"bufio"
	"container/list"
	"encoding/csv"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"strconv"
)

//sample	subject_id	cell_type	cloneCount	cloneFraction	targetSequences	bestVGene	aaSeqCDR3	nMutations.total	aaMutations.total	mut.n	mut.aa.n
type clone struct {
	sample string
	subject_id string
	cell_type string
	cloneCount int64
	cloneFraction float64 
	targetSequences string
	bestVGene string 
	aaSeqCDR3 string
	nMutations_total string
	aaMutations_total string
	mut_n int
	mut_aa_n int
}

//test whether *c1 and *c2 goes to the same combined clone
//the function is suppose to be symmetric, ifeqcl(c1,c2) == ifeqcl(c2,c1)
func ifeqcl(c1 *clone, c2 *clone, max_mismatches_share float64, max_terminal_del int) bool {
	return string_nonstrict_match(&c1.aaSeqCDR3, &c2.aaSeqCDR3, max_mismatches_share, max_terminal_del)
}

//shift the start, test the match
func string_nonstrict_match(s1 *string, s2 *string, maxmismatchshare float64, maxtermdel int) bool {
	l1 := len(*s1)
	l2 := len(*s2)
	//too shifted
	if l2-l1 > 2*maxtermdel || l1-l2 > 2*maxtermdel {
		return false
	}
	max_mismatch_no := int(maxmismatchshare * float64(l1))
	//shift s2 to left (e.g. start not from start)
	for shift := 0; shift <= maxtermdel; shift++ {
		shiftedstring := (*s2)[shift:]
		if string_nonstrict_match_from_start(s1, &shiftedstring, max_mismatch_no, maxtermdel) {
			return true
		}
	}
	//shift s1 to left (e.g. start not from start)
	for shift := 1; shift <= maxtermdel; shift++ {
		shiftedstring := (*s1)[shift:]
		if string_nonstrict_match_from_start(&shiftedstring, s2, max_mismatch_no, maxtermdel) {
			return true
		}
	}
	return false
}

func string_nonstrict_match_from_start(s1 *string, s2 *string, maxmismatch, maxlendiff int) bool {
	//the strings have common start, and they can differ in lenght no more than maxlendiff
	if len(*s1) > len(*s2) {
		s3 := s1
		s1 = s2
		s2 = s3
	}
	//s2 now longer or eq
	if len(*s2)-len(*s1) > maxlendiff {
		return false
	}
	mismatches := 0
	for i := 0; i < len(*s1); i++ {
		if (*s1)[i] != (*s2)[i] {
			mismatches++
			if mismatches > maxmismatch {
				return false
			}
		}
	}
	return true
}

//read from file : we want a buffer
//read appending to clones; it is a list.List of *clone
func readclones_from_file(filename string, clones *list.List) {
	f, err := os.Open(filename)
	if err != nil {
		log.Fatalf("Error opening file: %v", err)
	}
	defer f.Close()

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

		var record_clone clone
		record_clone.sample = record[0]
		record_clone.subject_id = record[1]
		record_clone.cell_type = record[2]
		record_clone.cloneCount,_ = strconv.ParseInt(record[3], 10, 64)
		record_clone.cloneFraction,_ = strconv.ParseFloat(record[4], 64)
		record_clone.targetSequences = record[5]
		record_clone.bestVGene = record[6]
		record_clone.aaSeqCDR3 = record[7]
		record_clone.nMutations_total= record[8]
		record_clone.aaMutations_total = record[9]
		record_clone.mut_n, _ = strconv.Atoi(record[10])
		record_clone.mut_aa_n, _ = strconv.Atoi(record[11])
		clones.PushBack(&record_clone)
	}
}

func write_list_of_clones_to_fasta (clonelist *list.List,input_file string, output_files_folder string,output_prefix string,output_suffix string) {
		var outfilename string
		front:=clonelist.Front().Value.(*clone)
		outfilename=fmt.Sprint(output_files_folder,"/",output_prefix,front.aaSeqCDR3,"_",front.subject_id,output_suffix)
		f, err := os.Create(outfilename)
    if(err!=nil){log.Fatalf("Error creating file: %v", err)}
		defer f.Close()
		//writing fasta
		var ws string
		ws=fmt.Sprint("; fasta for ",clonelist.Front().Value.(*clone).aaSeqCDR3," clone from ",input_file," file.")
		counter:=0
		for the_inner_clone := clonelist.Front(); the_inner_clone != nil; the_inner_clone = the_inner_clone.Next() {
			counter++
			ws=fmt.Sprint(">",strconv.Itoa(counter),"|",the_inner_clone.Value.(*clone).cell_type,"\n")
			_,err=f.WriteString(ws)
			if(err!=nil){log.Fatalf("Error writing to fasta file file: %v", err)}
			ws=fmt.Sprint(the_inner_clone.Value.(*clone).targetSequences,"\n","\n")
			_,err=f.WriteString(ws)
			if(err!=nil){log.Fatalf("Error writing to fasta file file: %v", err)}
		}
}


func main() {
	var output_files_folder string //where to write
	var output_prefix string //the name of output file prefix 
	var output_suffix string //the name of output file suffix 
	var max_terminal_del int
	var max_mismatches_share float64
	var input_file string
	//What is the input
	flag.StringVar(&output_files_folder, "out", "b-and-b", "folder with the output files")
	flag.StringVar(&output_prefix, "prefix", "vgenes_", "common prefix of all the output files")
	flag.StringVar(&output_suffix, "suffix", ".fasta", "common suffix for putput files")
	//How to combine clones, probably, 1 and 0.05, now we put 0
	flag.IntVar(&max_terminal_del, "terminal-del", 0, "Maximal terminal CDR3 deletion difference that is allowed inside one clone")
	flag.Float64Var(&max_mismatches_share, "mismatches-share", 0, "Maximal share of in-CDR3 mismatches that is allowed inside one clone")

	flag.Usage = func() {
  	fmt.Printf("Usage of %s:\n", os.Args[0])
    fmt.Printf("    %s inputfile ... \n",os.Args[0])
    flag.PrintDefaults()
	}

	flag.Parse()

	//read the filtered clone list from a file, it is be here 
	if flag.NArg() != 1 { 
			flag.Usage()
			os.Exit(11)
	}

	input_file=flag.Arg(0) // the first positional argumant

	_ = os.Mkdir(output_files_folder, 0755)
	
	all_clones := list.New()

	readclones_from_file(input_file, all_clones)


	//organise them to combined_clones list.List<*list.List<*clone>> ; we believe that ifeqcl symmetric
	//so, we take *clone one-by-one and present them to all clones already in the combined clones
	clonoteque := list.New()

	for the_clone := all_clones.Front(); the_clone != nil; the_clone = the_clone.Next() {
		var found bool
		for the_combined_clone := clonoteque.Front(); the_combined_clone != nil; the_combined_clone = the_combined_clone.Next() {
			found = false
			for the_inner_clone := the_combined_clone.Value.(*list.List).Front(); the_inner_clone != nil; the_inner_clone = the_inner_clone.Next() {
				//if the_clone is in combined with a clone from the the_combined_clone list.List, the_clone also goes to the_combined_clone
				if ifeqcl(the_inner_clone.Value.(*clone), the_clone.Value.(*clone), max_mismatches_share, max_terminal_del) {
					the_combined_clone.Value.(*list.List).PushBack(the_clone.Value.(*clone))
					found = true
					break
				}
			}
			//if we are here, and we added the_clone somewhere - break
			if found {
				break
			}
			// else, we go on to the next combined clone
		}
		if found {
			continue
		}
		//if we are here not by break - add new combined clone to the clonoteque and put the_clone there
		clonoteque.PushBack(list.New())
		clonoteque.Back().Value.(*list.List).PushBack(the_clone.Value.(*clone))
	}
	
	for the_combined_clone := clonoteque.Front(); the_combined_clone != nil; the_combined_clone = the_combined_clone.Next() {
		counter:=0
		subject_ids := make(map[string]bool)
		cell_types := make(map[string]bool)
		for the_inner_clone := the_combined_clone.Value.(*list.List).Front(); the_inner_clone != nil; the_inner_clone = the_inner_clone.Next() {
			counter++
			subject_ids[the_inner_clone.Value.(*clone).subject_id]=true
			cell_types[the_inner_clone.Value.(*clone).cell_type]=true
		} //how many subjects do we span? how many cell types do we span?
		if (counter<=2) {
			continue //2-point tree is dull
		}
		if	(len(cell_types)<2) {
			continue
		}
		if	(len(subject_ids)>1) {
			continue
		}
		//prepare fasta....
		write_list_of_clones_to_fasta(the_combined_clone.Value.(*list.List),input_file,output_files_folder,output_prefix,output_suffix)
	}
}
