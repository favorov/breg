package main

import (
	"bufio"
	"container/list"
	"encoding/csv"
	"fmt"
	"io"
	"log"
	"os"
	"regexp"
	"strconv"
	//	"strings"
)

const support_coverage = 500
const max_terminal_del = 2
const max_mismatches_no = 2

type clone struct {
	cdr3aa string
	cdr3nt string
	v      string
	d      string
	j      string
	VEnd   int
	DStart int
	DEnd   int
	JStart int
	count  int64
	freq   float64
	//actual order in file:
	//count	freq	cdr3nt	cdr3aa	v	d	j	VEnd	DStart	DEnd	JStart
	sample string
	//to know where it comes from
}

//test whether *c1 and *c2 goes to the same common clone
//the function is suppose to be symmetric, ifeqcl(c1,c2) == ifeqcl(c2,c1)
func ifeqcl(c1 *clone, c2 *clone) bool {
	return string_nonstrict_match(&c1.cdr3aa,&c2.cdr3aa,max_mismatches_no,max_terminal_del) 
}

func string_nonstrict_match (s1 *string, s2 *string, maxmismatch, maxtermdel int) bool {
	if (len(*s1)!=len(*s2)) {
		return false
	}
	return string_eq_len_nonstrict_match(s1,s2,maxtermdel) 
}

func string_eq_len_nonstrict_match (s1 *string, s2 *string, maxmismatch int) bool {
	if (len(*s1)!=len(*s2)) {
		log.Println("string_eq_len_nonstrict_match get unequal dting length\n")
		return false
	}
	mismatches := 0
	for i := 0; i < len(*s1); i++ {	
		if ((*s1)[i] != (*s2)[i]) {
			mismatches++
			if (mismatches>maxmismatch) {
				return false
			}
		}
	}
	return true
}

//read from file, sample name is sample_name,
//read appending to clones; it is a list.List of *clone
func readclones_from_file(filename string, sample_name string, clones *list.List) {
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
		record_clone.count, _ = strconv.ParseInt(record[0], 10, 64)
		record_clone.freq, _ = strconv.ParseFloat(record[1], 32)
		record_clone.cdr3nt = record[2]
		record_clone.cdr3aa = record[3]
		record_clone.v = record[4]
		record_clone.d = record[5]
		record_clone.j = record[6]
		record_clone.VEnd, _ = strconv.Atoi(record[7])
		record_clone.DStart, _ = strconv.Atoi(record[8])
		record_clone.DEnd, _ = strconv.Atoi(record[9])
		record_clone.JStart, _ = strconv.Atoi(record[10])
		record_clone.sample = sample_name
		clones.PushBack(&record_clone)
	}
}

func print_clones_header(sample_names []string) {
	fmt.Print("cdr3aa\t")
	//for counts
	for _, sample := range sample_names {
		fmt.Print(sample, "\t")
	}
	//for freqs
	for _, sample := range sample_names {
		fmt.Print(sample, "\t")
	}
	fmt.Println()
}

//print clones info; first, cdr is printed as key
func print_common_clone(sample_names []string, common_clone *list.List) {
	//it is to know the sample number in sample_names by the name
	sample_by_name := make(map[string]int)
	for n, name := range sample_names {
		sample_by_name[name] = n
	}

	//prepare what-to-print
	cdrs := make(map[string]bool)
	counts := make([]int64, len(sample_names))
	freqs := make([]float64, len(sample_names))

	//counts
	//the clone here is the list element with *clone as Value
	for the_clone := common_clone.Front(); the_clone != nil; the_clone = the_clone.Next() {
		clone_ptr := the_clone.Value.(*clone) //not to convert every step
		cdrs[clone_ptr.cdr3aa] = true
		counts[sample_by_name[clone_ptr.sample]] += clone_ptr.count
		freqs[sample_by_name[clone_ptr.sample]] += clone_ptr.freq
	}

	//print - names
	firstcdr := true
	//fmt.Print(common_clone.Front().Value.(*clone).cdr3aa, "\t")
	for cdr, _ := range cdrs {
		if !firstcdr {
			fmt.Print(", ")
		}
		fmt.Print(cdr)
		firstcdr = false
	}
	fmt.Print("\t")
	//print - counts
	for n, _ := range sample_names {
		fmt.Print(counts[n], "\t")
	}
	//print - freqs
	for n, _ := range sample_names {
		fmt.Printf("%6f\t", freqs[n])
	}

	fmt.Println()
	return
}

func main() {
	sample_files := []string{
		"vdj_.S22_clones.txt",
		"vdj_.S23_clones.txt",
		"vdj_.S24_clones.txt",
		"vdj_.S25_clones.txt",
		"vdj_.S26_clones.txt",
		"vdj_.S27_clones.txt",
		"vdj_.S28_clones.txt",
		"vdj_.S29_clones.txt",
	}

	//read all clone from file to clone table [samples][lines]
	all_clones := list.New()
	var sample_names []string

	//prepare to parse "S22", "S23", etc in strings
	sample_name_regstr := "S[0-9]*"
	sample_regexp, reerr := regexp.Compile(sample_name_regstr)
	if reerr != nil {
		log.Fatalf("Error compiling regexp %v: %v", sample_name_regstr, reerr)
	}

	//read clones from files
	for _, sample_file := range sample_files {
		sample_name := sample_regexp.FindString(sample_file)
		sample_names = append(sample_names, sample_name)
		readclones_from_file(sample_file, sample_name, all_clones)
	}

	//organise them to combined_clones list.List<*list.List<*clone>> ; we believe that ifeqcl symmetric
	//so, we take *clone one-by-one and present them to all clones already in the common clones
	clonoteque := list.New()

	for the_clone := all_clones.Front(); the_clone != nil; the_clone = the_clone.Next() {
		var found bool
		for the_common_clone := clonoteque.Front(); the_common_clone != nil; the_common_clone = the_common_clone.Next() {
			found = false
			for the_inner_clone := the_common_clone.Value.(*list.List).Front(); the_inner_clone != nil; the_inner_clone = the_inner_clone.Next() {
				//if the_clone is in common with a clone from the the_common_clone list.List, the_clone also goes to the_common_clone
				if ifeqcl(the_inner_clone.Value.(*clone), the_clone.Value.(*clone)) {
					the_common_clone.Value.(*list.List).PushBack(the_clone.Value.(*clone))
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

	print_clones_header(sample_names)
	for the_common_clone := clonoteque.Front(); the_common_clone != nil; the_common_clone = the_common_clone.Next() {
		var count int64
		for the_inner_clone := the_common_clone.Value.(*list.List).Front(); the_inner_clone != nil; the_inner_clone = the_inner_clone.Next() {
			count += the_inner_clone.Value.(*clone).count
		}
		if count >= support_coverage  {
			print_common_clone(sample_names, the_common_clone.Value.(*list.List))
		}
	}
}
