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
	"sort"
	// "math"
	//	"strings"
)

//output
const write_heavy = true
const support_coverage_for_heavy = 500
const heavy_prefix = "heavy_combined_clones"

//this is for combibnig clones
const max_terminal_del = 1 
const max_mismatches_share = 0.05 

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

//test whether *c1 and *c2 goes to the same combined clone
//the function is suppose to be symmetric, ifeqcl(c1,c2) == ifeqcl(c2,c1)
func ifeqcl(c1 *clone, c2 *clone) bool {
	return string_nonstrict_match(&c1.cdr3aa,&c2.cdr3aa,max_mismatches_share,max_terminal_del) 
}

//shift the start, test the match
func string_nonstrict_match (s1 *string, s2 *string, maxmismatchshare float64, maxtermdel int) bool {
	l1:=len(*s1)
	l2:=len(*s2)
	//too shifted
	if ( l2-l1 > 2*maxtermdel || l1-l2 > 2*maxtermdel) {return false}
	max_mismatch_no:=int(maxmismatchshare*float64(l1))
	//shift s2 to left (e.g. start not from start)
	for shift:=0; shift <= maxtermdel; shift++ {
		shiftedstring:=(*s2)[shift:]
		if (string_nonstrict_match_from_start(s1,&shiftedstring,max_mismatch_no,maxtermdel)) {return true}
	}
	//shift s1 to left (e.g. start not from start)
	for shift:=1; shift <= maxtermdel; shift++ {
		shiftedstring:=(*s1)[shift:]
		if (string_nonstrict_match_from_start(&shiftedstring,s2,max_mismatch_no,maxtermdel)) {return true}
	}
	return false
}

func string_nonstrict_match_from_start (s1 *string, s2 *string, maxmismatch, maxlendiff int) bool {
	//the strings have common start, and they can differ in lenght no more than maxlendiff
	if (len(*s1)>len(*s2)) {
		s3:=s1
		s1=s2
		s2=s3
	}
	//s2 now longer or eq
	if (len(*s2)-len(*s1) > maxlendiff) { return false }
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
	//for alleles
	fmt.Println("v\td\tj")
}


//form a string from alleles map
//map of allele->counters per sample
func string_from_allele_map(allele_map map[string][]int64, sample_names []string) string {
	out:=""
	//all the stuff with alleles is to get alleles in sorted order;
	//allele,counters:=range allele_map give then unsorted
	alleles := make([]string, 0, len(allele_map))
	for allele := range allele_map {
		alleles = append(alleles, allele)
	}
	sort.Strings(alleles)

	for i, allele := range alleles {
		counters:=allele_map[allele]
		//do not print ; before first 
		if i>1 {
			out+="; "
		}
		out=out+allele+": "
		//first_sample is bool to track 
		//whether to print , before the value or not
		first_sample := true
		for n,name:=range sample_names {
			if 0 == counters[n] {continue}
			if !first_sample {
				out+=", "
			}
			first_sample = false
			out+=fmt.Sprint(counters[n]," in ",name)
		}
	}
	return out
}

//print clones info
func print_combined_clone(sample_names []string, combined_clone *list.List) {
	//it is to know the sample number in sample_names by the name
	sample_by_name := make(map[string]int)
	for n, name := range sample_names {
		sample_by_name[name] = n
	}

	//prepare what-to-print
	cdrs := make(map[string]bool)
	counts := make([]int64, len(sample_names))
	freqs := make([]float64, len(sample_names))
	v_alleles := make(map[string][]int64) //map of allele->counters per sample
	d_alleles := make(map[string][]int64) //map of allele->counters per sample
	j_alleles := make(map[string][]int64) //map of allele->counters per sample

	//counts
	//the clone here is the list element with *clone as Value
	for the_clone := combined_clone.Front(); the_clone != nil; the_clone = the_clone.Next() {
		clone_ptr := the_clone.Value.(*clone) //not to convert every step
		cdrs[clone_ptr.cdr3aa] = true
		counts[sample_by_name[clone_ptr.sample]] += clone_ptr.count
		freqs[sample_by_name[clone_ptr.sample]] += clone_ptr.freq
		//init if it is the first mention
		if 0 == len(v_alleles[clone_ptr.v]) { v_alleles[clone_ptr.v]=make([]int64,len(sample_names)) }
		//add to counter
		v_alleles[clone_ptr.v][sample_by_name[clone_ptr.sample]]+=clone_ptr.count
		//init if it is the first mention
		if 0 == len(d_alleles[clone_ptr.d]) { d_alleles[clone_ptr.d]=make([]int64,len(sample_names)) }
		//add to counter
		d_alleles[clone_ptr.d][sample_by_name[clone_ptr.sample]]+=clone_ptr.count
		//init if it is the first mention
		if 0 == len(j_alleles[clone_ptr.j]) { j_alleles[clone_ptr.j]=make([]int64,len(sample_names)) }
		//add to counter
		j_alleles[clone_ptr.j][sample_by_name[clone_ptr.sample]]+=clone_ptr.count
	}

	//print - names
	first_cdr := true
	//fmt.Print(combined_clone.Front().Value.(*clone).cdr3aa, "\t")
	for cdr, _ := range cdrs {
		if !first_cdr {
			fmt.Print(", ")
		}
		first_cdr = false
		fmt.Print(cdr)
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

	fmt.Print(string_from_allele_map(v_alleles,sample_names),"\t")
	fmt.Print(string_from_allele_map(d_alleles,sample_names),"\t")
	fmt.Print(string_from_allele_map(j_alleles,sample_names))

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
	//so, we take *clone one-by-one and present them to all clones already in the combined clones
	clonoteque := list.New()

	for the_clone := all_clones.Front(); the_clone != nil; the_clone = the_clone.Next() {
		var found bool
		for the_combined_clone := clonoteque.Front(); the_combined_clone != nil; the_combined_clone = the_combined_clone.Next() {
			found = false
			for the_inner_clone := the_combined_clone.Value.(*list.List).Front(); the_inner_clone != nil; the_inner_clone = the_inner_clone.Next() {
				//if the_clone is in combined with a clone from the the_combined_clone list.List, the_clone also goes to the_combined_clone
				if ifeqcl(the_inner_clone.Value.(*clone), the_clone.Value.(*clone)) {
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

	if write_heavy {
		heavy_clones_file_name:=fmt.Sprint(heavy_prefix,"_",support_coverage_for_heavy,".tsv")
		outf, err := os.Create(heavy_clones_file_name)
		if err != nil {
        panic(err)
    }
		defer outf.Close()
		old_stdout := os.Stdout
		os.Stdout=outf
		print_clones_header(sample_names)
		for the_combined_clone := clonoteque.Front(); the_combined_clone != nil; the_combined_clone = the_combined_clone.Next() {
			var count int64
			for the_inner_clone := the_combined_clone.Value.(*list.List).Front(); the_inner_clone != nil; the_inner_clone = the_inner_clone.Next() {
				count += the_inner_clone.Value.(*clone).count
			}
			if count >= support_coverage_for_heavy  {
				print_combined_clone(sample_names, the_combined_clone.Value.(*list.List))
			}
		}
		os.Stdout=old_stdout
	}
}
