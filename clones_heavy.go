package main

import (
	"encoding/csv"
	"fmt"
	"io"
	"bufio"
	"log"
	"os"
	"strconv"
	"regexp"
)

type clone struct {
	cdr3aa string	
	cdr3nt string
	v	string
	d	string
	j	string
	VEnd int
	DStart int
	DEnd int
	JStart int
	count int64
	freq float64
//actual order in file:
//count	freq	cdr3nt	cdr3aa	v	d	j	VEnd	DStart	DEnd	JStart
	sample string
//to know where it comes from
}

type common_clone struct {
	cdr3aa_set [] string
	clones []*clone
	count int64
}

//func are_clones_common (c1 *clone, c2 *clone) bool {
//it will be more complicated later
//	eq := false
//	eq = (c1.cdr3aa == c2.cdr3aa)
//	return eq
//}

//func are_clones_common (cs []*clone, c2 *clone) bool {
//	for _,c1 := range cs {
//		if (are_clones_common(c1,c2)) {return true}
//	}
//	return false
//}

//read from file, sample name is sample_name,
//read appenfding to clones
func readclones_from_file(filename string, sample_name string, clones *[]clone) {
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
		record_clone.count, _ = strconv.ParseInt(record[0],10,64)
		record_clone.freq, _ = strconv.ParseFloat(record[1],32)
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
		*clones = append(*clones, record_clone)
	}
}

//extract cdr keys from map [cdr][]*clone
func keys(map_of_clones map[string][]*clone) []string {
	keys := make([]string, len(map_of_clones))
	i:=0
	for key,_ := range map_of_clones {
    keys[i] = key
		i++
	}
	return keys
}

func print_clones_header(sample_names []string) {
	print("cdr3aa\t")
	//for counts
	for _,sample:= range sample_names {
		print(sample,"\t")
	}
	//for freqs 
	for _,sample:= range sample_names {
		print(sample,"\t")
	}
	print("\n")
}


//print clones info; first, cdr is printed as key
func print_common_clone (sample_names []string,cclone common_clone){
	
	print(cclone.cdr3aa_set[0],"\t")
	//counts
	for _, sample := range sample_names {
		var count int64
		for _, clone := range cclone.clones {
			if (clone.sample == sample) {
				count += clone.count //sum is for future
			}
		print(count,"\t")
		}
	}

	//freqs
	for _, sample := range sample_names {
		var freq float64
		for _, clone := range cclone.clones {
			if (clone.sample == sample) {
				freq += clone.freq //sum is for future
			}
		}
		print(fmt.Sprintf("%6f\t",freq))
	}
	print("\n")
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
	var all_clones []clone
	var sample_names []string

	//prepare to parse "S22", "S23", etc in strings
	sample_name_regstr:="S[0-9]*"
	sample_regexp,reerr := regexp.Compile(sample_name_regstr)
	if reerr != nil {
			log.Fatalf("Error compiling regexp %v: %v",sample_name_regstr,reerr)
	}
	
	//read clones from files
	for _, sample_file := range sample_files {
		sample_name:=sample_regexp.FindString(sample_file)
		sample_names=append(sample_names,sample_name) 
		readclones_from_file(sample_file,sample_name,&all_clones)
	}

	//organise their &  to map
	map_of_clones := make(map[string][]*clone)
	for n, the_clone := range all_clones {
		//and refer all the cones from the sample into the map_of_clones
		map_of_clones[the_clone.cdr3aa]=append(map_of_clones[the_clone.cdr3aa],&(all_clones[n]))
		//looks simple... but if there is no cdr key in the map, map_of_clones[cdr] return zero []*clone, so we append and thus init
	}
	
	var common_clones []common_clone
	
	for cdr, clones :=range map_of_clones {
		var counter int64
		for _, clone := range clones {counter+=clone.count}
		if (counter >= 500) {
			var new_cclone common_clone
			new_cclone.cdr3aa_set = append(new_cclone.cdr3aa_set,cdr)
			//actually, it is just init, cclone.cdr3aa_set is empty 
			new_cclone.clones = clones
			new_cclone.count=counter
			common_clones = append (common_clones,new_cclone)
		}
	}

	print_clones_header(sample_names)
	for _, cclone :=range common_clones {
		print_common_clone(sample_names,cclone)
	}
}


