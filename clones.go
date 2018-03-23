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

func readclones(filename string) []clone {
 	f, err := os.Open(filename)
	if err != nil {
			log.Fatalf("Error opening file: %v", err)
	}
	defer f.Close()

	regstr:="S[0-9]*"
	sampleregexp,reerr := regexp.Compile(regstr)
	if reerr != nil {
			log.Fatalf("Error compiling regexp %v: %v",regstr,reerr)
	}
	sample:=sampleregexp.FindString(filename)

	var clones []clone

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
		record_clone.sample = sample
		clones = append(clones, record_clone)
	}
	return clones
}

//extract cdr3aa from each clone in clones []clone
//make a []srtring
func cdrs3aa(clones []clone) []string {
	cdrs := make([]string,len(clones))
	for n, clone := range clones {cdrs[n]=clone.cdr3aa}
	return cdrs
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
func print_clones_map_pair (key string,sample_names []string,clones []*clone){
	//counts
	counts:=make([]int64,len(sample_names)) 
	for n, sample := range sample_names {
		for _, clone := range clones {
			if (clone.sample == sample) {
				counts[n]=counts[n]+clone.count //sum is for future
			}
		}
	}
	var count_counts int
	//how mane samples gained to the clone more that one read?
	for _, count := range counts {
		if (count>1) {
			count_counts=count_counts+1
		}
	}
	//id 1, do not print
	if (count_counts<2) {return}

	print(key,"\t")
	for _,count := range counts {
		print(count,"\t")
	}
	//freqs
	for _, sample := range sample_names {
		var freq float64
		for _, clone := range clones {
			if (clone.sample == sample) {
				freq=freq+clone.freq //sum is for future
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
	var samples_clones [][]clone 
	for _, sample_file := range sample_files {
		samples_clones = append(samples_clones,readclones(sample_file))
	}

	var sample_names []string
	map_of_clones := make(map[string][]*clone)
	
	for _, sample_clones := range samples_clones {
		//get sample name
		sample_names=append(sample_names,sample_clones[0].sample) 
		//and code everything form sample into the map_of_clones
		cdrs := cdrs3aa(sample_clones)
		for n,cdr := range cdrs{
				map_of_clones[cdr]=append(map_of_clones[cdr],&sample_clones[n])
				//looks simple... but if there is no cdr key, map_of_clones[cdr] return zero []clones, so we append and thus init
		}
	}
	
	map_of_common_clones := make(map[string][]*clone)
	for cdr, clones :=range map_of_clones {
		if (len(clones) > 1) {
			map_of_common_clones[cdr] = clones 
		}
	}

	print_clones_header(sample_names)
	for cdr, clones :=range map_of_common_clones {
		print_clones_map_pair(cdr,sample_names,clones)
	}
}


