package main

import (
	"encoding/csv"
//	"fmt"
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

func cdrs3aa(clones []clone) []string {
	cdrs := make([]string,len(clones))
	for n, clone := range clones {cdrs[n]=clone.cdr3aa}
	return cdrs
}

func main() {
	samples_clones := [][]clone {
		readclones("vdj_.S22_clones.txt"),
		readclones("vdj_.S23_clones.txt"),
		readclones("vdj_.S24_clones.txt"),
		readclones("vdj_.S25_clones.txt"),
		readclones("vdj_.S26_clones.txt"),
		readclones("vdj_.S27_clones.txt"),
		readclones("vdj_.S28_clones.txt"),
		readclones("vdj_.S29_clones.txt"),
	}
	//we read everything in the collection
	
	//var common_clones [][] *clone
	cdrs:=cdrs3aa(samples_clones[0])
	println(cdrs[0])
	println(cdrs[1])

}


