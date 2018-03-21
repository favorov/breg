package main

import (
	"encoding/csv"
//	"fmt"
	"io"
	"bufio"
	"log"
	"os"
	"strconv"
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
	freq float32
//actual order in file:
//count	freq	cdr3nt	cdr3aa	v	d	j	VEnd	DStart	DEnd	JStart
}

func readclones(filename string) []clone {
 	f, err := os.Open(filename)
	if err != nil {
			log.Fatalf("Error opening file: %v", err)
	}
	defer f.Close()
	
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
		count, _ := strconv.ParseInt(record[0],10,64)
		freq, _ := strconv.ParseFloat(record[1],32)
		VEnd, _ := strconv.ParseInt(record[7],10,32)
		DStart, _ := strconv.ParseInt(record[8],10,32)
		DEnd, _ := strconv.ParseInt(record[9],10,32)
		JStart, _ := strconv.ParseInt(record[10],10,32)

		clones = append(clones, clone{
					count: count,
					freq: float32(freq),
					cdr3nt: record[2],
					cdr3aa: record[3],
					v: record[4],
					d: record[5],
					j: record[6],
					VEnd: int(VEnd),
					DStart: int(DStart),
					DEnd: int(DEnd),
					JStart: int(JStart),
    })
	}
	return clones
}

func main() {
	S22 := readclones("vdj_.S22_clones.txt")
	S22 := readclones("vdj_.S22_clones.txt")
	S22 := readclones("vdj_.S22_clones.txt")
}


