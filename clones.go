package main

import (
	"encoding/csv"
//	"fmt"
	"io"
	"log"
	"os"
	"strconv"
)

func main() {
	file := "vdj_.S22_clones.txt"
 	f, err := os.Open(file)
	if err != nil {
			log.Fatalf("Error opening file: %v", err)
	}
	defer f.Close()

	rdr := csv.NewReader(f)
	rdr.Comma = '\t'
	for {
		record, err := rdr.Read()
		if err != nil {
    	if err == io.EOF {
				break
			}
			log.Fatal(err)
		}
		print(record[2],"\n")
	}
}


