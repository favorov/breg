package main

import (
	"bufio"
	//"encoding/csv"
	//"flag"
	//"fmt"
	"io"
	"log"
	"os"
	//"strconv"
	"time"
)

//being wise, we wold make a package with the string_nonstrict_match dunction, but we are not
//and, it is a bit different
//shift the start, test the match
func string_nonstrict_match(s1 string, s2 string, max_mismatch_no int, maxtermdel int) bool {
	l1 := len(s1)
	l2 := len(s2)
	//too shifted
	if l2-l1 > 2*maxtermdel || l1-l2 > 2*maxtermdel {
		return false
	}
	//shift s2 to left (e.g. start not from start)
	for shift := 0; shift <= maxtermdel; shift++ {
		shiftedstring := s2[shift:]
		if string_nonstrict_match_from_start(s1, shiftedstring, max_mismatch_no, maxtermdel) {
			return true
		}
	}
	//shift s1 to left (e.g. start not from start)
	for shift := 1; shift <= maxtermdel; shift++ {
		shiftedstring := s1[shift:]
		if string_nonstrict_match_from_start(shiftedstring, s2, max_mismatch_no, maxtermdel) {
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

