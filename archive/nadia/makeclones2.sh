#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l hostname=hpc03


#Java 8 is requires for mixcr. It's installed (only?) on the hpc03 node, so we limit our request by "-l hostname=hpc03" to this node

cd /mnt/mapr/group/ms-solid/immuno/fastq/
path="fastq/"
pathout="outfastq/"


for f in *_R1_001.fastq.gz
do
	#remove _R*
	x=$(basename $f)
	echo $x
	y="${x%_R*}"
	#add new ends
	fin1="_R1_001.fastq.gz"
	fin2="_R2_001.fastq.gz"
	#run script
	cd /mnt/mapr/group/ms-solid/immuno/
	f1=$path$y$fin1
	f2=$path$y$fin2
	echo "mixcr-3.0.3/mixcr align -s hs -f $f1 $f2 $y.vdjca"
	mixcr-3.0.3/mixcr align -s hs -f $f1 $f2 $y.vdjca

	#here is version for CDR3 only assembling
	#mixcr-3.0.3/mixcr assemble -f $y.vdjca $y.clns
	#mixcr-3.0.3/mixcr exportClones -c IGH -f $y.clns $y.IGH.txt
	#java -jar vdjtools-1.2.1/vdjtools-1.2.1.jar Convert -S mixcr $y.IGH.txt vdj
	#this one is for CDR1-CDR3 assembling below
	mixcr-3.0.3/mixcr assemble -OassemblingFeatures="[CDR1+FR2+CDR2+FR3+CDR3]" -f $y.vdjca $y.all_cdrs.clns
	mixcr-3.0.3/mixcr exportClones -c IGH -f $y.all_cdrs.clns $y.all_cdrs.IGH.txt
	java -jar vdjtools-1.2.1/vdjtools-1.2.1.jar Convert -S mixcr $y.all_cdrs.IGH.txt vdj


done
