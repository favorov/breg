#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l hostname=hpc03


#Java 8 is requires for mixcr. It's installed (only?) on the hpc03 node, so we limit our request by "-l hostname=hpc03" to this node

#we start in immunoseq home

path="fastq/"
pathout="mixcr_results/"
mixcr="tools/mixcr-3.0.5/mixcr"

pushd $path
fastqs_001=*_R1_001.fastq.gz
popd

for f in $fastqs_001
do
	#remove _R*
	x=$(basename $f)
	echo $x
	y="${x%_R*}"
	#add new ends
	fin1="_R1_001.fastq.gz"
	fin2="_R2_001.fastq.gz"
	#run script
	f1=$path$y$fin1
	f2=$path$y$fin2
	y=$pathout$y
	echo "$mixcr align -s hs -f $f1 $f2 $y.vdjca"
	rm $y.align.report
	$mixcr align -s hs -f -r $y.align.report $f1 $f2 $y.vdjca

	#here is version for CDR3 only assembling
	rm $y.assemble.report
	$mixcr assemble -f -OaddReadsCountOnClustering=true -OseparateByV=true -r $y.assemble.report $y.vdjca $y.clns
	$mixcr exportClones -c IGH -f -cloneId -count -fraction -targetSequences -vGene -aaFeature CDR3 $y.clns $y.IGH.txt
	$mixcr exportClones -c IGL -f -cloneId -count -fraction -targetSequences -vGene -aaFeature CDR3 $y.clns $y.IGL.txt
	$mixcr exportClones -c IGK -f -cloneId -count -fraction -targetSequences -vGene -aaFeature CDR3 $y.clns $y.IGK.txt
	#java -jar vdjtools-1.2.1/vdjtools-1.2.1.jar Convert -S mixcr $y.IGH.txt vdj
	#this one is for CDR1-CDR3 assembling below
	rm $y.all_cdrs.assemble.report
	$mixcr assemble -f -OaddReadsCountOnClustering=true -OassemblingFeatures="[CDR1+FR2+CDR2+FR3+CDR3]" -r $y.all_cdrs.assemble.report $y.vdjca $y.all_cdrs.clns
	$mixcr exportClones -c IGH -f -cloneId -count -fraction -targetSequences -vGene -aaFeature CDR3 -nMutations CDR1 -nMutations FR2 -nMutations CDR2 -nMutations FR3 -aaMutations CDR1 -aaMutations FR2 -aaMutations CDR2 -aaMutations FR3  $y.all_cdrs.clns $y.all_cdrs.nMut.IGH.txt
	$mixcr exportClones -c IGL -f -cloneId -count -fraction -targetSequences -vGene -aaFeature CDR3 -nMutations CDR1 -nMutations FR2 -nMutations CDR2 -nMutations FR3 -aaMutations CDR1 -aaMutations FR2 -aaMutations CDR2 -aaMutations FR3  $y.all_cdrs.clns $y.all_cdrs.nMut.IGL.txt
	$mixcr exportClones -c IGK -f -cloneId -count -fraction -targetSequences -vGene -aaFeature CDR3 -nMutations CDR1 -nMutations FR2 -nMutations CDR2 -nMutations FR3 -aaMutations CDR1 -aaMutations FR2 -aaMutations CDR2 -aaMutations FR3  $y.all_cdrs.clns $y.all_cdrs.nMut.IGK.txt
	#java -jar vdjtools-1.2.1/vdjtools-1.2.1.jar Convert -S mixcr $y.all_cdrs.nMut.IGH.txt vdj


done
