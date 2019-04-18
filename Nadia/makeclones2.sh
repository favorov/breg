#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -l hostname=hpc03


#Java 8 is requires for mixcr. It's installed (only?) on the hpc03 node, so we limit our request by "-l hostname=hpc03" to this node

cd /mnt/mapr/group/ms-solid/immuno/fastq/
path="fastq/"
pathout="mixcr_results/"


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
	y=$pathout$y
	echo "mixcr-3.0.3/mixcr align -s hs -f $f1 $f2 $y.vdjca"
	rm $y.align.report
	mixcr-3.0.3/mixcr align -s hs -f -r $y.align.report $f1 $f2 $y.vdjca

	#here is version for CDR3 only assembling
	rm $y.assemble.report
	mixcr-3.0.3/mixcr assemble -f -OaddReadsCountOnClustering=true -OseparateByV=true -r $y.assemble.report $y.vdjca $y.clns
	mixcr-3.0.3/mixcr exportClones -c IGH -f -cloneId -count -fraction -targetSequences -vGene -aaFeature CDR3 $y.clns $y.IGH.txt
	mixcr-3.0.3/mixcr exportClones -c IGL -f -cloneId -count -fraction -targetSequences -vGene -aaFeature CDR3 $y.clns $y.IGL.txt
	mixcr-3.0.3/mixcr exportClones -c IGK -f -cloneId -count -fraction -targetSequences -vGene -aaFeature CDR3 $y.clns $y.IGK.txt
	#java -jar vdjtools-1.2.1/vdjtools-1.2.1.jar Convert -S mixcr $y.IGH.txt vdj
	#this one is for CDR1-CDR3 assembling below
	rm $y.all_cdrs.assemble.report
	mixcr-3.0.3/mixcr assemble -f -OaddReadsCountOnClustering=true -OassemblingFeatures="[CDR1+FR2+CDR2+FR3+CDR3]" -r $y.all_cdrs.assemble.report $y.vdjca $y.all_cdrs.clns
	mixcr-3.0.3/mixcr exportClones -c IGH -f -cloneId -count -fraction -targetSequences -vGene -aaFeature CDR3 -nMutations CDR1 -nMutations FR2 -nMutations CDR2 -nMutations FR3 -aaMutations CDR1 -aaMutations FR2 -aaMutations CDR2 -aaMutations FR3  $y.all_cdrs.clns $y.all_cdrs.nMut.IGH.txt
	mixcr-3.0.3/mixcr exportClones -c IGL -f -cloneId -count -fraction -targetSequences -vGene -aaFeature CDR3 -nMutations CDR1 -nMutations FR2 -nMutations CDR2 -nMutations FR3 -aaMutations CDR1 -aaMutations FR2 -aaMutations CDR2 -aaMutations FR3  $y.all_cdrs.clns $y.all_cdrs.nMut.IGL.txt
	mixcr-3.0.3/mixcr exportClones -c IGK -f -cloneId -count -fraction -targetSequences -vGene -aaFeature CDR3 -nMutations CDR1 -nMutations FR2 -nMutations CDR2 -nMutations FR3 -aaMutations CDR1 -aaMutations FR2 -aaMutations CDR2 -aaMutations FR3  $y.all_cdrs.clns $y.all_cdrs.nMut.IGK.txt
	#java -jar vdjtools-1.2.1/vdjtools-1.2.1.jar Convert -S mixcr $y.all_cdrs.nMut.IGH.txt vdj


done
