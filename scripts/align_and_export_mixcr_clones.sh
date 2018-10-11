#! /bin/bash
#we run it from the folder where the clone files will appear as result

path_to_mixcr="/home/sasha/gabibov/tools/mixcr-2.1.12"
path_to_vdjtools="/home/sasha/gabibov/tools/vdjtools-1.1.10"
mixcr="java -jar ${path_to_mixcr}/mixcr.jar"
vdjtools="java -jar ${path_to_vdjtools}/vdjtools-1.1.10.jar"
fastq_path="../fastq"

SAMPLES=(S1 S2 S3 S4 S5 S6 S7 S8 S22 S23 S24 S25 S26 S27 S28 S29)
SAMPLE=(S1)
CHAINS=(IGH IGL IGK)
CHAIN=(IGH)
#align -g (--save-reads) is to save reads

for sample in ${SAMPLES[@]}; do
	for chain in ${CHAINS[@]}; do
		fastq_pair=`find $fastq_path -name "*_${sample}_*"`
		vdjca_file=${sample}.${chain}.vdjca
		if [ ! -f $vdjca_file ]
		then 
			echo $mixcr align -c ${chain} -r ${sample}.${chain}.aln.report.txt ${fastq_pair} $vdjca_file
			$mixcr align -c ${chain} -r ${sample}.${chain}.aln.report.txt ${fastq_pair} $vdjca_file
		else
			echo "File $vdjca_file exists."
		fi
		
		clns_file=${sample}.${chain}.clns
		if [ ! -f $clns_file ]
		then 
			echo $mixcr assemble -r ${sample}.${chain}.assemb.report.txt $vdjca_file $clns_file
			$mixcr assemble -r ${sample}.${chain}.assemb.report.txt $vdjca_file $clns_file
		else
			echo "File $clns_file exists."
		fi
		clones_text_file=${sample}.${chain}.clones.txt
		if [ ! -f $clones_text_file ]
		then 
			echo $mixcr exportClones $clns_file $clones_text_file
			$mixcr exportClones $clns_file $clones_text_file
		else
			echo "File $clones_text_file exists."
		fi
		vdj_clones_text_file=vdj.${sample}.${chain}.clones.txt
		if [ ! -f $vdj_clones_text_file ]
		then 
			echo $vdjtools Convert -S mixcr $clones_text_file vdj
			$vdjtools Convert -S mixcr $clones_text_file vdj
		else
			echo "File $vdj_clones_text_file exists."
		fi
	done
done
