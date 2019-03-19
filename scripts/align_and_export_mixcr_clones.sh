#we run it from the folder where the clone files will appear as result
#we suppose that mixcr is at ../tools/mixcr.version

path_to_mixcr="../tools/mixcr-3.0.5"
mixcr="java -Xmx4g -Xms3g -jar ${path_to_mixcr}/mixcr.jar"
#path_to_vdjtools="../vdjtools-1.2.1"
#vdjtools="java -jar ${path_to_vdjtools}/vdjtools-1.2.1.jar"
fastq_path="../fastq"

SAMPLES=(S1 S2 S3 S4 S5 S6 S7 S8 S22 S23 S24 S25 S26 S27 S28 S29)
SAMPLE=(S1) #this is to test
CHAINS=(IGH IGL IGK)
CHAIN=(IGH)
#align -g (--save-reads) is to save reads

#for sample in ${SAMPLES[@]}; do
for sample in ${SAMPLE[@]}; do
	fastq_pair=`find $fastq_path -name "*_${sample}_*"`
	vdjca_file=${sample}.vdjca
	if [ ! -f $vdjca_file ]
	then 
		$mixcr align -s hs -r ${sample}.aln.report.txt ${fastq_pair} $vdjca_file
	else
		echo "File $vdjca_file exists."
	fi
	clns_file=${sample}.clns
	$mixcr assemble -OassemblingFeatures="[CDR1+FR2+CDR2+FR3+CDR3]" -f $vdjca_file $clns_file
	$mixcr exportClones -f $clns_file ${clns_file}.tsv
done
