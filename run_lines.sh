#mixcr align -c IGH -r S22.aln_report.txt Lom1_S22_L001_R1_001.fastq.gz Lom1_S22_L001_R2_001.fastq.gz S22.vdjca
#mixcr align -c IGH -r S23.aln_report.txt Lom2_S23_L001_R1_001.fastq.gz Lom2_S23_L001_R2_001.fastq.gz S23.vdjca

#mixcr align -c IGH -r S24.aln_report.txt Lom3_S24_L001_R1_001.fastq.gz Lom3_S24_L001_R2_001.fastq.gz S24.vdjca
#mixcr align -c IGH -r S25.aln_report.txt Lom4_S25_L001_R1_001.fastq.gz Lom4_S25_L001_R2_001.fastq.gz S25.vdjca
#mixcr align -c IGH -r S26.aln_report.txt Lom5_S26_L001_R1_001.fastq.gz Lom5_S26_L001_R2_001.fastq.gz S26.vdjca
#mixcr align -c IGH -r S27.aln_report.txt Lom6_S27_L001_R1_001.fastq.gz Lom6_S27_L001_R2_001.fastq.gz S27.vdjca
#mixcr align -c IGH -r S28.aln_report.txt Lom7_S28_L001_R1_001.fastq.gz Lom7_S28_L001_R2_001.fastq.gz S28.vdjca
#mixcr align -c IGH -r S29.aln_report.txt Lom8_S29_L001_R1_001.fastq.gz Lom8_S29_L001_R2_001.fastq.gz S29.vdjca


#mixcr assemble -f -r S22_assembly_report.txt S22.vdjca S22.clns
#mixcr assemble -f -r S23_assembly_report.txt S23.vdjca S23.clns
#mixcr assemble -f -r S24_assembly_report.txt S24.vdjca S24.clns
#mixcr assemble -f -r S25_assembly_report.txt S25.vdjca S25.clns
#mixcr assemble -f -r S26_assembly_report.txt S26.vdjca S26.clns
#mixcr assemble -f -r S27_assembly_report.txt S27.vdjca S27.clns
#mixcr assemble -f -r S28_assembly_report.txt S28.vdjca S28.clns
#mixcr assemble -f -r S29_assembly_report.txt S29.vdjca S29.clns

#mixcr exportClones -f S22.clns S22_clones.txt
#mixcr exportClones -f S23.clns S23_clones.txt
#mixcr exportClones -f S24.clns S24_clones.txt
#mixcr exportClones -f S25.clns S25_clones.txt
#mixcr exportClones -f S26.clns S26_clones.txt
#mixcr exportClones -f S27.clns S27_clones.txt
#mixcr exportClones -f S28.clns S28_clones.txt
#mixcr exportClones -f S29.clns S29_clones.txt

#java -jar vdjtools-1.1.7/vdjtools-1.1.7.jar Convert -S mixcr S22_clones.txt vdj_
#java -jar vdjtools-1.1.7/vdjtools-1.1.7.jar Convert -S mixcr S23_clones.txt vdj_
#java -jar vdjtools-1.1.7/vdjtools-1.1.7.jar Convert -S mixcr S24_clones.txt vdj_
#java -jar vdjtools-1.1.7/vdjtools-1.1.7.jar Convert -S mixcr S25_clones.txt vdj_
#java -jar vdjtools-1.1.7/vdjtools-1.1.7.jar Convert -S mixcr S26_clones.txt vdj_
#java -jar vdjtools-1.1.7/vdjtools-1.1.7.jar Convert -S mixcr S27_clones.txt vdj_
#java -jar vdjtools-1.1.7/vdjtools-1.1.7.jar Convert -S mixcr S28_clones.txt vdj_
#java -jar vdjtools-1.1.7/vdjtools-1.1.7.jar Convert -S mixcr S29_clones.txt vdj_

#java -jar vdjtools-1.1.7/vdjtools-1.1.7.jar  OverlapPair -p vdj_.S22_clones.txt vdj_.S23_clones.txt S22_S23
#java -jar vdjtools-1.1.7/vdjtools-1.1.7.jar  OverlapPair -p vdj_.S24_clones.txt vdj_.S25_clones.txt S24_S25
#java -jar vdjtools-1.1.7/vdjtools-1.1.7.jar  OverlapPair -p vdj_.S22_clones.txt vdj_.S24_clones.txt S22_S24
#java -jar vdjtools-1.1.7/vdjtools-1.1.7.jar  OverlapPair -p vdj_.S23_clones.txt vdj_.S25_clones.txt S23_S25
#java -jar vdjtools-1.1.7/vdjtools-1.1.7.jar  OverlapPair -p vdj_.S22_clones.txt vdj_.S25_clones.txt S22_S25
#java -jar vdjtools-1.1.7/vdjtools-1.1.7.jar  OverlapPair -p vdj_.S23_clones.txt vdj_.S24_clones.txt S23_S24
java -jar vdjtools-1.1.7/vdjtools-1.1.7.jar  OverlapPair -p vdj_.S26_clones.txt vdj_.S27_clones.txt S26_S27
java -jar vdjtools-1.1.7/vdjtools-1.1.7.jar  OverlapPair -p vdj_.S28_clones.txt vdj_.S29_clones.txt S28_S29
java -jar vdjtools-1.1.7/vdjtools-1.1.7.jar  OverlapPair -p vdj_.S22_clones.txt vdj_.S26_clones.txt S22_S26
java -jar vdjtools-1.1.7/vdjtools-1.1.7.jar  OverlapPair -p vdj_.S22_clones.txt vdj_.S27_clones.txt S22_S27
java -jar vdjtools-1.1.7/vdjtools-1.1.7.jar  OverlapPair -p vdj_.S22_clones.txt vdj_.S28_clones.txt S22_S28
java -jar vdjtools-1.1.7/vdjtools-1.1.7.jar  OverlapPair -p vdj_.S22_clones.txt vdj_.S29_clones.txt S22_S29
