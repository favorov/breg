#! /bin/bash
t_coffee_dir=../t_coffee
for fasta in *fasta
do
	${t_coffee_dir}/t_coffee ${fasta} -output=fasta_aln,html
done
