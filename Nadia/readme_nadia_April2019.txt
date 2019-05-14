############### I. using Mixcr to get BCR repertoires ##############################

makeclones2.sh
(details inside the script, but you need input /fasta dir, and you will get /mixcr_results dir)
(launch line: qsub -l hostname=hpc03 makeclones2.sh)

note that we use following options in MiXCR assemble:
-OaddReadsCountOnClustering=true (add the number of seqs with corrected PCR errors to parent target sequence count)
-OseparateByV=true (in CDR3 assemble, we separate V genes, after we will define best V gene ourselves)
-OassemblingFeatures="[CDR1+FR2+CDR2+FR3+CDR3]" (we use it to get something very close to full transcript of the clone)

we use -aaMutation option in exportClones to record mutations from germline 

as a result we get: 
	alignment and assemble reports,
	CDR3 clones
	CDR1-CDR3 clones, i.e full transcripts

we will analyze this in downstream analysis

## ALL scripts are in the /scripts dir and analyze data from /mixcr_results dir or /script_results, they write data and pics to script_results ######

############## II. analyze report script ###############

report_table.R

!!! specify yourself sample_info_file,mixcr_results_dir, script_results_dir FULL pathways at the beginning of the script !!!

it will generate stats.mixcr.txt with data collected from mixcr reports in out_dir
also 4 report pics report*.pdf

############## III. filter script ######################

filters_stats.R

!!! specify yourself sample_info_file,mixcr_results_dir,script_results_dir FULL pathways at the beginning of the script !!!

it will generate filter pic and 2 contamination pics, stat.filter.txt and 4 types of filtered data for downstream analysis:
-all.filtered.txt
-all.filtered.ST.txt
-all.filtered.clean.txt
-all.filtered.clean.ST.txt

!!! NOTE that all further scripts will read this one of this files as an input instead of mixcr_results files!!!
############## IV. distance to germline ################
############## V. clonality   ##########################
############## VI. shared clones #######################
############## VII. trees of shared clones #############
