library(data.table)
library(dplyr)

sample_info_file = "/Users/koshka/Documents/fav_new/mixcr_n_test/sample_names.txt"
mixcr_results_dir = "/Users/koshka/Documents/fav_new/mixcr_n_test/mixcr_results"
script_results_dir = "/Users/koshka/Documents/fav_new/mixcr_n_test/script_results"
temp.aln = list.files(mixcr_results_dir,"align.report",full.names=T)
temp.cdr123.ass = list.files(mixcr_results_dir,".all_cdrs.assemble.report",full.names=T)
temp.cdr3.ass = list.files(mixcr_results_dir,".L001.assemble.report",full.names=T)

get_value = function(text,lines){
  line = lines[grep(text,lines)]
  e = regexpr("[0-9]+\\.?[0-9]+",line)
  return(as.numeric(substr(line,e,(e+attr(e,"match.length")-1))))
}
get_value_p = function(text,lines){
  line = lines[grep(text,lines)]
  e = regexpr("[0-9]+\\.?[0-9]+%",line)
  return(as.numeric(substr(line,e,(e+attr(e,"match.length")-2))))
}

aln = mapply(function(x){
  lines = readLines(x)
  total.reads = get_value("Total sequencing reads:",lines)
  aln.reads = get_value("Successfully aligned reads:",lines)
  aln.reads.p = get_value_p("Successfully aligned reads:",lines)
  aln.pair.conflicts.p = get_value_p("Paired-end alignment conflicts eliminated:",lines)
  aln.igh.chains.p = get_value_p("IGH chains:",lines)
  aln.igk.chains.p = get_value_p("IGK chains:",lines)
  aln.igl.chains.p = get_value_p("IGL chains:",lines)
  sample = sub("/.*/","",x)
  sample = gsub("_.*","",sample)
  return(data.frame(sample=sample,total.reads=total.reads,aln.reads=aln.reads,
                    aln.reads.p=aln.reads.p,aln.pair.conflicts.p=aln.pair.conflicts.p,
                    aln.igh.chains.p=aln.igh.chains.p,aln.igk.chains.p=aln.igk.chains.p,
                    aln.igl.chains.p=aln.igl.chains.p))
},temp.aln,SIMPLIFY = F) %>% rbindlist

read_assemble = function(x){
  lines = readLines(x)
  ass.reads.used = get_value("Reads used in clonotypes, percent of total:",lines)
  ass.reads.used.p = get_value_p("Reads used in clonotypes, percent of total:",lines)
  ass.pcr.correction.p = get_value_p("Reads clustered in PCR error correction, percent of used:",lines)
  ass.no.target.seq.p = get_value_p("Reads dropped due to the lack of a clone sequence:",lines)
  ass.failed.mapping.p = get_value_p("Reads dropped due to failed mapping:",lines)
  ass.cln.count = get_value("Final clonotype count:",lines)
  ass.reads.per.cln = get_value("Average number of reads per clonotype:",lines)
  sample = gsub("/.*/","",x)
  sample = gsub("_.*","",sample)
  return(data.frame(sample=sample,ass.reads.used=ass.reads.used,
                    ass.reads.used.p=ass.reads.used.p,
                    ass.pcr.correction.p=ass.pcr.correction.p,
                    ass.no.target.seq.p=ass.no.target.seq.p,
                    ass.failed.mapping.p=ass.failed.mapping.p,
                    ass.cln.count=ass.cln.count,ass.reads.per.cln=ass.reads.per.cln))
}
ass.cdr123 = mapply(read_assemble,temp.cdr123.ass,SIMPLIFY = F) %>% rbindlist
ass.cdr3 = mapply(read_assemble,temp.cdr3.ass,SIMPLIFY = F) %>% rbindlist

ass = merge(ass.cdr3,ass.cdr123,by = "sample",suffixes = c(".cdr3",".cdr123"))

report.table = merge(aln,ass,by = "sample")
############# write ################
fwrite(report.table,paste(script_results_dir,"stats.mixcr.txt",sep="/"),sep="\t")

############## pics ################
sample_info = fread(sample_info_file) %>% mutate(sample = sub("_.*","",sample)) %>%
  select(sample,run_id,subject_id,cell_type)
report.table = report.table %>% merge(sample_info,by = "sample")

# check reads quantity
df1 = report.table %>% 
  select(sample,run_id,subject_id,cell_type,total.reads,aln.reads,ass.reads.used.cdr3,ass.reads.used.cdr123)
df1.p = df1 %>% mutate(aln.reads = aln.reads/total.reads,
                       ass.reads.used.cdr3 = ass.reads.used.cdr3/total.reads,
                       ass.reads.used.cdr123 = ass.reads.used.cdr123/total.reads,
                       total.reads = 1)
df1 = df1 %>% melt(id = c("sample","run_id","subject_id","cell_type"))
df1.p = df1.p %>% melt(id = c("sample","run_id","subject_id","cell_type"))

r1.1 = ggplot(df1,aes(variable,value,group=sample,color=subject_id,linetype=cell_type))+
  geom_line()+theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab("")
r1.2 = ggplot(df1.p,aes(variable,value,group=sample,color=factor(run_id)))+
  geom_line()+theme_minimal()+scale_y_continuous(limits = c(0,1))+
  scale_color_discrete(name="runID")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab("")

# IHG/IGK/IGL
df2 = report.table %>% 
  select(sample,subject_id,cell_type,aln.igh.chains.p,
         aln.igk.chains.p,aln.igl.chains.p) %>%
  mutate(sample.id = as.integer(sub("Lom","",sample)))
sample_levels = (df2 %>% arrange(sample.id))$sample
df2 = df2 %>% melt(id = c("sample","sample.id","subject_id","cell_type"))

r2 = ggplot(df2,aes(factor(sample,levels=sample_levels),value,fill=variable))+
  geom_col()+theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab("sample")+ylab("%")

# problems %
df3 = report.table %>% 
  select(sample,run_id,subject_id,cell_type,aln.pair.conflicts.p,
         ass.pcr.correction.p.cdr3,ass.pcr.correction.p.cdr123,
         ass.no.target.seq.p.cdr3,ass.no.target.seq.p.cdr123,
         ass.failed.mapping.p.cdr3,ass.failed.mapping.p.cdr123)
df3 = df3 %>% melt(id = c("sample","run_id","subject_id","cell_type"))

r3 = ggplot(df3,aes(variable,value,color=subject_id,shape=as.factor(run_id)))+
  geom_jitter(height=0)+theme_minimal()+scale_shape_manual(values = c(1,19),name="runID")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab("problems")

## final clonotypes
df4.1 = report.table %>% 
  select(sample,run_id,subject_id,cell_type,
         ass.cln.count.cdr3,ass.cln.count.cdr123)
df4.1 = df4.1 %>% melt(id = c("sample","run_id","subject_id","cell_type"))

r4.1 = ggplot(df4.1,aes(variable,value,color=cell_type,shape=as.factor(run_id)))+
  geom_jitter(height=0,width=0.2)+theme_minimal()+scale_shape_manual(values = c(1,19),name="runID")+
   theme(axis.text.x = element_text(angle = 90, hjust = 1))

r4.2 = ggplot(report.table,aes(total.reads,ass.cln.count.cdr3,color=cell_type,shape=as.factor(run_id)))+
  geom_point()+theme_minimal()+scale_shape_manual(values = c(1,19),name="runID")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

library(gridExtra)
r1 = grid.arrange(r1.1,r1.2,ncol=2)
r4 = grid.arrange(r4.1,r4.2,ncol=2)
ggsave(paste(script_results_dir,"report.r1.reads.pdf",sep="/"),r1,height=6,width=12)
ggsave(paste(script_results_dir,"report.r2.chains.pdf",sep="/"),r2,height=6,width=8)
ggsave(paste(script_results_dir,"report.r3.problems.pdf",sep="/"),r3,height=6,width=8)
ggsave(paste(script_results_dir,"report.r4.other.pdf",sep="/"),r4,height=6,width=12)

