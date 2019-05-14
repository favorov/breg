library(tibble)
library(dplyr)
library(data.table)
library(ggplot2)

### functions ###

write_with_CF_recount = function(df,file.name){
  df.sums = df %>% group_by(sample) %>% summarise(total = sum(cloneCount))
  df = df %>% merge(df.sums,by="sample") %>% mutate(cloneFraction = cloneCount/total)
  df = df %>% select(sample,subject_id,cell_type,cloneCount,cloneFraction,
                     targetSequences,bestVGene,aaSeqCDR3,
                     nMutations.total,aaMutations.total,mut.n,mut.aa.n)
  fwrite(df,file.name,sep="\t")
}
#################

sample_info_file = "/Users/koshka/Documents/fav_new/mixcr_n_test/sample_names.txt"
mixcr_results_dir = "/Users/koshka/Documents/fav_new/mixcr_n_test/mixcr_results"
script_results_dir = "/Users/koshka/Documents/fav_new/mixcr_n_test/script_results"
temp = list.files(mixcr_results_dir,"nMut.IGH.txt",full.names=T)

df = mapply(function(x){
  return(rownames_to_column(fread(x) %>%
                              mutate(sample = gsub("/.*/","",x) %>% gsub(pattern = "_.*",replacement = ""),
                                     nMutations.total = paste(nMutationsCDR1,nMutationsFR2,nMutationsCDR2,nMutationsFR3,sep=","),
                                     aaMutations.total = paste(aaMutationsCDR1,aaMutationsFR2,aaMutationsCDR2,aaMutationsFR3,sep=","),
                                     mut.n = sapply(nMutations.total,function(x){
                                       grex = gregexpr("[A-Z0-9]+",x)[[1]]
                                       return(ifelse(grex[1]==-1,0,length(grex)))
                                     }),
                                     mut.aa.n = sapply(aaMutations.total,function(x){
                                       grex = gregexpr("[A-Z0-9]+",x)[[1]]
                                       return(ifelse(grex[1]==-1,0,length(grex)))
                                     }))))
},temp,SIMPLIFY = F) %>% rbindlist

sample_info = fread(sample_info_file) %>% mutate(sample = sub("_.*","",sample)) %>%
  select(sample,subject_id,cell_type)
df = df %>% merge(sample_info,by = "sample")

reads1 = sum(df$cloneCount)
targets1 = df %>% nrow
clns1  = length(unique(df$aaSeqCDR3))

f1 = df %>% group_by(sample) %>% summarize(reads = sum(cloneCount),targets=n(),
                                           clones = length(unique(aaSeqCDR3)))

############## productive filter ##############
df = df %>% filter(!grepl("\\*",aaSeqCDR3))
df = df %>% filter(!grepl("_",aaSeqCDR3))
df = df %>% filter(aaSeqCDR3!="")
targets2 = df %>% nrow
reads2 = sum(df$cloneCount)
clns2  = length(unique(df$aaSeqCDR3))
print(paste("Reads after productive filter:",reads2/reads1))
print(paste("Transcripts after productive filter:",targets2/targets1))
print(paste("Clones after productive filter:",clns2/clns1))

f2= df %>% group_by(sample) %>% summarize(reads = sum(cloneCount),targets=n(),
                                           clones = length(unique(aaSeqCDR3)))


############## best V gene (anti-chimeric) filter ##############
df.vgene = df %>% group_by(sample,aaSeqCDR3) %>% arrange(desc(cloneCount)) %>% summarise(bestVGene = first(bestVGene))
df = df %>% merge(df.vgene,by = c("sample","aaSeqCDR3","bestVGene"))
targets3 = df %>% nrow
reads3 = sum(df$cloneCount)
clns3  = length(unique(df$aaSeqCDR3))
print(paste("Reads after best V gene (anti-chimeric) filter:",reads3/reads2))
print(paste("Transcripts after best V gene (anti-chimeric) filter:",targets3/targets2))
print(paste("Clones after best V gene (anti-chimeric) filter:",clns3/clns2))

f3= df %>% group_by(sample) %>% summarize(reads = sum(cloneCount),targets=n(),
                                          clones = length(unique(aaSeqCDR3)))

############## try singletones (for clone!) #########
df.no.singles = df %>% group_by(sample,aaSeqCDR3) %>% filter(sum(cloneCount)!=1) %>% ungroup
print(paste("Reads after best V gene (anti-chimeric) filter:",
            sum(df.no.singles$cloneCount)/sum(df$cloneCount)))
print(paste("Targets after best V gene (anti-chimeric) filter:",nrow(df.no.singles)/nrow(df)))
print(paste("Clones after best V gene (anti-chimeric) filter:",
            length(unique(df.no.singles$aaSeqCDR3))/length(unique(df$aaSeqCDR3))))

f3.s= df.no.singles %>% group_by(sample) %>% summarize(reads = sum(cloneCount),targets=n(),
                                          clones = length(unique(aaSeqCDR3)))

############ write filtered data ########################################
write_with_CF_recount(df,paste(script_results_dir,"all.filtered.txt",sep="/"))
write_with_CF_recount(df.no.singles,paste(script_results_dir,"all.filtered.ST.txt",sep="/"))

############ first conamination filter ####################
df.shared = df %>% group_by(aaSeqCDR3) %>% summarise(subjects = paste(unique(subject_id),collapse=","),
                                                     subjects.n = length(unique(subject_id)))
bad.clones.df = df.shared %>% filter(subjects.n>1) %>% arrange(desc(subjects.n))
bad.clones = unique(bad.clones.df$aaSeqCDR3)
print(paste("Number of bad clones (shared between subjects):",length(bad.clones)))

#ranks
ranks = df %>% group_by(sample) %>% arrange(desc(cloneCount)) %>% 
  mutate(rank = c(1:length(sample)))
ranks %>% filter(aaSeqCDR3 %in% bad.clones) %>% select(sample,aaSeqCDR3,rank,cloneCount,cloneFraction) %>%
  arrange(sample,desc(cloneCount)) %>% filter(rank<=20) %>% as.data.frame()

#pics
df.sh = bad.clones.df %>% select(aaSeqCDR3) %>% merge(df %>% group_by(sample,aaSeqCDR3) %>% 
                                                        summarise(cloneFraction = sum(cloneFraction)),by="aaSeqCDR3") %>%
  arrange(aaSeqCDR3,desc(cloneFraction))
df.sh.cdr = df.sh %>% group_by(aaSeqCDR3) %>% 
  summarize(cloneFraction =  sum(cloneFraction)) %>%
  arrange(desc(cloneFraction))
top_bad = df.sh.cdr[1:8,]$aaSeqCDR3
cont_cdr3_levels = c(top_bad,"Other")
df.sh = df.sh %>% mutate(aaSeqCDR3.2 = ifelse(aaSeqCDR3 %in% top_bad,aaSeqCDR3,"Other"),
                         cont.class = ifelse(aaSeqCDR3 %in% top_bad,1,2),
                         sample.id = as.integer(sub("Lom","",sample)))
sample_levels = (df.sh %>% arrange(sample.id))$sample %>% unique
cont.c1 = ggplot(df.sh,aes(x=factor(sample,levels = sample_levels),y=cloneFraction,
                             color = factor(aaSeqCDR3.2,levels = cont_cdr3_levels),
                                            shape=as.factor(cont.class)))+geom_point()+
  scale_color_discrete(name = "top 8 cont CDR3")+scale_y_log10()+
  scale_shape_manual(values = c(19,1),labels = c("top 8","Others"),
                     name="cont.class")+theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab("sample")+
  facet_wrap(~cont.class)

df.sh.sum = df.sh %>% group_by(sample) %>% summarise(cloneFraction.sum = sum(cloneFraction)) %>%
  mutate(sample.id = as.integer(sub("Lom","",sample)))
df.sh.exception = df.sh %>% filter(!(aaSeqCDR3=="CARGTYYDFWSGYLIHDYYYGMDVW" & sample %in% c("Lom5","Lom6","Lom11","Lom12")))
df.sh.exception.sum = df.sh.exception %>% group_by(sample) %>% summarise(cloneFraction.sum = sum(cloneFraction)) %>%
  mutate(sample.id = as.integer(sub("Lom","",sample)))
df.sh.sum$type="normal"
df.sh.exception.sum$type="withSidorovException\n(we know parent clone\nfor contamination)"
df.sh.sum = rbind(df.sh.sum,df.sh.exception.sum)
cont.c2 = ggplot(df.sh.sum)+geom_point(aes(x=factor(sample,levels=sample_levels),y=cloneFraction.sum,shape=type))+
  theme_minimal()+scale_shape_manual(values = c(1,19))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +xlab("sample")+ylab("Fraction of contaminating clones")

ggsave(paste(script_results_dir,"contamination.c1.pdf",sep="/"),cont.c1,height=6,width=12)
ggsave(paste(script_results_dir,"contamination.c2.pdf",sep="/"),cont.c2,height=6,width=8)

############ conamination filter from old paper ####################
e = c("b7","IGHV1-2*04","VRGSTYSPSGYFEY",
   "g10","IGHV1-3*01","ARIFEGLSGIAAPFDY",
   "h1","IGHV1-8*01","AREVSDYSDYGDVYYMDV",
     "e11","IGHV1-8*01","AREVSDYSDYGDVYYMDV",
   "h11","IGHV1-18*04","ARREEGLYTTSPGYFGV",
   "e12","IGHV1-46*01","ARRGFDY",
   "c11","IGHV1-46*03","AKDLRPRDIGDMDV",
   "c12","IGHV1-69*06","ARCGILRSHYFYGMDV",
   "a6","IGHV3-7*01","VRGGLGAGADY",
   "c3","IGHV4-b*01","AGLTQSSHNDAN")
names = e[seq(1, length(e), 3)]
vgene = e[seq(2, length(e), 3)]
seq =  e[seq(3, length(e), 3)]
exp = data.frame(names=names,vgene = vgene,seq = seq)
exp = exp %>% mutate(seq2 = paste("C",seq,"W",sep=""))

df.old = df %>% filter(aaSeqCDR3 %in% exp$seq2)
bad.old.clones = unique(df.old$aaSeqCDR3)
print(paste("Number of bad clones from old paper:",length(bad.old.clones)))

df = df %>% filter(!(aaSeqCDR3 %in% bad.clones),!(aaSeqCDR3 %in% bad.old.clones))
df.no.singles = df.no.singles %>% filter(!(aaSeqCDR3 %in% bad.clones),!(aaSeqCDR3 %in% bad.old.clones))

#pics
f4.s= df.no.singles %>% group_by(sample) %>% summarize(reads = sum(cloneCount),targets=n(),
                                                       clones = length(unique(aaSeqCDR3)))

f = merge(f1,f2,by="sample",suffixes = c("", ".productive")) %>% 
  merge(f3,by="sample",suffixes = c("", ".v.chimeras")) %>%
  merge(f3.s,by="sample",suffixes = c("", ".singltone.cdr3")) %>%
  merge(f4.s,by="sample",suffixes = c("", ".contaminations"))
fm = melt(f,id="sample") %>% mutate(type= sapply(as.character(variable),function(x){
  strsplit(x,"\\.")[[1]][1]}))
fm = fm %>% mutate(filter=mapply(function(x,y){
  res=gsub(x,"",y)
  ifelse(res=="","orig",res)
  },type,variable)) %>% mutate(filter = sub("\\.","",filter))

f_levels = c("orig","productive","v.chimeras","singltone.cdr3","contaminations")
fm = fm %>% merge(sample_info,by="sample")
f.orig= melt(f1,id="sample")
colnames(f.orig) = c("sample","type","total")
fm.p = fm %>%select(-variable) %>% merge(f.orig,by=c("sample","type")) %>% mutate(value.p=value/total)
fm.p.m = melt(fm.p,measure.vars=c("value","value.p"))
filt.f1a = ggplot(fm.p)+geom_line(aes(factor(filter,levels=f_levels),value,
                         group=sample,color=subject_id,linetype=cell_type))+
  facet_wrap(~factor(type,levels=c("reads","targets","clones")),scales = "free")+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab("filters")+
  theme(legend.position="none")
filt.f1b = ggplot(fm.p)+geom_line(aes(factor(filter,levels=f_levels),value.p,
                           group=sample,color=subject_id,linetype=cell_type))+
  facet_wrap(~factor(type,levels=c("reads","targets","clones")))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+xlab("filters")
library(gridExtra)
filt.f1 = grid.arrange(filt.f1a,filt.f1b,nrow=2)
ggsave(paste(script_results_dir,"filter.f1.pdf",sep="/"),filt.f1,height=10,width=12)


############ write filtered data without contamination ###########################
write_with_CF_recount(df,paste(script_results_dir,"all.filtered.clean.txt",sep="/"))
write_with_CF_recount(df.no.singles,paste(script_results_dir,"all.filtered.ST.clean.txt",sep="/"))

############ write stats #########################################################
fwrite(f,paste(script_results_dir,"stats.filters.txt",sep="/"),sep="\t")

####### END ########
####################
## p.s ############# 
# library(stringdist)
# all = unique(df$aaSeqCDR3)
# d = stringdistmatrix(exp$seq2,all,method="lv")
# d.ind = t(which(d<=2,arr.ind=T))
# exp$seq2[d.ind[1,]]
# zasor = all[d.ind[2,]]
# sum((df %>% filter(aaSeqCDR3 %in% zasor))$cloneFraction)

############ top 20 - all nMut big? ####################

############ shared clone Bcell/Breg// Sirokin otdelno/ Sidirov with time - dynamic population ####################
############ control!

