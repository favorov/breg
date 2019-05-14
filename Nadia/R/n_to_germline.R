library(data.table)
library(ggplot2)
library(dplyr)

in_file = "/Users/koshka/Documents/fav_new/mixcr_n_test/script_results/all.filtered.txt"
in_file = "/Users/koshka/Documents/fav_new/mixcr_n_test/script_results/all.filtered.ST.txt"
in_file = "/Users/koshka/Documents/fav_new/mixcr_n_test/script_results/all.filtered.clean.txt"
in_file = "/Users/koshka/Documents/fav_new/mixcr_n_test/script_results/all.filtered.ST.clean.txt"
script_results_dir = "/Users/koshka/Documents/fav_new/mixcr_n_test/script_results"

df = fread(in_file)

df.prop.n = df %>% group_by(sample) %>% mutate(sample.sum = n()) %>% ungroup %>%
  group_by(sample,subject_id,cell_type,mut.n) %>% summarise(mut.n.p = n()/first(sample.sum))

df.prop.aa = df %>% group_by(sample) %>% mutate(sample.sum = n()) %>% ungroup %>%
  group_by(sample,subject_id,cell_type,mut.aa.n) %>% summarise(mut.aa.p = n()/first(sample.sum))

#n1.1 = ggplot(df,aes(mut.aa.n,fill=cell_type))+geom_histogram(position = "dodge")+
# facet_wrap(~subject_id,ncol=3)+scale_y_log10()+theme_minimal()

#n1.2 = ggplot(df.prop.aa,aes(x = mut.aa.n,y = mut.aa.p,group=cell_type,color = cell_type))+geom_line()+
#facet_wrap(~subject_id,ncol=3)+scale_y_log10()+theme_minimal()

n1 = ggplot(df.prop.aa,aes(x = mut.aa.n,y = mut.aa.p,group=cell_type,color = cell_type))+geom_point()+
  facet_wrap(~subject_id,nrow = 1)+theme_minimal()+theme(legend.position = "bottom")+
  xlab("# of AA substitutions to gremline (in V transcript)")+
  ylab("")

# n2 = ggplot(df.prop.aa,aes(mut.aa.n,mut.aa.p,group=cell_type,color=cell_type))+
#   stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),geom="pointrange")+
#   theme_minimal()+theme(legend.position = "bottom")+
#   xlab("# of AA substitutions to gremline (in V transcript)")+
#   ylab("")

df.prop.aa.minus.Solovjev = df.prop.aa %>% filter(subject_id!="Solovjev")

# n2.minus.Solovjev = ggplot(df.prop.aa.minus.Solovjev,aes(mut.aa.n,mut.aa.p,group=cell_type,color=cell_type))+
#   stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),geom="pointrange")+
#   theme_minimal()+theme(legend.position = "bottom")+
#   xlab("# of AA substitutions to gremline (in V transcript)")+
#   ylab("")

df.prop.aa.compare = rbind(df.prop.aa %>% mutate(type="all"),
                           df.prop.aa.minus.Solovjev %>% mutate(type="without.Solovjev"))

n2 = ggplot(df.prop.aa.compare,aes(mut.aa.n,mut.aa.p,group=cell_type,color=cell_type))+
  stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),geom="pointrange")+
  facet_wrap(~type)+
  theme_minimal()+theme(legend.position = "bottom")+
  xlab("# of AA substitutions to gremline (in V transcript)")+
  ylab("fraction of all unique transcripts")

naive.table.0 = df.prop.aa.compare %>% filter(mut.aa.n==0) %>% group_by(sample,cell_type,type,subject_id) %>%
  summarise(naive.p = sum(mut.aa.p))
naive.table.0$naive.threshold="SHM==0"
naive.table.1 = df.prop.aa.compare %>% filter(mut.aa.n<=1) %>% group_by(sample,cell_type,type,subject_id) %>%
  summarise(naive.p = sum(mut.aa.p))
naive.table.1$naive.threshold="SHM<=1"
naive.table=rbind(naive.table.0,naive.table.1)

n3 = ggplot(naive.table,aes(cell_type,naive.p,color=cell_type))+geom_boxplot()+geom_point()+
  facet_grid(type~naive.threshold)

naive.tbl = naive.table.1 %>% filter(type=="all",naive.threshold=="SHM<=1") %>% 
  dcast(type+subject_id+naive.threshold~cell_type,value.var="naive.p",
        fun.aggregate = mean)
fwrite(naive.tbl,paste(script_results_dir,"distGermline.naive.p.txt",sep="/"),sep="\t")

ggsave(paste(script_results_dir,"distGermline.n1.pdf",sep="/"),n1,height=6,width=14)
ggsave(paste(script_results_dir,"distGermline.n2.pdf",sep="/"),n2,height=6,width=12)
ggsave(paste(script_results_dir,"distGermline.n3.pdf",sep="/"),n3,height=6,width=6)


#n4.1 = ggplot(df,aes(mut.aa.n,cloneFraction,group=cell_type,color=cell_type))+
#  geom_point(alpha=0.5)+
#  facet_wrap(~subject_id,nrow = 2)
n4.2 = ggplot(df,aes(mut.aa.n,cloneFraction,group=cell_type,color=cell_type))+
  geom_point(alpha=0.5)+
  facet_wrap(~subject_id,nrow = 2)+scale_y_log10()+
  xlab("# of AA substitutions to gremline (in V transcript)")
n4.2.2 = ggplot(df,aes(mut.aa.n,cloneFraction,group=cell_type,color=cell_type))+
  geom_point(alpha=0.5)+
  facet_grid(subject_id~cell_type)+scale_y_log10()+
  xlab("# of AA substitutions to gremline (in V transcript)")


# df.s = df %>% filter(subject_id!="Solovjev")
# df.compare = rbind(df %>% mutate(type="all"),df.s %>% mutate(type="without.Solovjev"))
# n4.3 = ggplot(df.compare,aes(mut.aa.n,cloneFraction,group=cell_type,color=cell_type))+
#   stat_summary(fun.data=mean_sdl, fun.args = list(mult=1),geom="pointrange")+scale_y_log10()+
#   facet_wrap(~type)


#ggsave(paste(script_results_dir,"distGermline.n4.1.pdf",sep="/"),n4.1,height=8,width=10)
ggsave(paste(script_results_dir,"distGermline.n4.2.pdf",sep="/"),n4.2,height=8,width=10)
ggsave(paste(script_results_dir,"distGermline.n4.2.2.pdf",sep="/"),n4.2.2,height=10,width=8)

############################################################################################
###################### new contamination investigation #####################################
######## we leave that for later! ########### doesnt really change our conclusions##########
ggplot(df %>% filter(subject_id=="Sidorov"),aes(mut.aa.n,cloneFraction,group=cell_type,color=cell_type))+
  geom_point(alpha=0.5)+
  facet_wrap(~sample,nrow = 2)+scale_y_log10()

ggplot(df %>% filter(subject_id=="Bumba"),aes(mut.aa.n,cloneFraction,group=cell_type,color=cell_type))+
  geom_point(alpha=0.5)+
  facet_wrap(~sample,nrow = 2)+scale_y_log10()

df %>% filter(subject_id=="Sidorov",mut.aa.n==11) %>% select(aaSeqCDR3,aaMutations.total)
df %>% filter(subject_id=="Bumba",mut.aa.n==11) %>% select(aaSeqCDR3,aaMutations.total)
df %>% filter(subject_id=="Engalitseva",mut.aa.n==11) %>% select(aaSeqCDR3,aaMutations.total)
 
df %>% filter(aaMutations.total == 'SD1G,SG6S,S*0W,SI1S,SH1Y,SY0N,SI9V,SM11I,ST15K,SY21S') %>%
  select(subject_id,aaSeqCDR3,aaMutations.total)

#ggplot(df,aes(mut.aa.n,cloneFraction,group=cell_type,color=cell_type))+geom_jitter(alpha=0.5)+
  #facet_grid(subject_id~cell_type)+scale_y_log10()

