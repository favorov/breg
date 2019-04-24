library(ape)
library(phangorn)

plottree<-function(fastaname,pdfname=NA){
  if (is.na(pdfname))
    pdfname<-paste0(gsub("\\..*$", "",fastaname),".pdf")
  clone<-read.FASTA(fastaname)
  cloneDat<-phyDat(clone)
  tree<-nj(dist.hamming(cloneDat))
  pdf(pdfname)
  plot(tree,main=fastaname,direction = "rightwards")
  edgelabels(tree$edge.length, bg="black", col="white", font=.5)
  dev.off()
}

sapply(list.files(pattern="*.fasta_aln"),plottree)