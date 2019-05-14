library(ape)
library(phangorn)

plottree<-function(fastaname,pdfname=NA){
  if (is.na(pdfname))
    pdfname<-paste0(gsub("\\..*$", "",fastaname),".pdf")
  clone<-read.FASTA(fastaname)
  #cloneDat<-phyDat(clone)
  #tree<-nj(dist.hamming(cloneDat))
  tree<-nj(adist(as.character(clone)))
  pdf(pdfname)
  plot(tree,main=fastaname,direction = "rightwards",cex=.5)
  edgelabels(round(tree$edge.length), bg="blue", col="white",font=2,frame="circle",cex=.5)
  dev.off()
}

sapply(list.files(pattern="*.fasta_aln"),plottree)