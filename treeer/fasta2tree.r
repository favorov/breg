library(ape)
library(phangorn)

plottree<-function(fastaname,pdfname){
clone<-read.FASTA(fastaname)
cloneDat<-phyDat(clone)
pdf(pdfname)
plot(nj(dist.ml(cloneDat)),main=fastaname,direction = "rightwards")
#edgelabels(t$edge.length, bg="black", col="white", font=2)
dev.off()
}
