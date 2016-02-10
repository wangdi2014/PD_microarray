# THIS CODE extracts the expresssion level (in log 2 base) of each probe using RMA method
# of given .CLE files

PD.gene.exp.process <- function(studyNumber){

setwd(paste("/home/wangdi/PD/042914_ParkinsonsProject_for_Validation/Validation_Test_Data/Data/",studyNumber,"_RAW/",sep=""))
  
source("http://bioconductor.org/biocLite.R")
#biocLite('hgu133a.db')
#biocLite('hgu133plus2.db')    #annotation database for Affymetrix Human Genome U133 Plus 2.0 Array
#browser()

###### install affy from bioconductor
library (affy)
######
Data <- ReadAffy(verbose = TRUE) 
#A##### will download cdf files automatically
#Data
#######
#Normalization default rma
eset <- rma(Data)

#remove '.CEL' extension from row name
dataMatrix=exprs(eset)
x=dimnames(dataMatrix)[[2]]
y=strsplit(x,split=".CEL|.cel")
COLCHG=unlist(y)
dimnames(dataMatrix)[[2]]= COLCHG


#remove Affymetrix's probe controls from dataset 
y=grep("^AFFX",dimnames(dataMatrix)[[1]])
dataMatrix2=dataMatrix[-c(y),]
dataMatrix=dataMatrix2

#save dataMatrix and image/project file in .Rna
save(dataMatrix,file=paste(studyNumber,".Rda",sep=""))



library(affyPLM)
dataPLM=fitPLM(Data)
dev = pdf(paste(studyNumber,'_RLE.pdf',sep=""))
Mbox(dataPLM, main = "RLE", ylim = c(-0.4,0.4),outline = FALSE, col = "lightblue",las = 3, whisklty = 0,cex.axis=0.6,col.axis="azure4", staplelty = 0,pars = list(boxwex = 0.6, staplewex = 0.5, outwex = 0.5))
abline(h=0.05,col="red")
abline(h=-0.05,col="red")
dev.off()

dev = pdf(paste(studyNumber,'_NUSE.pdf',sep=""))
boxplot(dataPLM, main = "NUSE", ylim = c(0.95,1.22), outline = FALSE, col = "lightblue",las = 3,cex.axis=.6,col.axis="azure4", whisklty = 0, staplelty = 0,pars = list(boxwex = 0.6, staplewex = 0.5, outwex = 0.5))
dev.off()
save.image(file=paste(studyNumber,".Rdata",sep=""))


}

