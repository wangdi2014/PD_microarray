# -----Parkinson's Disease microassay analysis-----
# This is the main function of PD.analysis that runs all functions.
#
# Methods: rma normalization on raw data and unpaired t-test
#
# required libraries: biocLite - affy, (affyPLM), qvalue, qtools, preprocessCore 
#                     RColorBrewer
# Created: Nov-05-13, Modified: Nov-11-13

#####INPUTS######
#####import miRNAlist (gene list) from the researcher
miRNAlist <- read.csv(file="/home/wangdi/PD/042914_ParkinsonsProject_for_Validation/Validation_Test_Data/Data/miRNAofInterest.csv",head=TRUE)

#path <- "/home/wangdi/PD/042914_ParkinsonsProject_for_Validation/Validation_Test_Data/"
#miRNAlistLoc <- paste(path,"Data/miRNAofInterest.csv",sep="")
#miRNAlist <- read.csv(file=miRNAlistLoc,head=TRUE)

##### GSE20292: Zhang study ##### 
#
# the GEO accession number of the study and its platform ("HGU133Plus2" or "HG133A")
studyNumber <- "GSE20292" 
studyName <-"Zhang"   
platform <- "HGU133A"

# control and case group sample names:
controls <- c("GSM508708","GSM508717","GSM508720","GSM508721","GSM508722","GSM508723","GSM508724","GSM508725","GSM508726",
              "GSM508729","GSM508730","GSM508733","GSM508734","GSM508735","GSM521253","GSM606624","GSM606625","GSM606626")

cases <- c( "GSM508710","GSM508711","GSM508712","GSM508713","GSM508714","GSM508715","GSM508716","GSM508718","GSM508728",
            "GSM508731","GSM508732")


#####--------------Code starts here-------------------#####

### Set directory to import functions
#codesLoc <- paste(path,"Codes/",sep="")
#setwd(codesLoc)
setwd("/home/wangdi/PD/042914_ParkinsonsProject_for_Validation/Validation_Test_Data/Codes/")
source("PD.gene.exp.process.R")
source("pvfoldchange.R")
source("PD.miscellaneous.R")

### Process raw data using RMA method to get log 2 base values of the expression levels
# raw data is saved as .Rda
PD.gene.exp.process(studyNumber)
dataLocation <- paste("/home/wangdi/PD/042914_ParkinsonsProject_for_Validation/Validation_Test_Data/Data/",studyNumber,"_RAW/",sep="")
setwd(dataLocation)
file=paste(studyNumber,".Rda",sep="")
load(file)  #load dataMatrix

### Annotate the data with correspoding Gene symbol 
# (conversion of probe labels to gene symbols)
dataMatrixGene <- annote.platform(dataMatrix, platform)

### Statistical analysis 
# Qunatile normalizations on dataset and t-test are conducted
# Produce and save the result as a csv file
outputfile <- paste(studyNumber,"_analysis.csv",sep="")
pvfoldchange(dataMatrixGene,cases,controls,quantiles.normalize = TRUE,gene.level = TRUE, paired.ttest = FALSE, outputfile)

### Make a csv file of the stat. analysis results only for the target gene (from miRNAlist)
result.all <- read.csv(file=paste(dataLocation,outputfile,sep=""),head=TRUE)
result.target <- target.gene.result(result.all, miRNAlist, studyNumber, dataLocation)


### Make Volcanot plot of all genes with the labels of the target genes
library(RColorBrewer)
myVocalnoPlot.PD(outputfile, miRNAlist, dataLocation,studyNumber,studyName, labelFontsize=3)
