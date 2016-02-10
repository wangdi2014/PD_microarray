#Miscellaneous functions used for PD.microarray analysis Nov. 2013

annote.platform <- function(dataMatrix, platform){
  #####Annotate the data with Gene symbol and get rid of data from probes with non-corresponding genes
  if(platform == "HGU133Plus2"){
    # [HG-U133_Plus_2] Affymetrix Human Genome U133 Plus 2.0 Array
    annotation.file <- ("/home/wangdi/PD/Data/platform_annotation/HG-U133_Plus_2.na32.annot/modifiedHG-U133_Plus_2.na32.annot.csv")
  }
  if(platform == "HGU133A" ){
    # [HG-U133A] Affymetrix Human Genome U133A Array
    annotation.file <- ("/home/wangdi/PD/Data/platform_annotation/HG-U133A/modifiedHG-U133A.na33.annot.csv")
    
  }
  if(platform == "HGU133B" ){
    # [HG-U133B] Affymetrix Human Genome U133B Array
    annotation.file <- ("/home/wangdi/PD/Data/platform_annotation/HG-U133B/modifiedHG-U133B.na33.annot.csv")
  }
  if(platform == "U133X3P" ){
      # [U133-X3P] Affymetrix Human Genome X3P Array
      annotation.file <- ("/home/wangdi/PD/Data/platform_annotation/U133_X3P/modifiedU133_X3P.na33.annot.csv")
  }
  if(platform == "HGFocus" ){
    # [HG-Focus] Affymetrix Human HG-Focus Target Array
    annotation.file <- ("/home/wangdi/PD/Data/platform_annotation/HG-Focus/modifiedHG-Focus.na31.annot.csv")
  }
  
  annotations <- read.csv(file=annotation.file,head=TRUE)
  positions <- match(rownames(dataMatrix),annotations$Probe.Set.ID)
  annotations.mat<-as.matrix(annotations)
  
   
  annoteGene <- matrix(data = NA, nrow=nrow(dataMatrix), ncol=1)
  for(i in 1:nrow(dataMatrix)){
    annoteGene[i] <- annotations.mat[positions[i],2]
  }
  rownames(dataMatrix) <- annoteGene
  dataMatrixGene<-dataMatrix[!rownames(dataMatrix) %in% "---",]  #get rid of no gene symbol listed in annotation file
  
  return(dataMatrixGene)
}

### Make a csv file of the stat analysis of the miRNA genes of interest (target)
target.gene.result <- function(result, miRNAlist, studyNumber, dataLocation){
    
  result.target = data.frame()
  for(i in 1:length(miRNAlist$Target)){
    pos <- match(miRNAlist$Target[i],result$X)
    result.target[i,'Target/Gene Sym'] = miRNAlist$Target[i]
    result.target[i,'pvalue'] = result$pvalue[pos]
    result.target[i,'qvalue'] = result$qvalue[pos]
    result.target[i,'Rawfoldchange'] = result$Rawfoldchange[pos]
    result.target[i,'comparison.group.mean'] = result$comparison.group.mean[pos]
    result.target[i,'base.group.mean'] = result$base.group.mean[pos]
  }
  outputfile.target <- paste(studyNumber,"_analysis.TARGET.csv",sep="")
  write.csv(result.target,file=outputfile.target)
  return(result.target)
}

#Plot volcano plot with the specific genes labeled.
myVocalnoPlot.PD <- function(file.name, miRNAlist, dataLocation,studyNumber,studyName,labelFontsize, log2fc.thr = 0.5, pvalue.thr = 0.05, xmin = -2.5, xmax = 2.5, ymin = 0, ymax = 6){
    
    table = as.matrix(read.csv(paste(dataLocation,file.name,sep=""), check.names = TRUE))
    row.names = table[,1]
    table = table[,-1]
    rownames(table) = row.names
    storage.mode(table) = 'double'
    p.value = table[,'pvalue']
    raw.fc = table[,'Rawfoldchange']
    idx = raw.fc < 0
    raw.fc[idx] = -1/raw.fc[idx]
    log2fc = log2(raw.fc)
    ID = row.names
    
    require(ggplot2)
    pv.fc.list = data.frame(p.value, log2fc, ID)
    #  pv.fc.list$threshold = as.factor(abs(pv.fc.list$log2fc) > log2fc.thr & pv.fc.list$p.value < pvalue.thr)
    pv.fc.list.idx = rep(0, length(pv.fc.list$log2fc))
    pv.fc.list.idx[pv.fc.list$log2fc > log2fc.thr & pv.fc.list$p.value < pvalue.thr] = 1
    pv.fc.list.idx[pv.fc.list$log2fc < -log2fc.thr & pv.fc.list$p.value < pvalue.thr] = -1
    
    pv.fc.list$threshold = as.factor(pv.fc.list.idx)
    
    myColors = brewer.pal(3,"Set1")
    cidx = c(3,2,1)
    myColors = myColors[cidx]
    myColors[2]='slategray2'
    myColors[1]='green'
    names(myColors) <- levels(pv.fc.list$threshold)
    colScale <- scale_colour_manual(name = "threshold",values = myColors)
    
    idx = abs(pv.fc.list$log2fc) > log2fc.thr & pv.fc.list$p.value < pvalue.thr
    mlist = table[idx,]
    write.csv(mlist, file = paste(dataLocation,gsub('.csv', '.pdf', file.name),sep=""))
    
    #For labeling the target gene only
    labeled.pv.fc.list = pv.fc.list[pv.fc.list[,3] %in% miRNAlist$Target ,]
    
    
    ##Construct the plot object
    dev = pdf(paste(dataLocation,gsub('.csv', '.pdf', file.name),sep=""))
    
    g<-ggplot(data=pv.fc.list, aes(x=log2fc, y=-log10(p.value), colour=threshold))+
    geom_point()+
    theme(legend.position = "none")+
    xlim(c(xmin, xmax)) + ylim(c(ymin, ymax))+
    xlab("log2 fold change") + ylab("-log10 p-value")
    
    g <- g + colScale
    
    #Following g definition:
    #Add for PD microarray analysis (for Dr.Su) to label the only targeted gene -Nov.2013 mn
    g<- g +
      ggtitle(paste(studyNumber,studyName,sep=" ")) +
      geom_point(data=labeled.pv.fc.list, aes(x=log2fc, y=-log10(p.value)), shape=10, colour='#FFFFFF') +
      geom_text(data=labeled.pv.fc.list, aes(x=log2fc, y=-log10(p.value)),
                label = rownames(labeled.pv.fc.list),
                hjust = 1, vjust=1, colour="black", size=labelFontsize)
    
    
    #  browser()
    #g + geom_text(aes(x=log2fc, y=-log10(p.value), label=ID, size=1.2), colour="black")
    print(g)
    dev.off()
    #  boxplot(main="All Samples' boxplot",AVGSignal.annot.noControl[,clinical.info[,1]], outline = FALSE, las =2, xpd=NA, xaxt='n',ylim=c(2,10) );
    #	text(x=1:ncol(AVGSignal.annot.noControl),y=1.0,labels=colnames(AVGSignal.annot.noControl[,clinical.info[,1]]),srt=90,cex=0.8,xpd=NA);
    #	text(x=1:ncol(AVGSignal.annot.noControl),y=9,labels=labels,srt=90,cex=0.8,xpd=NA);dev.off()
    
}

TargetGenes.Heatmap <- function(target.list, dataMatrix, cases, controls, resultLocation,studyNumber){
  #  browser()
  library(gplots)
  
  Target.dataMatrix = dataMatrix[target.list,c(cases, controls)]
  colors.col = c(rep('red', length(cases)), rep('blue',length(controls)))
  
  library(som)
  names = dimnames(Target.dataMatrix)
  normalize(Target.dataMatrix)
  dimnames(Target.dataMatrix) = names
  dev = pdf(file = paste(resultLocation,studyNumber, '_Clustering.pdf'))
  
  heatmap.2(Target.dataMatrix,ColSideColors=colors.col,col='greenred',
            dendrogram='row', scale='row', trace='none',
            main=paste("Targets of miRNA: ",studyNumber,sep=""), 
            labRow = rownames(Target.dataMatrix),cexRow = 0.9)
            
  
  dev.off()
  
}

