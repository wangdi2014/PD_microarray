 #pvfoldchange(dataMatrix,comp.group.names=recalcitrant,base.group.names=healing,FALSE,FALSE,FALSE,file.name='Recalcitrant_VS._Healing_Visit_First.csv')

'pvfoldchange' <- function(dataMatrix, comp.group.names, base.group.names, 
													quantiles.normalize = TRUE,gene.level = TRUE, paired.ttest = FALSE, file.name){
#	browser()
	library(preprocessCore)
	if(quantiles.normalize){
		names = dimnames(dataMatrix)
		dataMatrix = normalize.quantiles(dataMatrix)
		dimnames(dataMatrix) = names
	}
	if(gene.level)
		dataMatrix = probesTogene(dataMatrix)

  col.names = colnames(dataMatrix)
	idx = match(comp.group.names, col.names)
	comp.group.names = comp.group.names[!is.na(idx)]
	idx = match(base.group.names, col.names)
	base.group.names = base.group.names[!is.na(idx)]
	subDataMatrix = dataMatrix[,c(comp.group.names, base.group.names)]
 	n.comp.group = length(comp.group.names)
	n.base.group = length(base.group.names)

  
  library(gtools)
	library(qvalue)
	
	pvfoldchange = matrix(0, nrow(subDataMatrix), 5)
	colnames(pvfoldchange) = c('pvalue','qvalue', 'Rawfoldchange','comparison.group.mean', 'base.group.mean')
	rownames(pvfoldchange) = rownames(subDataMatrix)
	row.names = rownames(subDataMatrix)		
	for(k in 1:nrow(subDataMatrix)){
		if(k %% 5000 == 1)
			cat(k, ' ')
		comp.group = subDataMatrix[k,1:n.comp.group]
		base.group = subDataMatrix[k,(n.comp.group+1):ncol(subDataMatrix)]
  
    res.ttest = t.test(comp.group, base.group, paired = paired.ttest)
 		if(is.nan(res.ttest$p.value)){
			pvfoldchange[k,'pvalue'] = 1
			next
		}
		else	
			pvfoldchange[k,'pvalue'] = res.ttest$p.value

		mean.comp.group = mean(comp.group)
		mean.base.group = mean(base.group)
		if(mean.comp.group > mean.base.group)
			pvfoldchange[k,'Rawfoldchange'] = 2^(mean.comp.group - mean.base.group)
		else
			pvfoldchange[k,'Rawfoldchange'] = -2^(mean.base.group - mean.comp.group)
		
#		pvfoldchange[k,'log2foldchange'] = foldchange(mean(comp.group), mean(base.group))
		pvfoldchange[k,'comparison.group.mean'] = mean.comp.group
		pvfoldchange[k,'base.group.mean'] = mean.base.group
	
	}
	cat('\n')
	pvfoldchange[,'qvalue'] = qvalue(p=pvfoldchange[,'pvalue'])$qvalues
	
	write.csv(pvfoldchange, file = file.name)

}
#this returns average probe intenisteis of the same-gene elements
'probesTogene' <- function(dataMatrix.probes){
#	browser()
	row.names = rownames(dataMatrix.probes)
	unique.genes = unique(row.names)
	dataMatrix.genes = matrix(0, length(unique.genes), ncol(dataMatrix.probes))
	colnames(dataMatrix.genes) = colnames(dataMatrix.probes)
	rownames(dataMatrix.genes) = unique.genes
	for(i in 1:length(unique.genes)){
		if(i %% 1000 == 0)
			cat(i, ' ')
		cur.gene.probes = dataMatrix.probes[row.names == unique.genes[i], ,drop=FALSE]
		dataMatrix.genes[i,] = colMeans(cur.gene.probes)
	}
	cat('\n')
	return (dataMatrix.genes)
}
#myVocalnoPlot('Recalcitrant_VS._Healing_Visit_First.csv',log2fc.thr=1)
'myVocalnoPlot' <- function(file.name, log2fc.thr = 0.5, pvalue.thr = 0.05, xmin = -5, xmax = 5, ymin = 0, ymax = 10){
#	browser()
	table = as.matrix(read.csv(file.name, check.names = TRUE))
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
#	pv.fc.list$threshold = as.factor(abs(pv.fc.list$log2fc) > log2fc.thr & pv.fc.list$p.value < pvalue.thr)
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
	write.csv(mlist, file = gsub('.csv', '.sig.list.csv',file.name))
 
##Construct the plot object
	
    dev = pdf(file = gsub('.csv', '.pdf', file.name))
	g<-ggplot(data=pv.fc.list, aes(x=log2fc, y=-log10(p.value), colour=threshold))+
  geom_point()+
  theme(legend.position = "none")+
  xlim(c(xmin, xmax)) + ylim(c(ymin, ymax))+
  xlab("log2 fold change") + ylab("-log10 p-value")
  g <- g + colScale
  #  browser()
#g + geom_text(aes(x=log2fc, y=-log10(p.value), label=ID, size=1.2), colour="black")
    dev.off()
#	boxplot(main="All Samples' boxplot",AVGSignal.annot.noControl[,clinical.info[,1]], outline = FALSE, las =2, xpd=NA, xaxt='n',ylim=c(2,10) );
#	text(x=1:ncol(AVGSignal.annot.noControl),y=1.0,labels=colnames(AVGSignal.annot.noControl[,clinical.info[,1]]),srt=90,cex=0.8,xpd=NA);
#	text(x=1:ncol(AVGSignal.annot.noControl),y=9,labels=labels,srt=90,cex=0.8,xpd=NA);dev.off()
	

}

#sigGenes.Heatmap('Recalcitrant_VS._Healing.sig.list.csv',dataMatrix,comp.group.names=c(samples1$Recalcitrant,samples2$Recalcitrant),base.group.names=c(samples1$Healing,samples2$Healing))
'sigGenes.Heatmap' <- function(sig.list.file.name, dataMatrix, 
																comp.group.names, base.group.names) {													
#	browser()
	library(gplots)
	sig.genes.list = as.matrix(read.csv(file = sig.list.file.name, check.names = FALSE))[,1]
#	dataMatrix = dataMatrixWithCategory$dataMatrix
	sigGenes.dataMatrix = dataMatrix[sig.genes.list,c(comp.group.names, base.group.names)]
	colors.col = c(rep('red', length(comp.group.names)), rep('blue',length(base.group.names)))
	library(som)
	names = dimnames(sigGenes.dataMatrix)
	normalize(sigGenes.dataMatrix)
	dimnames(sigGenes.dataMatrix) = names
	dev = pdf(file = gsub('.sig.list.csv', '_Clustering.pdf', sig.list.file.name))
	
	heatmap.2(sigGenes.dataMatrix,ColSideColors=colors.col,col='greenred',dendrogram='both', scale='row', trace='none',
	main=gsub('.genes.sig.list.csv', ' Clustering', sig.list.file.name), labRow = rownames(sigGenes.dataMatrix),cexRow = 0.9)
	
	dev.off()

}

