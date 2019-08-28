#Copyright Government of Canada 2019
#
#"Written by: Julie Chih-yu Chen @ National Microbiology Laboratory, Public Health Agency of Canada"
#
#"Licensed under the Apache License, Version 2.0 (the ""License""); you may not use"
#this work except in compliance with the License. You may obtain a copy of the
#License at:
#
#http://www.apache.org/licenses/LICENSE-2.0
#
#"Unless required by applicable law or agreed to in writing, software distributed"
#"under the License is distributed on an ""AS IS"" BASIS, WITHOUT WARRANTIES OR"
#"CONDITIONS OF ANY KIND, either express or implied. See the License for the"
#specific language governing permissions and limitations under the License.


#######################################################
####
#### proteomics data: from maxQuant + CARD search
#### 
#######################################################

require(reshape2)
require(gplots)
require(limma)
require(ggplot2)

### load files
maxQuantDir="MaxQuant"
mfile=paste(maxQuantDir,"proteinGroups.txt",sep="\\")
expFile=paste(maxQuantDir,"ExperimentalSetup.txt",sep="\\")
eviFile=paste(maxQuantDir,"evidence.txt",sep="\\")
maxQuantDir=paste0(maxQuantDir,"/")

### param & color definition
perspective="expression" ;datType="labeled"
cbbPalette <- c( "#E69F00", "#56B4E9", "#CC79A7", "#009E73", "#0072B2", "#D55E00","#000000") 
heatmapcol=colorRampPalette(c("white","yellow","darkgoldenrod1","red"),bias=1.5)(50)
heatmapcolVal=seq(2,17.5,length=51)
	

### expSetup where "Name" is the name from maxQuant, groupid is the group name
expSetup<-read.table(expFile, sep="\t", header=T, quote = "", as.is=T)
expSetup$Name=gsub(" ","\\.",expSetup$Name)
rownames(expSetup)=expSetup$Name

colPrefix="Reporter.intensity.corrected\\."
modname=""
fileNameAdd=paste(datType,perspective,modname,sep="_")

modTable<-read.table(mfile, sep="\t", header=T, quote = "", as.is=T)

### getting column ids & sample names
sampDatCols<- unique(sort(unlist(lapply(expSetup$Name,function(y){ grep(y,colnames(modTable))}))))
intenColID<-intersect(grep(colPrefix,colnames(modTable)),sampDatCols)
sampleNames<-gsub(colPrefix,"",colnames(modTable)[intenColID])
pheno=expSetup[sampleNames,]


### getting revData for distribution screening
revDat<-modTable[modTable$Reverse=="+",]  ## saving revData
revDats<- data.matrix(revDat[,intenColID])

### setting the protein presence threshold to 95 percentile of abundance for reversed proteins
revThreshPerc=0.95
revThresh=quantile(as.vector(revDats), probs=revThreshPerc)
print(revThresh)

### filtering for contaminants and reverse
print("Filtering for contaminants and reverse")
modTable<-jFiltContRev(modTable )

### examine distribution of reversed protein abundance versus real protein abundance
pdf(paste0(maxQuantDir,"Density_distribution_reverseVSproteins.pdf"), width=6, height=6)
	idat<-data.frame(cond=c(rep("Reversed",length(as.vector(revDats))), rep("Protein", length(as.vector(data.matrix(modTable[,intenColID]))))), log2Abundance=c(log2(as.vector(revDats)), log2(as.vector(data.matrix(modTable[,intenColID])))))
	# Density plots with semi-transparent fill
	gp<-ggplot(idat, aes(x=log2Abundance, fill=cond)) + geom_density(alpha=.3) +geom_vline(data=idat, aes(xintercept=log2(revThresh)), linetype="dashed", size=1, colour="red")+ labs(x = "Log2 Abundance", y = "Density")+ theme(legend.title = element_blank()) 
	print(gp)
	
dev.off()



##transform to modDatList format
## adding presence/absence matrix datP
datP=data.matrix(modTable[,intenColID])>revThresh 
colnames(datP)=paste0(colnames(datP),"Present")
modDatList<-list(rannot=modTable[,-intenColID], cannot=expSetup[sampleNames,],dat=data.matrix(modTable[,intenColID]), datP=datP)
colnames(modDatList$dat)=gsub(colPrefix,"",colnames(modDatList$dat))



### log2 transformation and filter
modDatList<-jlog2(modDatList)
modDatList.s<-jFilterValid(modDatList,prop=1)  

modDatList.s$dat0<-modDatList.s$dat

## Within batch/run normalization for labeled data
## Adjust total intensity between runs
## watch out if batches are from completely different conditions
## normalize batches, col and row 3 times
for(a in 1:3){
	modDatList.s$dat<-jMedianNormBatch(modDatList.s$dat, batch=modDatList.s$cannot$Run)
	modDatList.s$dat<-jMedianNormBatchRow(modDatList.s$dat, batch=modDatList.s$cannot$Run)
}
	
	
### additional loess normalization, skip if the number of proteins are low due to assumption
if(nrow(modDatList.s$dat)>=100){
	normalization2use="loess"
	normeddat<-jNormalizeLoess(dat4stat=modDatList.s$dat, dat4statFull=data.frame(modDatList.s$rannot,modDatList.s$dat,modDatList.s$datP))
	dat4stat<-normeddat[[1]]; dat4statFull<-normeddat[[2]]; rm(normeddat)
}else{
	normalization2use=""
	dat4stat=modDatList.s$dat
	dat4statFull=data.frame(modDatList.s$rannot,modDatList.s$dat,modDatList.s$datP)
}

	
###  for AMR proteins 
 amrs<-sapply(strsplit(dat4statFull$Majority.protein.IDs,";"), function(y){
	length(y)!=length(grep("gi|^O[0-9]|^P[0-9]|^Q[0-9]",y))
 })


##################
###### Protein abundance
##################
amrDat<-dat4statFull[which(amrs),]
amrDat$Majority.protein.IDs<-gsub("tetO;gi\\|756095536\\|gb\\|AJK83773.1\\|", "gi\\|756095536\\|gb\\|AJK83773.1\\;tetO",gsub("Q03470","gyrA",gsub("CampylobacterjejunigyrAconferringresistancetofluoroquinolones","gyrA mutant",amrDat$Majority.protein.IDs)))

### clean ID
amrDat$Majority.protein.IDs2<-sapply(strsplit(amrDat$Majority.protein.IDs,";"), function(y)y[length(y)])
rownames(amrDat)=amrDat$Majority.protein.IDs2

### presence absence
amrDatPresent<-amrDat[,grep("Present",colnames(amrDat))]
colnames(amrDatPresent)<-gsub("Reporter\\.intensity\\.corrected\\.|Present","",colnames(amrDatPresent))

amrDat.s<-data.matrix(amrDat[,colnames(amrDat)%in%c(pheno$Name, paste0("X",pheno$Name))])
colnames(amrDat.s)=gsub("^X","",colnames(amrDat.s))
amrDat.s<-amrDat.s[, order(pheno$groupid)]


#reordering columns
amrDat.s2<-amrDat.s[,c(13:24,1:12)] ## 6 rep, LM
snames=gsub("i","",pheno[colnames(amrDat.s2),"groupid"])



########################
### same order as RGI
#####################
colnames(rgih$carpet)

amrDat.s3<-merge(data.frame(id=colnames(rgih$carpet),order=1:ncol(rgih$carpet)), data.frame(id=rownames(amrDat.s2),amrDat.s2,amrDat), all=T)
#amrDat.s3<-merge(data.frame(id= gsub(" mutant","",colnames(rgih$carpet)),order=1:ncol(rgih$carpet)), data.frame(id=rownames(amrDat.s2),amrDat.s2), all=T)
amrDat.s3<-amrDat.s3[order(amrDat.s3$order),]
colnames(amrDat.s3)=gsub("X|.1$","",colnames(amrDat.s3))
amrDat.s3<-amrDat.s3[nrow(amrDat.s3):1,]
amrDat.s3Annot<-amrDat.s3[,c( grep("Protein.IDs",colnames(amrDat.s3)):ncol(amrDat.s3))]

amrDat.s4<-amrDat.s3[,-c(1,2, grep("Protein.IDs",colnames(amrDat.s3)):ncol(amrDat.s3))]
rownames(amrDat.s4)=amrDat.s3[,1]
#amrDat.s4[is.na(amrDat.s4)]=0
amrDat.s4<-data.matrix(amrDat.s4)
amrDat.s40<-amrDat.s4;amrDat.s40[is.na(amrDat.s40)]=0

### remove those with Sequence.coverage.... <20% 
amrDat.s5<-amrDat.s4[is.na(amrDat.s3Annot$Unique...razor.sequence.coverage)|amrDat.s3Annot$Sequence.coverage....>20,]


### anova & Tukeys
aovpListProt<-apply(amrDat.s5, 1, function(y,grpid=pheno[colnames(amrDat.s5),"groupid"]){
	if(length(y)!=sum(is.na(y))){
		tmp<-data.frame(y, grpid=factor(grpid))
		res.aov<-aov(y~grpid, data=tmp)
		p0<-summary(res.aov)[[1]][1,5]
		tuk<-TukeyHSD(res.aov)$grpid
	}else{
		p0=NA; tuk=NA
	}
	list(p0, tuk)
})
names(aovpListProt)=rownames(amrDat.s5)
aovpProt<-sapply(aovpListProt, function(y)y[[1]])
aovpBHProt <-p.adjust(aovpProt, method="BH")
tukProt<-lapply(aovpListProt, function(y)y[[2]])
tukProt[which(aovpProt<=0.05)]

tukProtPadj<-sapply(tukProt[which(aovpProt<=0.05)],function(y){y[,"p adj"]})
tukProtDiff<-sapply(tukProt[which(aovpProt<=0.05)],function(y){y[,"diff"]})

write.table(t(tukProtPadj), file=paste0(maxQuantDir,"Table_proteinAbundance_tukeys_padj.txt"),row.names=T, col.names=T, quote=F, sep="\t")
write.table(t(tukProtDiff), file=paste0(maxQuantDir,"Table_proteinAbundance_tukeys_Diff.txt"),row.names=T, col.names=T, quote=F, sep="\t")


siglab<-ifelse(aovpBHProt<=0.05, "*",ifelse(aovpProt<=0.05, "^","")) ### * for significant (<=0.05) with one way ANOVA after BH adjustment; ^ for raw p<0.05, but not significant after adjustment.
siglab[is.na(siglab)]=""

###testing 1597 against others
pvaluesProt<-apply(amrDat.s5,1, function(y, z=pheno[colnames(amrDat.s5),"groupid"]=="i1597"){
	if(sum(is.na(y))!=length(y)){
	 toret=t.test(y[z],y[!z], alternative="greater")$p.value 
	 }else toret=NA
	 toret
})
 adjpvaluesProt<-p.adjust(pvaluesProt,method="BH")

 
 ###testing 1512 against others
pvaluesProt1512<-apply(amrDat.s5,1, function(y, z=pheno[colnames(amrDat.s5),"groupid"]=="i1512"){
	if(sum(is.na(y))!=length(y)){
	 toret=t.test(y[z],y[!z], alternative="greater")$p.value 
	 }else toret=NA
	 toret
})
 adjpvaluesProt1512<-p.adjust(pvaluesProt1512,method="BH")

 
 
 
#### Fig 1B
amrDatPresent.s5=amrDatPresent[rownames(amrDat.s5), colnames(amrDat.s5)]
pdf(paste0(maxQuantDir,"Fig1B_heatmap_AMR_protein_levels_novariants_size8_reorder_matchRGI_rmvScov20_201807.pdf"), width=8)	
	for ( jj in 1:2){
		roundi= 0
		snames=gsub("i","",pheno[colnames(amrDat.s5),"groupid"])
		heatmap.2(amrDat.s5, distfun=dist.pear, keysize=0.9, key.title="log2 Int.",
			ColSideColors= cbbPalette[factor(pheno[colnames(amrDat.s5),"groupid"])], 
			trace="none", scale="none", 
			labRow=paste(rownames(amrDat.s5), siglab),
			labCol=paste(ifelse(snames=="1512","01-","00-"),snames, sep=""), margins=c(7,9), 
			cexCol=1.2,cexRow=1.2,
			cellnote=round(amrDat.s5,roundi),
			  notecex=1.0,
			  notecol="black",
			 sepwidth=c(0.1,0.05),sepcolor=c("white"),colsep=c(6,12,18),
			main="Protein levels matching RGI order", col= heatmapcol,
			breaks=heatmapcolVal,
		 dendrogram = "none", Colv = FALSE, Rowv = FALSE)

		heatmap.2(amrDat.s5, distfun=dist.pear, keysize=0.9, key.title="log2 Int.",
				ColSideColors= cbbPalette[factor(pheno[colnames(amrDat.s5),"groupid"])], 
				trace="none", scale="none", 
				labRow=paste(rownames(amrDat.s5), siglab),
				labCol=paste(ifelse(snames=="1512","01-","00-"),snames, sep=""), margins=c(7,9), 
				cexCol=1.2,cexRow=1.2,
				 sepwidth=c(0.1,0.05),sepcolor=c("white"),colsep=c(6,12,18),
				main="Protein levels matching RGI order", col= heatmapcol,
			 dendrogram = "none", Colv = FALSE, Rowv = FALSE)
			 
		### This is set so when jj==2, it plots only abundance defined to be present according to the 95 percentile reversed protein abundance.
		amrDat.s5[!(amrDatPresent.s5)]=NA	
		heatmapcolVal=seq(12, 17.5,length=51)
	}
dev.off()



#############################
#### peptide analysis
#############################


evi<-read.table(eviFile, sep="\t", header=T, quote = "", as.is=T)
evidat<-evi[,grep("Reporter.intensity.corrected",colnames(evi))]
expNames<- unique(evi$Experiment) # this is the number of replicates

### getting revData for distribution screening
revDatP<-evi[evi$Reverse=="+"&evi$Potential.contaminant!="+",]  ## saving revData
revDatPs<- data.matrix(evidat[evi$Reverse=="+"&evi$Potential.contaminant!="+",])
#revDatPs<-revDatPs[!duplicated(revDatP$Sequence),]


datColNames=c("Reporter.intensity.corrected.0", "Reporter.intensity.corrected.1", "Reporter.intensity.corrected.2", "Reporter.intensity.corrected.3")
eviMedInt<-apply(evi[,datColNames], 1, median, na.rm=T)


### getting median abundance in each arrach for each reverse peptide and obtain the 95th percentile abundance per peptide length
revThreshPercP=0.95
tmp2<-aggregate(eviMedInt[evi$Reverse=="+"], list(evi$Length[evi$Reverse=="+"]),quantile,probs=revThreshPercP)
tmp2<-cbind(tmp2,table(evi$Length[evi$Reverse=="+"]),log2abundance=log2(1+tmp2[,2]))
	
		
### fitting abundance loess line for reversed peptides
lo <- loess(x~Group.1, tmp2 ,span = 1, degree = 2); 
tmp2pred <- predict(lo, data.frame( Group.1=as.numeric(names(table(evi$Length[evi$Reverse!="+"])))))
names(tmp2pred)=as.numeric(names(table(evi$Length[evi$Reverse!="+"])))
(tmp2predl2<-log2(1+tmp2pred))  ### log2 abundance thresholds (from reversed peptides) for varying peptide length
tmp2predl2[c(min(which(is.nan(tmp2predl2))): length(tmp2predl2))]=0  ### making the rest of the lengths 0

### visualizing the abundance thresholds determined w.r.t. peptide lengths & highlights higher reported abundance when peptide lengths are shorter
pdf(paste0(maxQuantDir,"Scatter_peptide_scores_deltascores_reverse_vs_hits_Abundance.pdf"), width=6,height=6)
	plot(jitter(evi$Length[evi$Reverse!="+"], 2),eviMedInt[evi$Reverse!="+" ], col="#F8766D",  pch=".", ylim=range(eviMedInt, na.rm=T), xlab="Peptide length", ylab="Median abundance of each peptide in an array")
	points(jitter(evi$Length[evi$Reverse=="+"], 2),eviMedInt[evi$Reverse=="+"], col="skyblue3",pch=20)
	lines(names(tmp2pred),tmp2pred, lwd=3,col="darkblue")
	legend("topright", legend=c("Peptide", "Reversed"), col=c("#F8766D","skyblue3") , pch=20)
dev.off()



protDatP<-evi[evi$Reverse!="+" | evi$Potential.contaminant!="+",]  ## peptide data with reversed and contaminant filtered out.


### the few peptides only
targetPep<-c("IMAIIPTTDFDESK","YHPHGDTAVYDALVR","YHPHGDIAVYDALVR","IALDNIDEVIALIK","IALDNIDK") 
evi.s<-protDatP[protDatP$Sequence%in% targetPep,]
pepl2i<-t(sapply(unique(evi.s$Sequence), toPeptideTable, dat0=evi.s, expNames=expNames))

pepl2i<-pepl2i[,order((colnames(pepl2i)))]

pepColOrder=c(7:18,1:6,19:24)
pepl2i<-pepl2i[targetPep,pepColOrder]#c("i1597","i6200","i0949","i1512")) ##  6 rep


#### function to determine absence according to loess fitted 95 percentile
jScreenAbs<-function(datMat, lengthThresholds){	
	dimnames0 <- dimnames(datMat)
	datMat <- cbind(thresh=lengthThresholds[as.character(nchar(rownames(datMat)))], datMat)
	ret<-t(apply(datMat,1, function(y){
		y[y< y[1]]=NA
		y[-1]
	}))
	dimnames(ret)=dimnames0
	ret
}
pepl2iAbs <- jScreenAbs(pepl2i, lengthThresholds=tmp2predl2)

pepsname<-gsub("i","",pheno[colnames(pepl2i),"groupid"])

pdf(paste0(maxQuantDir,"Fig2_heatmap_AMR_peptide_levels_new_201807_scalematch.pdf"), width=8, height=4)	
	roundi= 0
		heatmap.2(pepl2i, 
		ColSideColors= cbbPalette[factor(pheno[colnames(pepl2i),"groupid"])], 
		trace="none", scale="none", 
		labRow=rownames(pepl2i), 
		labCol=paste0(ifelse(pepsname=="1512","01-","00-"),pepsname), margins=c(5,10), 
		  cellnote=round(pepl2i,roundi),
			  notecex=1.0,
			  notecol="black",
		main="Peptide log2 intensity: 3 examples for gyrA", col= heatmapcol,
		#breaks=heatmapcolVal,
		cexCol=1,cexRow=1,
		sepwidth=c(0.1,0.05),sepcolor=c("white"),colsep=c(6,12,18),rowsep=c(1,3),
		dendrogram = "none", Colv = FALSE, Rowv = FALSE, 
		lhei = c(1.2,1.8))
	
		### with absence filter
		heatmap.2(pepl2iAbs, 
		ColSideColors= cbbPalette[factor(pheno[colnames(pepl2iAbs),"groupid"])], 
		trace="none", scale="none", 
		labRow=rownames(pepl2iAbs), 
		labCol=paste0(ifelse(pepsname=="1512","01-","00-"),pepsname), margins=c(5,10), 
		  cellnote=round(pepl2iAbs,roundi),
			  notecex=1.0,
			  notecol="black",
		main="Peptide log2 intensity: 3 examples for gyrA", col= heatmapcol,
		#breaks=heatmapcolVal,
		cexCol=1,cexRow=1,
		sepwidth=c(0.1,0.05),sepcolor=c("white"),colsep=c(6,12,18),rowsep=c(1,3),
		dendrogram = "none", Colv = FALSE, Rowv = FALSE, 
		lhei = c(1.2,1.8))
	

dev.off()

	
			
# ANOVA
aovpList<-apply(pepl2i, 1, function(y,grpid=pheno[colnames(pepl2i),"groupid"]){
	tmp<-data.frame(y, grpid=factor(grpid))
	res.aov<-aov(y~grpid, data=tmp)
	p0<-summary(res.aov)[[1]][1,5]
	tuk<-TukeyHSD(res.aov)$grpid
	list(p0, tuk)
})
aovp<-sapply(aovpList, function(y)y[[1]])
aovpBH <-p.adjust(aovp, method="BH")
tuk<-lapply(aovpList, function(y)y[[2]])
tuk[aovp<=0.05]


### two sided, compare i1597 to others
id1597 <- pheno[colnames(pepl2i),"groupid"]=="i1597"
for (pepname in targetPep){
	cat(pepname, " ")
	cat(t.test(pepl2i[pepname,id1597],pepl2i[pepname,!id1597])$p.value , "\n")
}
#YHPHGDIAVYDALVR  0.002563197 
#YHPHGDTAVYDALVR  7.444988e-05 
#IALDNIDK  0.001976885 
#IALDNIDEVIALIK  0.05284511 
#IMAIIPTTDFDESK  0.4801575 

### testing 1597 against others, one sided
pvalues<-apply(pepl2i,1, function(y, z=id1597){
	t.test(y[z],y[!z], alternative="greater")$p.value 

})
pvalues
# IMAIIPTTDFDESK YHPHGDTAVYDALVR YHPHGDIAVYDALVR  IALDNIDEVIALIK        IALDNIDK 
#  0.7599212596    0.9999627751    0.0012815985    0.9735774448    0.0009884423
adjpvalues<-p.adjust(pvalues,method="BH")



2^( mean(pepl2i["YHPHGDIAVYDALVR",id1597], na.rm=T)- mean(pepl2i["YHPHGDIAVYDALVR",!id1597], na.rm=T)) # 12.69232
2^( mean(pepl2i["IALDNIDK",id1597], na.rm=T)- mean(pepl2i["IALDNIDK",!id1597], na.rm=T)) #13.03566


