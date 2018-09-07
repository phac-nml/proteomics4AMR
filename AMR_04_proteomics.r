#Copyright Government of Canada 2018
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
#### proteomics data: from maxQuant +CARD search
#### 
#######################################################

require(reshape2)
require(gplots)
require(limma)


maxQuantDir="MaxQuant"
mfile=paste(maxQuantDir,"proteinGroups.txt",sep="\\")
expFile=paste(maxQuantDir,"ExperimentalSetup.txt",sep="\\")
eviFile=paste(maxQuantDir,"evidence.txt",sep="\\")


perspective="expression" ;datType="labeled"
cbbPalette <- c( "#E69F00", "#56B4E9", "#CC79A7", "#009E73", "#0072B2", "#D55E00","#000000") #, "#F0E442" yellow



heatmapcol=colorRampPalette(c("white","yellow","darkgoldenrod1","red"),bias=1.5)(50)
heatmapcolVal=seq(2,17.5,length=51)
	

### expSetup where "Name" is the name from maxQuant, groupid is the group name
expSetup<-read.table(expFile, sep="\t", header=T, quote = "", as.is=T)
expSetup$Name=gsub(" ","\\.",expSetup$Name)
rownames(expSetup)=expSetup$Name

colPrefix="Reporter.intensity.corrected\\."
modname=""
modTable<-read.table(mfile, sep="\t", header=T, quote = "", as.is=T)


### filtering for contaminants and reverse
print("Filtering for contaminants and reverse")
modTable<-jFiltContRev(modTable )


sampDatCols<- unique(sort(unlist(lapply(expSetup$Name,function(y){ grep(y,colnames(modTable))}))))
intenColID<-intersect(grep(colPrefix,colnames(modTable)),sampDatCols)
sampleNames<-gsub(colPrefix,"",colnames(modTable)[intenColID])

##transform to modDatList format
modDatList<-list(rannot=modTable[,-intenColID], cannot=expSetup[sampleNames,],dat=data.matrix(modTable[,intenColID]) )
colnames(modDatList$dat)=gsub(colPrefix,"",colnames(modDatList$dat))


boxplot(modDatList$dat) 
colSums(modDatList$dat)

### log2 a
modDatList<-jlog2(modDatList)
modDatList.s<-jFilterValid(modDatList,prop=1)  


## Within batch/run normalization for labeled data
## Adjust total intensity between runs
modDatList.s$dat0<-modDatList.s$dat

## watch out if batches are from completely different conditions
## normalize batches, col and row 3 times
for(a in 1:3){
	modDatList.s$dat<-jMedianNormBatch(modDatList.s$dat, batch=modDatList.s$cannot$Run)
	modDatList.s$dat<-jMedianNormBatchRow(modDatList.s$dat, batch=modDatList.s$cannot$Run)
}
	
### boxplot
boxplot(modDatList.s$dat,  main="Batch normed") 
pheno=modDatList.s$cannot
fileNameAdd=paste(datType,perspective,modname,sep="_")

### additional loess normalization, skip if the number of proteins are low!
if(nrow(modDatList.s$dat)>=100){
	normalization2use="loess"
	normeddat<-jNormalizeLoess(dat4stat=modDatList.s$dat, dat4statFull=data.frame(modDatList.s$rannot,modDatList.s$dat))
	dat4stat<-normeddat[[1]]; dat4statFull<-normeddat[[2]]; rm(normeddat)
}else{
	normalization2use=""
	dat4stat=modDatList.s$dat
	dat4statFull=data.frame(modDatList.s$rannot,modDatList.s$dat)
}

		
		
###  for AMR!

 amrs<-sapply(strsplit(dat4statFull$Majority.protein.IDs,";"), function(y){
	length(y)!=length(grep("gi|^O[0-9]|^P[0-9]|^Q[0-9]",y))
 })

 
amrDat<-dat4statFull[amrs,]
amrDat$Majority.protein.IDs2<-gsub("CampylobacterjejunigyrAconferringresistancetofluoroquinolones","gyrA mutant",
	gsub("EnterococcusfaeciumEF-TumutantsconferringresistancetoGE2270A","EF-Tu mutant",
	gsub("MycobacteriumtuberculosisrpsLmutationsconferringresistancetoStreptomycin","rpsL mutant",
	gsub("PlanobisporaroseaEF-TumutantsconferringresistancetoinhibitorGE2270A;StreptomycescinnamoneusEF-Tumutantsconferringresistancetoelfamycin","EF-Tu mutant",
	gsub("StaphylococcusaureusrpoCconferringresistancetodaptomycin","rpoC mutant",
	gsub("HaemophilusparainfluenzaegyrAconferringresistancetofluoroquinolones","gyrA mutant Haemophiluspara",
	gsub("StaphylococcusaureusmurAwithmutationconferringresistancetofosfomycin","murA mutant",
	gsub("tetO;gi\\|756095536\\|gb\\|AJK83773.1\\|","tetO",
	amrDat$Majority.protein.IDs))))))))
amrDat$Majority.protein.IDs2<-sapply(strsplit(amrDat$Majority.protein.IDs2,";"), function(y)y[length(y)])

amrDat.s<-data.matrix(amrDat[,colnames(amrDat)%in%c(pheno$Name, paste0("X",pheno$Name))])
colnames(amrDat.s)=gsub("^X","",colnames(amrDat.s))



##################
###### Protein abundance
##################
amrDat<-dat4statFull[c(which(amrs& nchar(dat4statFull$Majority.protein.IDs)<40)),]
amrDat$Majority.protein.IDs<-gsub("tetO;gi\\|756095536\\|gb\\|AJK83773.1\\|", "gi\\|756095536\\|gb\\|AJK83773.1\\;tetO",gsub("Q03470","gyrA",gsub("CampylobacterjejunigyrAconferringresistancetofluoroquinolones","gyrA mutant",amrDat$Majority.protein.IDs)))

amrDat$Majority.protein.IDs2<-sapply(strsplit(amrDat$Majority.protein.IDs,";"), function(y)y[length(y)])
rownames(amrDat)=amrDat$Majority.protein.IDs2


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
amrDat.s3=amrDat.s3[nrow(amrDat.s3):1,]
amrDat.s3Annot<-amrDat.s3[,c( grep("Protein.IDs",colnames(amrDat.s3)):ncol(amrDat.s3))]

amrDat.s4<-amrDat.s3[,-c(1,2, grep("Protein.IDs",colnames(amrDat.s3)):ncol(amrDat.s3))]
rownames(amrDat.s4)=amrDat.s3[,1]
#amrDat.s4[is.na(amrDat.s4)]=0
amrDat.s4<-data.matrix(amrDat.s4)
amrDat.s40<-amrDat.s4;amrDat.s40[is.na(amrDat.s40)]=0


### using! remove those with Sequence.coverage.... <20% 
amrDat.s5<-amrDat.s4[is.na(amrDat.s3Annot$Unique...razor.sequence.coverage)|amrDat.s3Annot$Sequence.coverage....>20,]


mean(amrDat.s5["tetO",13:24])-mean(amrDat.s5["tetO",1:12])
#[1] 3.149055
t.test(amrDat.s5["tetO",13:24], amrDat.s5["tetO",1:12], alternative="greater")
# 1.265e-12


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

write.table(t(tukProtPadj), file="Table_proteinAbundance_tukeys_padj.txt",row.names=T, col.names=T, quote=F, sep="\t")
write.table(t(tukProtDiff), file="Table_proteinAbundance_tukeys_Diff.txt",row.names=T, col.names=T, quote=F, sep="\t")


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

 
 
 
#### using!
pdf("Fig1B_heatmap_AMR_protein_levels_novariants_size8_reorder_matchRGI_rmvScov20_201807.pdf", width=8)	
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
dev.off()




#############################
#### peptide analysis
#############################



evi<-read.table(eviFile, sep="\t", header=T, quote = "", as.is=T)
evidat<-evi[,grep("Reporter.intensity.corrected",colnames(evi))]
expNames<- unique(evi$Experiment) # this is the number of replicates



### the few peptides only
targetPep<-c("YHPHGDIAVYDALVR","YHPHGDTAVYDALVR","IALDNIDK","IALDNIDEVIALIK","IMAIIPTTDFDESK")
#targetPep<-c("YHPHGDIAVYDALVR","IALDNIDK","IMAIIPTTDFDESK")
evi.s<-evi[evi$Sequence%in% targetPep,]
pepl2i<-t(sapply(unique(evi.s$Sequence), toPeptideTable, dat0=evi.s, expNames=expNames))

#pepl2i<-pepl2i[c(3,1,2),order((colnames(pepl2i)))]
pepl2i<-pepl2i[,order((colnames(pepl2i)))]

#pepColorder=c(4:9,1:3,10:12)#c("i1597","i6200","i0949","i1512")) ##  3 rep only
pepColOrder=c(7:18,1:6,19:24)
pepl2i<-pepl2i[targetPep,pepColOrder]#c("i1597","i6200","i0949","i1512")) ##  6 rep


pepsname<-gsub("i","",pheno[colnames(pepl2i),"groupid"])

pdf("Fig2_heatmap_AMR_peptide_levels_new_201807_scalematch.pdf", width=8, height=4)	
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
		breaks=heatmapcolVal,
		cexCol=1,cexRow=1,
		sepwidth=c(0.1,0.05),sepcolor=c("white"),colsep=c(6,12,18),rowsep=c(2,4),
		dendrogram = "none", Colv = FALSE, Rowv = FALSE, 
		lhei = c(1.2,1.8))
	
	
	 heatmap.2(pepl2i, 
		ColSideColors= cbbPalette[factor(pheno[colnames(pepl2i),"groupid"])], 
		trace="none", scale="none", 
		labRow=rownames(pepl2i), 
		labCol=paste0(ifelse(pepsname=="1512","01-","00-"),pepsname), margins=c(5,10), 
		main="Peptide log2 intensity: 3 examples for gyrA", col= heatmapcol,
		cexCol=1,cexRow=1,
		sepwidth=c(0.1,0.05),
         sepcolor=c("darkgray"),#,
		colsep=c(6,12,18),
		rowsep=c(2,4),
		dendrogram = "none", Colv = FALSE, Rowv = FALSE, 
		lhei = c(1.2,1.8))
		
dev.off()


	
		
### all gyra peptides
gyra<-evi[grep("CampylobacterjejunigyrAconferringresistancetofluoroquinolones|Q03470",evi$Proteins),] #Q03470
unique(gyra$Sequence)
pepl2i<-t(sapply(unique(gyra$Sequence), toPeptideTable, dat0=gyra, expNames=expNames))
#pepl2i<-pepl2i[c(3,1,2),order((colnames(pepl2i)))]
pepl2i<-pepl2i[,order((colnames(pepl2i)))]
pepl2i<-pepl2i[,pepColOrder]#c("i1597","i6200","i0949","i1512"))
pepsname<-gsub("i","",pheno[colnames(pepl2i),"groupid"])



### filter for at least 50% sampels with data
pepl2i2<- pepl2i[rowSums(!is.na(pepl2i))>=(0.5*ncol(pepl2i)),]

# ANOVA
aovpList<-apply(pepl2i2, 1, function(y,grpid=pheno[colnames(pepl2i2),"groupid"]){
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


for (pepname in targetPep){
	cat(pepname, " ")
	cat(t.test(pepl2i2[pepname,1:6],pepl2i2[pepname,7:ncol(pepl2i2)])$p.value , "\n")
}
#YHPHGDIAVYDALVR  0.002563197 
#YHPHGDTAVYDALVR  7.444988e-05 
#IALDNIDK  0.001976885 
#IALDNIDEVIALIK  0.05284511 
#IMAIIPTTDFDESK  0.4801575 

### testing 1597 against others, one sided
pvalues<-apply(pepl2i2,1, function(y, z=pheno[colnames(pepl2i2),"groupid"]=="i1597"){
	t.test(y[z],y[!z], alternative="greater")$p.value 

})
adjpvalues<-p.adjust(pvalues,method="BH")


2^( mean(pepl2i2["YHPHGDIAVYDALVR",1:6], na.rm=T)- mean(pepl2i2["YHPHGDIAVYDALVR",7:24], na.rm=T)) # 12.69232
2^( mean(pepl2i2["IALDNIDK",1:6], na.rm=T)- mean(pepl2i2["IALDNIDK",7:24], na.rm=T)) #13.03566



### sequence to sort by location
gyrA="MENIFSKDSDIELVDIENSIKSSYLDYSMSVIIGRALPDARDGLKPVHRRILYAMQNDEAKSRTDFVKSARIVGAVIGRYHPHGDIAVYDALVRMAQDFSMRYPSITGQGNFGSIDGDSAAAMRYTEAKMSKLSHELLKDIDKDTVDFVPNYDGSESEPDVLPSRVPNLLLNGSSGIAVGMATNIPPHSLNELIDGLLYLLDSKDASLEEIMQFIKGPDFPTGGIIYGKKGIIEAYRTGRGRVKVRAKTHIEKKTNKDVIVIDELPYQTNKARLIEQIAELVKEKQIEGISEVRDESNKEGIRVVIELKREAMSEIVLNNLFKSTTMESTFGVIMLAIYNKEPKIFSLLELLNLFLTHRKTVIIRRTIFELQKARARAHILEGLKIALDNIDKVIALIKNSSDNNTARDSLVAKFGLSELQANAILDMKLGRLTGLEREKIENELAELMKEIARLEEILKSETLLENLIRDELKEIRSKFDVPRITQIEDDYDDIDIEDLIPNENMVVTITHRGYIKRVPSKQYEKQKRGGKGKLAVTTYDDDFIESFFTANTHDTLMFVTDRGQLYWLKVYKIPEGSRTAKGKAVVNLINLQAEEKIMAIIPTTDFDESKSLCFFTKNGIVKRTNLSEYQNIRSVGVRAINLDENDELVTAIIVQRDEDEIFATGGEENLENQEIENLDDENLENEESVSTQGKMLFAVTKKGMCIKFPLAKVREIGRVSRGVTAIKFKEKNDELVGAVVIENDEQEILSISAKGIGKRTNAGEYRLQSRGGKGVICMKLTEKTKDLISVVIVDETMDLMALTSSGKMIRVDMQSIRKAGRNTSGVIVVNVENDEVVSIAKCPKEENDEDELSDENFGLDLQ"
swissp="MENIFSKDSDIELVDIENSIKSSYLDYSMSVIIGRALPDARDGLKPVHRRILYAMQNDEAKSRTDFVKSARIVGAVIGRYHPHGDTAVYDALVRMAQDFSMRYPSITGQGNFGSIDGDSAAAMRYTEAKMSKLSHELLKDIDKDTVDFVPNYDGSESEPDVLPSRVPNLLLNGSSGIAVGMATNIPPHSLNELIDGLLYLLDNKDASLEEIMQFIKGPDFPTGGIIYGKKGIIEAYRTGRGRVKVRAKTHIEKKTNKDVIVIDELPYQTNKARLIEQIAELVKERQIEGISEVRDESNKEGIRVVIELKREAMSEIVLNNLFKSTTMESTFGVIMLAIHNKEPKIFSLLELLNLFLTHRKTVIIRRTIFELQKARARAHILEGLKIALDNIDEVIALIKNSSDNNTARDSLVAKFGLSELQANAILDMKLGRLTGLEREKIENELAELMKEIARLEEILKSETLLENLIRDELKEIRSKFDVPRITQIEDDYDDIDIEDLIPNENMVVTITHRGYIKRVPSKQYEKQKRGGKGKLAVTTYDDDFIESFFTANTHDTLMFVTDRGQLYWLKVYKIPEGSRTAKGKAVVNLINLQAEEKIMAIIPTTDFDESKSLCFFTKNGIVKRTNLSEYQNIRSVGVRAINLDENDELVTAIIVQRDEDEIFATGGEENLENQEIENLDDENLENEESVSTQGKMLFAVTKKGMCIKFPLAKVREIGRVSRGVTAIKFKEKNDELVGAVVIENDEQEILSISAKGIGKRTNAGEYRLQSRGGKGVICMKLTEKTKDLISVVIVDETMDLMALTSSGKMIRVDMQSIRKAGRNTSGVIVVNVENDEVVSIAKCPKEENDEDELSDENFGLDLQ"

aaposs<-t(sapply(rownames(pepl2i2), function(y){
	tmp<-c(nchar(strsplit(gyrA, y)[[1]][1]), nchar(strsplit(swissp, y)[[1]][1]))
	c(ifelse(tmp[1]==nchar(gyrA), NA, tmp[1]), ifelse(tmp[2]==nchar(swissp), NA, tmp[2]))
}))+1
colnames(aaposs)=c("card","swissp")

### both:1 black; CARD only :2 red; swissprot only:4 blue
aawhich<-apply(aaposs,1,function(y){
	if(sum(is.na(y))==0){
		ret=1
	}else{
		ret=which(!is.na(y))+1
		ret=ifelse(ret==3,4,2)
	}
	ret})
	
	
aapos<-apply(aaposs, 1, min, na.rm=T)
aaorder<-order(aapos)

pdf("FigS4_heatmap_AMR_peptide_levels_gyra_filt50perc_201807_sortedByLocation.pdf", width=8, height=6.2)	
roundi= 0
heatmap.2(pepl2i2[aaorder,], distfun=dist.pear, 
	ColSideColors= cbbPalette[factor(pheno[colnames(pepl2i2),"groupid"])], 
	trace="none", scale="none", 
	labRow=paste(rownames(pepl2i2), ifelse(aovpBH<=0.05, "*",ifelse(aovp<=0.05, "^","")))[aaorder],
	colRow=aawhich[aaorder], 
	labCol=paste0(ifelse(pepsname=="1512","01-","00-"),pepsname), margins=c(5,14), 
	  cellnote=round(pepl2i2[aaorder,],roundi),
		  notecex=0.9,
		  notecol="black",
	sepwidth=c(0.2, 0.1),sepcolor=c("white"),colsep=c(6,12,18),rowsep=c(5,7,20,22),
	main="Peptide log2 intensity for gyrA 50%filt", 
	col= heatmapcol,
	breaks=heatmapcolVal,
	cexCol=1,cexRow=1,
	dendrogram = "none", Colv = FALSE, Rowv = FALSE)
dev.off()



###My order
aaposMine=aapos+5
aaposMine[targetPep[1:4]]=1:4
aaorderMine<-order(aaposMine)
pdf("heatmap_AMR_peptide_levels_gyra_filt50perc_201807_sortedByLocation_highlight.pdf", width=8, height=6.2)	
	roundi= 0
	heatmap.2(pepl2i2[aaorderMine,], distfun=dist.pear, 
		ColSideColors= cbbPalette[factor(pheno[colnames(pepl2i2),"groupid"])], 
		trace="none", scale="none", 
		labRow=paste(rownames(pepl2i2), ifelse(aovpBH<=0.05, "*",ifelse(aovp<=0.05, "^","")))[aaorderMine],
		colRow=aawhich[aaorderMine], 
		labCol=paste0(ifelse(pepsname=="1512","01-","00-"),pepsname), margins=c(5,14), 
		  cellnote=round(pepl2i2[aaorderMine,],roundi),
			  notecex=0.9,
			  notecol="black",
		sepwidth=c(0.2, 0.05),sepcolor=c("white"),colsep=c(6,12,18),rowsep=c(2,4),
		main="Peptide log2 intensity for gyrA 50%filt", 
		col= heatmapcol,
		breaks=heatmapcolVal,
		cexCol=1,cexRow=1,
		dendrogram = "none", Colv = FALSE, Rowv = FALSE)
dev.off()



