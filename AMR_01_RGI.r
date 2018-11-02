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

#############################################################
### new rgi 4.1.0 2018
########################################################
library(plyr)
library(reshape2)
library(gplots)

rgiFolder="RGI"
jsons<-list.files(rgiFolder,pattern=".fasta.txt|.json.txt|rgiSummary_|coding.txt.txt$")
jsonsNames<-gsub(".fasta.txt|.json.txt|.fasta.json.txt|rgiSummary_|.txt|_coding","",jsons)
jsons<-paste(rgiFolder,jsons,sep="\\")
cbbPalette <- c( "#E69F00", "#56B4E9", "#CC79A7", "#009E73", "#0072B2", "#D55E00","#000000") #, "#F0E442" yellow

	

#### read jsons
rgi<-lapply(jsons, function(json){
	print(".")
	read.table(json, sep="\t", header=T, as.is=T,quote="", comment.char="") #, comment.char="" as the new output from CARD includes '#' character
})
names(rgi)=jsonsNames
	#rgi$'01-1512'<-rbind(rgi$'01-1512',rgi$'01-1512_plasmid_pCj1')
	#rgi$'00-0949'<-rbind(rgi$'00-0949',rgi$'00-0949_plasmid_pTet')
	#rgi<-rgi[!names(rgi)%in%c('01-1512_plasmid_pCj1','00-0949_plasmid_pTet')]

#### to data frame, long format
rgidf<- ldply(rgi, data.frame)

#### reshape
bestID<-dcast(rgidf, .id~Best_Hit_ARO, value.var = "Best_Identities", fun.aggregate = mean)
bestID[is.na(bestID)]=0
rownames(bestID)=bestID[,1]; bestID<-data.matrix(bestID[,-1])
bestIDcolnamesO<-colnames(bestID)


colnames(bestID)=gsub("Campylobacter jejuni gyrA conferring resistance to fluoroquinolones","gyrA mutant",colnames(bestID))
		#,gsub("Mycobacterium tuberculosis gyrB mutant conferring resistance to fluoroquinolone","gyrB mutant",
		#gsub("Staphylococcus aureus rpoB mutants conferring resistance to rifampicin","rpoB mutant",
		#gsub("Streptomyces cinnamoneus EF-Tu mutants conferring resistance to elfamycin","EF-Tu mutant",
		#gsub("antibiotic resistant fabI","fabI",
		#gsub("with mutation conferring resistance to|with mutation conferring antibiotic resistance|mutant conferring antibiotic resistance", "mutant", gsub("Escherichia coli","E. coli",colnames(bestID))))))))
bestID<-t(bestID)

bestE<-dcast(rgidf, .id~Best_Hit_ARO, value.var = "Best_Hit_Bitscore", fun.aggregate = mean)
bestE[is.na(bestE)]=0
rownames(bestE)=bestE[,1]; bestE<-data.matrix(bestE[,-1])
colnames(bestE)=rownames(bestID)
bestE<-t(bestE)
	 	 
bestCate<-t(dcast(rgidf, .id~Best_Hit_ARO, value.var = "Cut_Off", fun.aggregate = paste,collapse="_"))
colnames(bestCate)=bestCate[1,]
bestCate<-bestCate[-1,]
bestCate=gsub("Strict","*",gsub("Perfect","**",bestCate))

	 
pdf("Fig1A_Heatmap_summary_AMR_CARD_RGI_201806.pdf", width=6) ## margin changed from 5 to 7 with 00-, 01- for campy
	rgih<-heatmap.2(bestID, scale="none",trace="none", margin=c(7,9), ColSideColors= cbbPalette[factor(colnames(bestID), levels=c("00-0949","01-1512","00-1597","00-6200"))], 
		col=colorRampPalette(c("white","lightpink","pink","red"),bias=0.1)(30), keysize=0.9,key.title="% Ident.", cexRow=1.2, cexCol=1.2,
		, cellnote=matrix( paste0(round(bestID,1), bestCate),ncol=ncol(bestID)), notecex=1.2, notecol="black")	
dev.off()




#### Output a table
rgiSumTable<-merge(data.frame(Best_Hit_ARO=rownames(bestCate), bestID, stringsAsFactors=F), rgidf) ## changed from rownames(bestID) to rownames(bestCate) as bestID row names has been changed.
### the following info can differ between samples
###.id ORF_ID Contig Start Stop Orientation Cut_Off Pass_Bitscore Best_Hit_Bitscore Best_Identities
rgiSumTable.s<-rgiSumTable[!duplicated(rgiSumTable$Best_Hit_ARO),]
rgiSumTable2<-sapply(unique(rgiSumTable$Best_Hit_ARO), function(y){
	ytmp<-rgiSumTable[y==rgiSumTable$Best_Hit_ARO,] #c(".id","ORF_ID")
	origGeneID<-sapply(strsplit(ytmp[,"ORF_ID"]," "), function(z)z[1])
	for(x in 1:nrow(ytmp)){
		origGeneID[[x]]=paste(paste(ytmp$.id[x],origGeneID[[x]],sep=","), collapse=";")
	}
	origGeneID
})
rgiSumTable.s<-data.frame(rgiSumTable.s, sampIDgeneID=sapply(rgiSumTable2,paste,collapse=";"), stringsAsFactors=F)
write.table(rgiSumTable.s, file="Table_summary_AMR_CARD_RGI.txt",row.names=F, col.names=T, quote=F, sep="\t")

