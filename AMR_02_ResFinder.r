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
### ResFinder
########################################################
library(plyr)
library(reshape2)
library(gplots)
 
rfFolder="ResFinder" 
 
rfgene<-list.files(rfFolder,pattern="_results_tab.txt")
rfpoint<-list.files(rfFolder,pattern="PointFinder_results.txt")
fNames<-gsub("_results_tab.txt","",rfgene)
rfgene<-paste(rfFolder,rfgene,sep="\\")
rfpoint<-paste(rfFolder,rfpoint,sep="\\")

	
	
	
	
#### read rfgene
rfgene2<-lapply(rfgene, function(file0){
	print(".")
	tmp0<-readLines(file0)
	nrs<-which(tmp0=="No resistance genes found.")
	tmp0<-tmp0[-c(nrs,nrs-1, nrs+1)]
	tmp1<-strsplit(tmp0,"\t")
	tmp1<-do.call(rbind,tmp1[sapply(tmp1,length)==7])
	tmp1<-tmp1[!duplicated(tmp1),]
	colnames(tmp1)=tmp1[1,]
	if(!is.null(tmp1)){
		toret<-(tmp1[-1,,drop=F])
	}else{
		toret<-rep(NA, 7)
	}
	toret
})
names(rfgene2)=fNames


#### to data frame, long format
rfgene2<- ldply(rfgene2, data.frame)
rfgene2<-data.frame(Isolate=sapply(strsplit(rfgene2$.id,"_"),function(y)y[1]), rfgene2, stringsAsFactors=F)
rfgene2$Identity <-as.numeric(as.character(rfgene2$Identity))
rfgene2$Isolate=factor(rfgene2$Isolate, levels=c("00-0949","01-1512","00-1597","00-6200")) ## didn't help

rfID<-dcast(rfgene2, Isolate~Resistance.gene, value.var = "Identity", fun.aggregate = mean)
rfID[is.na(rfID)]=0
rownames(rfID)=rfID[,1]; 
rfID<-data.matrix(rfID[,-c(1, which(colnames(rfID)=="NA"))])
rfID[levels(rfgene2$Isolate),]




#### read rfpoint
rfpoint2<-lapply(rfpoint, function(file0){
	suppressWarnings(tmp0<-read.table(file0, sep="\t", header=T, as.is=T, fill=T) )
	if(nrow(tmp0)==0){
		holder=t(rep(NA,ncol(tmp0)))
		colnames(holder)=colnames(tmp0)
		tmp0=rbind(tmp0,holder)
	}
	tmp0
})
names(rfpoint2)=gsub(paste0(rfFolder,"\\\\"),"",gsub("_PointFinder_results.txt","",rfpoint))
rfpoint2<- ldply(rfpoint2, data.frame)
rfpoint2<-data.frame(Isolate=sapply(strsplit(rfpoint2$.id,"_"),function(y)y[1]), rfpoint2)
rfpoint2

rfpoint2<-dcast(rfpoint2, Isolate~Mutation, value.var = "Amino.acid.change", fun.aggregate = length)
rfpoint2<-rfpoint2[,colnames(rfpoint2)!="NA"]
rownames(rfpoint2)=rfpoint2[,1]; rfpoint2<-as.matrix(rfpoint2[,-1, drop=F])
rfpoint2[rfpoint2==1]=100

 rfpoint2<- rfpoint2[rownames(rfID),, drop=F]
rfResult<-cbind(rfID,rfpoint2)
colnames(rfResult)=c(colnames(rfID), colnames(rfpoint2))
rfResult<-t(rfResult)

pdf("FigS2_Heatmap_summary_ResFinder_201806_celllabel.pdf", width=9, height=7) ## margin changed from 5 to 7 with 00-, 01- for campy
	rfh<-heatmap.2(rfResult, scale="none",trace="none", margin=c(6,8), ColSideColors= cbbPalette[factor(colnames(rfResult), levels=c("00-0949","01-1512","00-1597","00-6200"))], 
		col=colorRampPalette(c("white","lightpink","pink","red"),bias=0.1)(30), keysize=0.8,key.title="% Ident.", cexRow=1.2, cexCol=1.2, main="ResFinder result"
		, cellnote=round(rfResult,1), notecex=1.2, notecol="black")	
dev.off()
