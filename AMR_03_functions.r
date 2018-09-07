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


### functions
require(gplots)
dist.pear <- function(x) as.dist(1-cor(t(x)))
	

jFiltContRev<-function(mfile){
	#mfile[mfile$Reverse!="+"&mfile$Potential.contaminant!="+",] ## edited on 20180403
	mfile[(is.na( mfile$Reverse)|mfile$Reverse!="+")& (is.na(mfile$Potential.contaminant)|mfile$Potential.contaminant!="+"),] ## edited on 20180403
	
}
jlog2<-function(y){
	y$dat=log2(y$dat)
	y$dat[y$dat==-Inf]=NA
	y
}
jFilterValid<-function(y, prop=0.7){
	dat<-y$dat
	tokeep<-rowSums(!is.na(dat))>=(ncol(dat)*prop)
	y$dat<-dat[tokeep,]
	y$rannot<-y$rannot[tokeep,]
	y
	
}
########
### Median norm of a matrix, columns/samples
########
jMedianNormBatch <-function(mat0, batch=NULL){
	ret<-mat0
	if(is.null(batch)){
		omed<-median(mat0, na.rm = TRUE)
		for (j in 1:ncol(mat0)){
			ret[, j] <- mat0[,j] - median(mat0[, j], na.rm = TRUE) +omed
		}
	}else{
		print("medNormColBatch")
		totbatch<-sapply(unique(batch), function(y){
			median(mat0[,batch==y])
		})
		print(totbatch)
		omed<-median(totbatch)
		print(omed)
		for(i in 1:length(unique(batch))){
			ret[,batch==unique(batch)[i]]=mat0[,batch==unique(batch)[i]]-totbatch[i]+omed
		}
	}
	ret
}

### when design in parallel, jMedianNormBatchRow , normalize the medians of each protein between Runs/Batches
jMedianNormBatchRow <-function(mat0, batch=NULL){
	ret<-mat0
	print("medNormRowBatch")
	totbatch<-sapply(unique(batch), function(y){
		apply(mat0[,batch==y], 1, median)
	}) #protein x run matrix
	#print(totbatch)
	omed<-apply(totbatch,1,median) #1*protein
	#print(omed)
	for(j in 1:length(omed)){
		for(i in 1:length(unique(batch))){
			ret[j,batch==unique(batch)[i]]=mat0[j,batch==unique(batch)[i]]-totbatch[j,i]+omed[j]
		}
	}

	ret
}


jNormalizeLoess<-function(dat4stat, dat4statFull, save2pdf=T){
		require(limma)
		#jMAplot(dat4stat, pheno=pheno)
		dat4stat<-normalizeCyclicLoess(dat4stat) ## not taking into account of groups, which makes sense
		   boxplot(dat4stat, main=paste("Loess normed"), col= cbbPalette[pheno$groupid], ylab="Log2 normed Intensity")
                        abline(h= median(t(dat4stat)),lty=2)


        #### replacing data in dat4statFull with the normalized dat4stat data
        for(i in colnames(dat4stat)){
                dat4statFull[,(colnames(dat4statFull)==i|colnames(dat4statFull)==paste0("X",i))]=dat4stat[,i]
        }
        ### remove stdev column: Std Dev
        if(length(grep("Std Dev",colnames(dat4statFull)))>0){
                dat4statFull<-dat4statFull[,-grep("Std Dev",colnames(dat4statFull))]
        }
        list(dat4stat=dat4stat, dat4statFull=dat4statFull)
}

toPeptideTable<-function(y,dat0, expNames){
	dat=dat0[dat0$Sequence==y,]
	iid=grep("Reporter.intensity.corrected",colnames(dat))
	ret=NULL
	for(expn in expNames){
		if(sum(dat$Experiment ==expn)>0){
			dat1=dat[dat$Experiment==expn,iid];
			colnames(dat1)=paste(gsub("Reporter.intensity.corrected.","",colnames(dat1)), ".",expn,sep="")
			ret=c(ret,log2(colSums(dat1)+1))
			ret[ret==0]=NA
		}else{
			tmp0<-rep(NA,length(iid)); names(tmp0)=paste(0:3,expn,sep=".") #0:3 for iTRAQ
			ret=c(ret,tmp0)
		}
	}
	ret
}
