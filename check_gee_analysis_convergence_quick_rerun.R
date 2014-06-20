##############################################################################
# University of Texas Health Science Center at San Antonio
# Department of Epidemiology and Biostatistics                        
##############################################################################
# Filename:   check_nmr.R                            
# Author:    Jon Gelfond                                             
# Project Name:  
# Input:       
# Output:   Source Filename directory
#
# Modification History:
# v 0.1 Creation
# v 0.2 add weighted glmm
# v 0.3 add permuations
##############################################################################


set.seed(2013)
rm(list=ls())
require(ggplot2)
require(DESeq)
library(reshape2)
library(lme4)
library(ClassComparison)
library(glmmADMB)

source.file <- gsub("\\.R","/","check_gee_analysis_convergence_quick_rerun.R")
source0 <- source.file
# define the base directory	

basedir <- "/Users/jonathangelfond/Documents/Projects/Buffenstein/NMR_Genomics/"

# Creates necessary directories

analysisdir <- paste(basedir,"Programs/",sep="") # where the programs are
datadir <- paste(basedir,"Data/Hibbs/countFiles/",sep="")  # where the data are
resultsdir <- paste(basedir,"Documents/Results/",source.file,sep="") # Standard output
tex.dir <- paste(resultsdir,"texdir/",sep="") # Publication quality output


mus2nmr.file <- "/Users/jonathangelfond/Documents/Projects/Buffenstein/NMR_Genomics/Data/Hibbs/Transcript_Mappings/mmu2nmr.hits.txt"
nmr2mus.file <- "/Users/jonathangelfond/Documents/Projects/Buffenstein/NMR_Genomics/Data/Hibbs/Transcript_Mappings/nmr2mmu.hits.txt"

#loads libraries and adds common functions to workspace
source(paste(analysisdir,"support_functions.R",sep="")) 


dir.create(resultsdir)
dir.create(tex.dir)

tex.list <- list()
workspace.file <- "/Users/jonathangelfond/Documents/Projects/Buffenstein/NMR_Genomics/Documents/Results/check_gee_analysis_1/all.R.Obj.Rdata"

#save(list=ls(),file=workspace.file)
load(workspace.file)


source.file <- source0
basedir <- "/Users/jonathangelfond/Documents/Projects/Buffenstein/NMR_Genomics/"

# Creates necessary directories

analysisdir <- paste(basedir,"Programs/",sep="") # where the programs are
datadir <- paste(basedir,"Data/Hibbs/countFiles/",sep="")  # where the data are
resultsdir <- paste(basedir,"Documents/Results/",source.file,sep="") # Standard output
tex.dir <- paste(resultsdir,"texdir/",sep="") # Publication quality output




source("/Users/jonathangelfond/Documents/Projects/Joe_Lab/NMR_Methods/Programs/support_functions.R")
source("/Users/jonathangelfond/Documents/Projects/sudheer/Analysis/table_functions.R")

library(geepack)
library(gee)
library(MGLM)




#install.packages("gplots")
#install.packages("blockmodeling")
#install.packages("/Users/jonathangelfond/Documents/Projects/Buffenstein/NMR_Genomics/Documents/Background/EBSeq_1.1.6.tar.gz", repos=NULL, type="source")
library(gplots)
library(blockmodeling)
library(EBSeq)




tissue.counts.0 <- expression.mat

pheno.factor <- as.factor(grep.to.label(colnames(tissue.counts.0),"Liver",1,2))

matches.2.1tomany <- subset(matches.2,(Count.Freq>1)&(Count.Freq<5))

tissue.counts <- tissue.counts.0[matches.2.1tomany$V2,]

match.table <- table(rownames(tissue.counts))

match.table <-match.table[match.table<30]

temp <- hist(match.table,breaks=seq(0.5,max(match.table,na.rm=TRUE)+0.5,by=1),main="NMR2MUS 1 to Many Matches")

matplot(rbind(temp$mids,temp$mids),rbind(0,temp$counts),lty=1,col=1,type="l",xlab="# of MUS matches",ylab="Count",main="NMR2MUS 1 to Many Matches")
points(temp$mids,temp$counts,pch=16)
text(temp$mids,temp$counts-100,temp$counts,col=2)






tissue.counts.nodups <- tissue.counts[unique(matches.2.1tomany$V2),]

#fit.gee.out <- fit.gee.4.array(tissue.counts,matches.2.1tomany$V1,pheno.factor,crash.list=c(4300:4400,5300:5400,8300:8400))

tissue.counts.norm <- normalize.quantiles(tissue.counts)

dimnames(tissue.counts.norm) <- dimnames(tissue.counts)

n.transcripts <- 1000


#heatmap.2(log2(tissue.counts.norm[1:n.transcripts,]+0.1))

families <- list(GEE="GEE",quasipossion=quasipoisson,nbinom1="nbinom1")


tissue.counts.norm <- tissue.counts.norm + 0.1

analysis.out <- lapply(families,function(x){fit.gee.5.array(tissue.counts.norm,matches.2.1tomany,pheno.factor,crash.list=NULL,family=x)})


analysis.out.perm <- lapply(families,function(x){fit.gee.5.array(tissue.counts.norm,matches.2.1tomany,pheno.factor,crash.list=1:n.transcripts,family=x,perm=4)})


single.analysis.out <- fit.singleton.5.array(tissue.counts.norm,matches.2.1tomany,pheno.factor,crash.list=NULL,weight=FALSE,perm=0,nb=TRUE)
single.analysis.out.perm <- fit.singleton.5.array(tissue.counts.norm,matches.2.1tomany,pheno.factor,crash.list=NULL,weight=FALSE,perm=1,nb=TRUE)

corrected.singleton.ps <- ecdf(single.analysis.out.perm[,"pvalue"])(single.analysis.out[,"pvalue"])

pdf(paste0(resultsdir,"comparing_count_models_prod.pdf"),height=8,width=11)


par(mfrow=c(2,3))

vars <- c("FC","pvalue")

for(var.iter in vars){

done <- names(analysis.out)

for(fam.iter in names(analysis.out)){for(fam.iter.2 in setdiff(done,fam.iter)){
		
plot(log10(analysis.out[[fam.iter]][[1]][,var.iter]),		log10(analysis.out[[fam.iter.2]][[1]][,var.iter]),main=paste(var.iter,fam.iter,fam.iter.2),
					xlab=paste("log10",var.iter,fam.iter),ylab=paste("log10",var.iter,fam.iter.2))

#done <- setdiff(done,fam.iter.2)

abline(0,1)
		
}} #loop over fam.iter
} #loop over var.iter		



perm.combo <- list()
for(namer in names(analysis.out)){perm.combo[[namer]] <- list(orig=analysis.out[[namer]][[1]],perm=analysis.out.perm[[namer]][[1]])}

x <- perm.combo[[1]]

par(mfrow=c(2,3))
for(namer in names(analysis.out)){
	nas <- table(is.na(perm.combo[[namer]]$orig[,"pvalue"]))
	na.out <- try(paste(nas["TRUE"],"/",sum(nas)))
	hist(perm.combo[[namer]]$orig[,"pvalue"],main=paste("Original p-value",namer,"\n # Convergence Failures =",na.out))}
for(namer in names(analysis.out)){
	hist(perm.combo[[namer]]$perm[,"pvalue"],main=paste("Perumuted p-value",namer))
	}


par(mfrow=c(2,3))
corrected.ps <- lapply(perm.combo,function(x){return(p.bum.correct.right.truncate(x$orig[,"pvalue"],x$perm[,"pvalue"]))})
for(namer in names(analysis.out)){
	
	hist(corrected.ps[[namer]],main=paste("Corrected p-value",namer))
	
	
	}

par(mfrow=c(1,1))




plot(sort(-log10(perm.combo[[namer]]$orig[,"pvalue"])),main=paste("Original p-value"),type="n",ylab="-log10")
step <- 1
for(namer in names(analysis.out)){
	
	lines(sort(-log10(perm.combo[[namer]]$orig[,"pvalue"]),decreasing=TRUE),main=paste("Original p-value",namer),col=step)
	step <- step + 1
	
	}
legend("topleft",legend=names(analysis.out),lty=1,col=1:length(corrected.ps))




plot(sort(-log10(corrected.ps[[namer]])),ylim=c(0,50),main=paste("Corrected p-value"),type="n",ylab="-log10")
step <- 1
for(namer in names(analysis.out)){
	
	lines(sort(-log10(corrected.ps[[namer]]),decreasing=TRUE),main=paste("Corrected p-value",namer),col=step)
	step <- step + 1
	
	}
legend("topleft",legend=names(analysis.out),lty=1,col=1:length(corrected.ps))




dev.off()



lapply(corrected.ps,function(x){sum(x<0.05,na.rm=TRUE)})

lapply(corrected.ps,function(x){sum(!is.na(x))})
lapply(corrected.ps,function(x){mean(x<0.05,na.rm=TRUE)})
temp <- lapply(corrected.ps,qvalue)


save(list=ls(),file=file.path(resultsdir,"results_workspace.RData"))



print(paste("EOF:",source.file))






