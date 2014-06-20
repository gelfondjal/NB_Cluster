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

source.file <- gsub("\\.R","/","check_gee_analysis_convergence.R")

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
source("/Users/jonathangelfond/Documents/Projects/Joe_Lab/NMR_Methods/Programs/support_functions.R")
source("/Users/jonathangelfond/Documents/Projects/sudheer/Analysis/table_functions.R")

library(geepack)
library(gee)




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

n.transcripts <- 3000


#heatmap.2(log2(tissue.counts.norm[1:n.transcripts,]+0.1))

families <- list(poisson=poisson,quasipoisson=quasipoisson,nbinom1="nbinom1")

analysis.out <- lapply(families,function(x){fit.gee.5.array(tissue.counts.norm,matches.2.1tomany,pheno.factor,crash.list=1:n.transcripts,family=x)})


analysis.out.perm <- lapply(families,function(x){fit.gee.5.array(tissue.counts.norm,matches.2.1tomany,pheno.factor,crash.list=1:n.transcripts,family=x,perm=4)})

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
	hist(perm.combo[[namer]]$orig[,"pvalue"],main=paste("Original p-value",namer))}
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
	
	lines(sort(-log10(perm.combo[[namer]]$orig[,"pvalue"])),main=paste("Original p-value",namer),col=step)
	step <- step + 1
	
	}
legend("topleft",legend=names(analysis.out),lty=1,col=1:length(corrected.ps))




plot(sort(-log10(corrected.ps[[namer]])),main=paste("Corrected p-value"),type="n",ylab="-log10")
step <- 1
for(namer in names(analysis.out)){
	
	lines(sort(-log10(corrected.ps[[namer]])),main=paste("Corrected p-value",namer),col=step)
	step <- step + 1
	
	}
legend("topleft",legend=names(analysis.out),lty=1,col=1:length(corrected.ps))




dev.off()




gee.corrected.p <- p.bum.correct.right.truncate(fit.gee.out$pvalue,fit.gee.out.perm$pvalue)
glm.corrected.p <- p.bum.correct.right.truncate(fit.glm.out$pvalue,na.exclude(fit.glm.out.perm$pvalue))


fit.poisson.out <- fit.gee.5.array(tissue.counts.norm,matches.2.1tomany,pheno.factor,crash.list=1:n.transcripts)


fit.glmmPQL.out <- fit.gee.5.array(tissue.counts.norm,matches.2.1tomany,pheno.factor,crash.list=1:n.transcripts,family=quasipoisson)

fit.glmmNB.out <- fit.gee.5.array(tissue.counts.norm,matches.2.1tomany,pheno.factor,crash.list=1:n.transcripts,family="nbinom1")







fit.gee.out.qp <- fit.gee.5.array(tissue.counts.norm,matches.2.1tomany,pheno.factor,crash.list=1:n.transcripts,family=quasipoisson)


fit.gee.out.perm <- fit.gee.5.array(tissue.counts.norm,matches.2.1tomany,pheno.factor,crash.list=1:n.transcripts,perm=4)

with(data.frame(fit.gee.out[[1]]),plot(log(FC),-log2(pvalue),pch=16))


fit.glm.out <-fit.singleton.5.array(tissue.counts.norm,matches.2.1tomany,pheno.factor,crash.list=1:n.transcripts,perm=0,nb=TRUE)

fit.glm.out.perm <-fit.singleton.5.array(tissue.counts.norm,matches.2.1tomany,pheno.factor,crash.list=1:n.transcripts,perm=2,nb=TRUE)


pdf(paste0(resultsdir,"Histograms_of_pvalues_nb.pdf"),height=8,width=11)

par(mfrow=c(2,2))


hist(na.exclude(fit.gee.out$pvalue),main="Original Skin vs Liver GEE-cluster")
hist(fit.gee.out.perm$pvalue,main="Permuted Skin vs Liver GEE")


hist(na.exclude(fit.glm.out$pvalue),main="Original Skin vs Liver GLM-singleton")
hist(fit.glm.out.perm$pvalue,main="Permuted Skin vs Liver GLM-singleton")


gee.corrected.p <- p.bum.correct.right.truncate(fit.gee.out$pvalue,fit.gee.out.perm$pvalue)
glm.corrected.p <- p.bum.correct.right.truncate(fit.glm.out$pvalue,na.exclude(fit.glm.out.perm$pvalue))
par(mfrow=c(1,2))

hist(gee.corrected.p)
hist(glm.corrected.p)

par(mfrow=c(1,2))

plot.new()

legend.out <- to.char.matrix(matrix(c(mean.x(gee.corrected.p<0.01),mean.x(glm.corrected.p<0.01)),2,1,dimnames=list(c("GEE","GLM"),"P<0.01")))
mtext("Proportion of Significant Calls")
legend("center",legend= legend.out,ncol=ncol(legend.out))

dev.off()


highly_sig_genes <- subset(fit.glm.out,pvalue==0)$V1
pdf(paste0(resultsdir,"most_significant_transcripts.pdf"),h=8,w=11)
plot.singleton.5.array(tissue.counts.norm,matches.2.1tomany,pheno.factor,transcipt.list=highly_sig_genes )
dev.off()



pdf(paste0(resultsdir,"most_significant_perm_transcripts.pdf"),h=8,w=11)
plot.singleton.5.array(tissue.counts.norm,matches.2.1tomany,pheno.factor,transcipt.list=subset(fit.glm.out,pvalue<=1e-16)$V1,perm=1 )
dev.off()

var.transcripts <- apply(tissue.counts.norm,1,var,na.rm=TRUE)

mean.transcripts <- apply(tissue.counts.norm,1,mean,na.rm=TRUE)
par(mfrow=c(1,1))
plot(log(mean.transcripts),log(var.transcripts),main="Count Variance vs. Mean log scale")
abline(0,1,lwd=2,col="red")

pheno.treat <- c(1,1,1,2,2)

fit.gee.out.liver <- fit.gee.5.array(tissue.counts.norm[,1:5],matches.2.1tomany,pheno.treat,crash.list=1:n.transcripts)
fit.gee.out.liver.perm <- fit.gee.5.array(tissue.counts.norm[,1:5],matches.2.1tomany,pheno.treat,crash.list=1:n.transcripts,perm=2)

fit.gee.out.skin <- fit.gee.5.array(tissue.counts.norm[,6:10],matches.2.1tomany,pheno.treat,crash.list=1:n.transcripts)
fit.gee.out.skin.perm <- fit.gee.5.array(tissue.counts.norm[,6:10],matches.2.1tomany,pheno.treat,crash.list=1:n.transcripts,perm=2)


par(mfrow=c(1,2))

hist(na.exclude(fit.gee.out$pvalue),main="Original")
hist(fit.gee.out.perm$pvalue,main="Permuted")


hist(na.exclude(fit.gee.out.liver$pvalue),main="Liver Original")
hist(fit.gee.out.liver.perm$pvalue,main="Liver Permuted")



hist(na.exclude(fit.gee.out.skin$pvalue),main="Skin Original")
hist(fit.gee.out.skin.perm$pvalue,main="Skin Permuted")


fit.glm.out <-fit.singleton.5.array(tissue.counts.norm,matches.2.1tomany,pheno.factor,crash.list=1:n.transcripts,perm=0)

fit.glm.out.perm <-fit.singleton.5.array(tissue.counts.norm,matches.2.1tomany,pheno.factor,crash.list=1:n.transcripts,perm=2)


hist(fit.glm.out$pvalue,main="Singleton Original")

hist(fit.glm.out.perm$pvalue,main="Singleton Perm")






#fit.gee.out.w <- fit.gee.5.array(tissue.counts.norm,matches.2.1tomany,pheno.factor,crash.list=1:1000,weight=TRUE)


plot(log2(fit.gee.out.w$FC),log2(fit.gee.out$FC),main="log2(FC) Weighted vs NotWeighted")
abline(0,1)


plot(-log10(fit.gee.out.w$pvalue),-log10(fit.gee.out$pvalue),main="log10(Pval) Weighted vs NotWeighted")
abline(0,1)


#fit.gee.perm.out <- fit.gee.2.array.perm(tissue.counts,n.clusters,pheno.factor)

fit.merge <- merge(fit.gee.out,cor.calc.out,by="V1")

fit.merge$log10pvalue <- -log10(fit.merge$pvalue+1e-15)

hist(fit.merge$pvalue,breaks=50)

ggplot2::qplot(mean.cor,log10pvalue,data=na.exclude(fit.merge),geom=c("point","smooth")) + labs(title="Log10 p-value DE vs Mean Correlation of Cluster")



ggplot2::qplot(mean.cor,log2(FC),data=na.exclude(fit.merge),geom=c("point","smooth")) + labs(title="Log2 FC DE vs Mean Correlation of Cluster")


subset(matches.2.1tomany,V2=="JH602108.222")

subset(fit.merge,V1=="uc007lip.1")


Sizes=MedianNorm(tissue.counts.nodups)


EBOut <- EBTest(Data=tissue.counts.nodups,Conditions=pheno.factor,sizeFactors=Sizes, maxround=5)

PP <- GetPPMat(EBOut)
eb.fc <- PostFC(EBOut)[[1]]
post.fc <- data.frame(V2=names(eb.fc),EB.FC=eb.fc)

hist(PP[,"PPDE"])

pp.df <- data.frame(V2=rownames(PP),PPDE=PP[,"PPDE"])



pp.df$EB.FC <- ifelse(post.fc$V2==pp.df$V2,post.fc$EB.FC,NA)

pp.df$log.EB.FC <- log2(pp.df$EB.FC)


head(pp.df[order(abs(pp.df$log.EB.FC),decreasing=TRUE),])

pp.df.mg <- merge(pp.df,matches.2.1tomany,by="V2")

agg.matches <- aggregate(subset(pp.df.mg,select=c("PPDE","log.EB.FC")),by=subset(pp.df.mg,select="V1"),mean,na.rm=TRUE)

merge.comp <- merge(agg.matches,fit.gee.out,by="V1")

merge.comp.w <- merge(agg.matches,fit.gee.out.w,by="V1")



par(mfrow=c(1,1))

with(merge.comp,plot(PPDE,-log10(pvalue),main="Cluster mean PPDE (EBSeq) vs. -log10(pvalue) Cluster (GEE)"))

legend("topleft",legend=paste("Cor =",round(with(merge.comp,cor.test(PPDE,-log10(pvalue)))$estimate,2)))

par(mfrow=c(1,1))

with(merge.comp.w,plot(PPDE,-log10(pvalue),main="Cluster mean PPDE (EBSeq) vs. -log10(pvalue) Cluster (GEE)"))

legend("topleft",legend=paste("Cor =",round(with(merge.comp.w,cor.test(PPDE,-log10(pvalue)))$estimate,2)))





#with(merge.comp,cor.test(PPDE,-log10(pvalue+1e-15)))

hist(log2(pp.df$EB.FC),xlim=c(0,10))

par(mfrow=c(1,2))

with(merge.comp,plot(log2(FC),log.EB.FC,main="Cluster mean Fold-Change (EBSeq) vs. Fold-Change  Cluster (GEE)"))

cor.fc <- with(merge.comp,cor.test(log2(FC),log.EB.FC))$estimate
legend("topleft",legend=paste("Cor =",round(cor.fc,2)))
abline(0,1)



with(merge.comp.w,plot(log2(FC),log.EB.FC,main="Weighted Cluster mean Fold-Change (EBSeq) vs. Fold-Change  Cluster (GEE)"))

cor.fc <- with(merge.comp.w,cor.test(log2(FC),log.EB.FC))$estimate
legend("topleft",legend=paste("Cor =",round(cor.fc,2)))
abline(0,1)




print(paste("EOF:",source.file))






