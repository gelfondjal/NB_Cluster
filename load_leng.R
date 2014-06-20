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
##############################################################################


set.seed(2013)
rm(list=ls())

library(DESeq)
library(reshape2)
library(geepack)
library(cummeRbund)
library(ClassComparison)

source.file <- gsub("\\.R","/","load_leng.R")

# define the base directory	

basedir <- "/Users/jonathangelfond/Documents/Projects/Buffenstein/NMR_Genomics/"

# Creates necessary directories

analysisdir <- paste(basedir,"Programs/",sep="") # where the programs are
datadir <- paste(basedir,"Data/Hibbs/countFiles/",sep="")  # where the data are
resultsdir <- paste(basedir,"Documents/Results/",source.file,sep="") # Standard output
tex.dir <- paste(resultsdir,"texdir/",sep="") # Publication quality output

source("/Users/jonathangelfond/Documents/Projects/Joe_Lab/NMR_Methods/Programs/support_functions.R")

#loads libraries and adds common functions to workspace
source(paste(analysisdir,"support_functions.R",sep="")) 


dir.create(resultsdir)
dir.create(tex.dir)

tex.list <- list()

# Assign inpute filenames




install.packages("gplots")
install.packages("blockmodeling")
install.packages("/Users/jonathangelfond/Documents/Projects/Buffenstein/NMR_Genomics/Documents/Background/EBSeq_1.1.6.tar.gz", repos=NULL, type="source")
library(gplots)
library(blockmodeling)
library(EBSeq)



#demo(EBSeq)


data(GeneMat)

 data(IsoList)
 str(IsoList)
  IsoMat=IsoList$IsoMat
   str(IsoMat)
 IsoSizes=MedianNorm(IsoMat)
 IsoNames=IsoList$IsoNames
 
  IsosGeneNames=IsoList$IsosGeneNames
 
  NgList=GetNg(IsoNames, IsosGeneNames)
  
   IsoNgTrun=NgList$IsoformNgTrun
   
    IsoNgTrun[c(1:3,201:203,601:603)]
   
   
   
   
table(table(IsosGeneNames))
 
IsoEBOut=EBTest(Data=IsoMat, NgVector=IsoNgTrun, Conditions=as.factor(rep(c("C1","C2"),each=5)),sizeFactors=IsoSizes, maxround=5)
 
order.gn <- order(IsosGeneNames)
 
iso.sort <- IsoMat[order.gn,]
 
gn.sort <- IsosGeneNames[order.gn]
 
gene.name.table <-  table(gn.sort)
 
table(gene.name.table) 
 
pdf(pastes(resultsdir,"leng_sims.pdf")) 
 
for(isonumber in 2:3){ 
 
gene.2 <- names(gene.name.table[which(gene.name.table==isonumber)])
 
iso.sort.2 <- iso.sort[ gn.sort%in% gene.2,] 

table(gn.sort[gn.sort %in% gene.2])

 
gene.sort.2 <- gn.sort[gn.sort %in% gene.2]
 
# Compute correlations

cors <- matrix(NA,length(unique(gene.sort.2)),1)
rownames(cors) <- unique(gene.sort.2) 

gene.iter <- unique(gene.sort.2)[1]

for(gene.iter in unique(gene.sort.2)){
	
		cors[gene.iter,] <- cor(t(iso.sort.2[gene.sort.2==gene.iter,]))[1,2]
	
	} 
 
 
hist(cors,main=paste("Correlations of Isoforms within gene, n_i=",isonumber) ) 
 
pheno.factor <-  as.factor(rep(c("1","2"),each=5))

iso.sort.2.generows <- iso.sort.2
rownames(iso.sort.2.generows) <- gene.sort.2
 
fit.gee.out <- fit.gee.3.array(iso.sort.2,gene.sort.2,pheno.factor)
 
head(fit.gee.out)


plot(cors,-log10(fit.gee.out[,"pvalue"]),main=paste("Significance vs. Correlation within Isoforms N_i =",isonumber))

log.significance <- -log10(fit.gee.out[,"pvalue"]+1e-15)

lm.fit <- lm(log.significance~poly(cors[1:length(cors)],4))

points(cors,predict(lm.fit),col=2,pch=16)

}

dev.off()

library(cummeRbund)

cuff <- readCufflinks(dir=system.file("extdata", package="cummeRbund"))


count(genes(cuff))

fpkm.matrix <- fpkmMatrix(isoforms(cuff))

gene.fpkm.matrix <- fpkmMatrix(genes(cuff))

genes <- rownames(gene.fpkm.matrix)

#annot.iso <- annotation(isoforms(cuff))

myGeneId<-genes[2]
myGene<-getGene(cuff,myGeneId)

head(fpkm(isoforms(myGene)))


cor.fpkm <- matrix(NA,length(genes),2,dimnames=list(genes,c("cor","n.isoforms")))

gene <- genes[1]

for(gene in genes){
	myGene<-getGene(cuff,gene)
	iso.data <- fpkm(isoforms(myGene))
	
	little.df<- dcast(subset(iso.data,select=c("isoform_id","sample_name","fpkm")),sample_name~isoform_id,value.var="fpkm")
	
	little.mat <- as.matrix(little.df[,-1])
	
	cor.fpkm[gene,"n.isoforms"] <- dim(little.mat)[2]
	
	if(dim(little.mat)[2]>1){
	
	cor.mat <- cor(little.mat)
	
	diag(cor.mat) <- NA
	
	mean.cor.fpkm <- mean(cor.mat,na.rm=TRUE)
	
	cor.fpkm[gene,"cor"] <- mean.cor.fpkm 
	
	}
	
	
	}

par(mfrow=c(1,1))	
hist(cor.fpkm[,"cor"],main="Cufflinks: Isoform correlations within gene",breaks=50)

plot(cor.fpkm[,2],cor.fpkm[,1],main="Number of Isoforms vs. Correlation within Isoforms",xlab="# of Isoforms",ylab="Correlation among isoforms")

df.fpkm <- na.exclude(data.frame(cor.fpkm))

lm.fit <- lm(cor~poly(n.isoforms,1),data=df.fpkm)

pred <- predict(lm.fit)

points(df.fpkm$n.isoforms,pred,col=2,pch=16)
lines(df.fpkm$n.isoforms,pred,col=2,pch=16)



rhos <- seq(-0.9,0.9,length=6)

rho <- rhos[2]

results.cat <- NULL

for(rho in rhos){

nb.size <- 10000

mu <- 10

p <- 2

sigma.mat <- matrix(rho,p,p) + diag(p)*(1-rho)


n.arrays <- 5
delta <- 10

n.clusters <- 500
p.delta <- 0.2


array.data.out <- simulate.2n.arrays(n.arrays,delta,p.delta,p,sigmat.mat,mu,nb.size)

tissue.counts <- array.data.out$array
de.info <- array.data.out$shifter

de.info <- cbind(de.info,de.info[,1]/de.info[,2])

pheno.factor <- as.factor(rep(c(1,2),each=n.arrays))

fit.gee.out <- fit.gee.2.array(tissue.counts,n.clusters,pheno.factor)
fit.gee.perm.out <- fit.gee.2.array.perm(tissue.counts,n.clusters,pheno.factor)

bum.perm <- Bum(fit.gee.perm.out[,2])
alpha.hat <- bum.perm@ahat
l.hat <- bum.perm@lhat
pi0.hat <- bum.perm@pihat

abline(h=pi0.hat,lty=2,col=2)

x <- seq(0,1,length=1000)

y <- l.hat + (1-l.hat)*alpha.hat *x^(alpha.hat-1)

empirical.p.function <- function(x,l.hat,a.hat){
	
	return(l.hat*x+(1-l.hat)*x^a.hat)
	
	}


GEEseq.cor.results <- aggregate(list(DE=empirical.p.function(fit.gee.out[,2],l.hat,alpha.hat)<0.05),list(FC=de.info[,4][seq(1,n.clusters*p,by=p)]),mean)

results.cat <- rbind(results.cat,data.frame(rho=rho,GEEseq.cor.results))

} #loop over rho

library(plotrix)

par(mfrow=c(1,2))
with(subset(results.cat,FC!=1),plot(rho,DE,xlab="Correlation",ylab="Power for DE",main="Power for DE vs. Correlation",pch=19,col=2))
with(subset(results.cat,FC==1),plot(rho,DE,xlab="Correlation",ylab="Type I Error DE",main="Type I Error  vs. Correlation",ylim=c(0,0.1),pch=19,col=2))
abline(h=0.05)

print(paste("EOF:",source.file))






