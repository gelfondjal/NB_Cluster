##############################################################################
# University of Texas Health Science Center at San Antonio
# Department of Epidemiology and Biostatistics                        
##############################################################################
# Filename:   simulate_neg_binom_v0.1.R                            
# Author:    Jon Gelfond                                             
# Project Name:  
# Input:       
# Output:   Source Filename directory
#
# Modification History:
# v 0.1 Creation 
# v 0.2 add FDR adjustment
# v 0.3 look at real data bloomy
##############################################################################



library(geepack)
library(ClassComparison)
source("http://bioinformatics.mdanderson.org/OOMPA/oompaLite.R")
oompaLite()
library(BMA)

set.seed(2013)
rm(list=ls())


source.file <- gsub("\\.R","/","simulate_neg_binom_v0.3.R")

# define the base directory	

basedir <- "/Users/jonathangelfond/Documents/Projects/Joe_Lab/NMR_Methods/"

# Creates necessary directories

analysisdir <- paste(basedir,"Programs/",sep="") # where the programs are
datadir <- paste(basedir,"Data/",sep="")  # where the data are
resultsdir <- paste(basedir,"Documents/Results/",source.file,sep="") # Standard output
tex.dir <- paste(resultsdir,"texdir/",sep="") # Publication quality output


#loads libraries and adds common functions to workspace
source(paste(analysisdir,"support_functions.R",sep="")) 

bottomly.dir <- "/Users/jonathangelfond/Documents/Projects/Joe_Lab/NMR_Methods/Data/Bottomly_Data/"


dir.create(resultsdir)
dir.create(tex.dir)

tex.list <- list()


# Assign input filenames


n <- 100000
sizer <- 200
mu <- 10

temp <- rnbinom(n,size=sizer,mu=mu)

mean(temp)
var(temp)
mu + mu^2/sizer


p <- 2

sigma.mat <- matrix(1,p,p) + diag(p)


temp <- rmvnorm(n,mean=rep(0,p),sigma=sigma.mat)



# Construct Correlated negative binomial



r.cor.nb <- function(n,p,sigma.mat,mu,nb.size){
	
	# n = # of samples
	# p = dimension of NB vector
	# sigma.mat correlation
	# mu = mean of NB distribution
	# nb.size = size parameter of NB distribution
	
	
	scaler <- exp(rmvnorm(n,mean=rep(0,p),sigma=sigma.mat))
	
	x.nb <- round(matrix(rnbinom(n*p,size=nb.size,mu=mu),n,p)*scaler,0)


	return(x.nb)
	
	
	
	} # END: r.cor.nb
	

n <- 1e2
p <- 2
rho <- 0.6
sigma.mat <- matrix(rho,p,p) + diag(p)*(1-rho)
nb.size <- 10000
mu <- 10

cor.nb.out <- r.cor.nb(n,p,sigma.mat,mu,nb.size)

tall.out <- as.vector(cor.nb.out)
cluster.id <- as.vector(matrix(1:n,n,p))

n.arrays <- 5
delta <- 10

n.clusters <- 1000
p.delta <- 0.2


array.data.out <- simulate.2n.arrays(n.arrays,delta,p.delta,p,sigmat.mat,mu,nb.size)

tissue.counts <- array.data.out$array
de.info <- array.data.out$shifter


de.info <- cbind(de.info,de.info[,1]/de.info[,2])

pheno.factor <- as.factor(rep(c(1,2),each=n.arrays))

colnames(tissue.counts) <- paste("c",pheno.factor,rep(1:n.arrays,2),sep=".")
heatmap.2(log(tissue.counts+0.5),col=redblue(256),trace="none",title="Simulated Data")

hist(log10(tissue.counts))

vars <- apply(tissue.counts,1,var)
means <- apply(tissue.counts,1,mean)

plot(log(means),log(vars))

library(gplots)
library(blockmodeling)
library(EBSeq)

Sizes=MedianNorm(tissue.counts)

itermax <- 5


library(biomaRt)
library(devtools)
temp <- digest("/Users/jonathangelfond/Documents/Projects/Ireland/Documents/Results/R2SAS/LPG_STYLE.sas",file=TRUE,algo="sha1",serialize=FALSE)

nchar(temp)

bottomly.counts <- read.table(pastes(bottomly.dir,"bottomly_count_table.txt"),header=TRUE)

bottomly.mat <- as.matrix(bottomly.counts[,-1])
colnames(bottomly.mat) <- names(bottomly.counts)[-1]
rownames(bottomly.mat) <- bottomly.counts[,1]

bottomly.pheno <- read.table(pastes(bottomly.dir,"bottomly_phenodata.txt"),header=TRUE)

identical(colnames(bottomly.mat),as.character(bottomly.pheno$sample.id))

mart<- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))

ensembl_genes<- rownames(bottomly.mat)

bottomly.gene.info <- getBM(
  filters= "ensembl_gene_id", 
  attributes= c("ensembl_gene_id", "external_gene_id", "entrezgene", "description"),
  values= ensembl_genes,
  mart= mart)


listdups <- bottomly.gene.info$external_gene_id[duplicated(bottomly.gene.info$external_gene_id)]

dup.ids <- unique(subset(bottomly.gene.info,external_gene_id %in% listdups)$ensembl_gene_id)

dup.counts <- bottomly.mat[dup.ids,]

bottomly.dup <- subset(bottomly.gene.info,ensembl_gene_id%in%dup.ids)

bottomly.dup <- subset(bottomly.dup ,!duplicated(bottomly.dup$ensembl_gene_id))

listdups <- bottomly.dup$external_gene_id[duplicated(bottomly.dup$external_gene_id)]


dup.ids <- unique(subset(bottomly.dup,external_gene_id %in% listdups)$ensembl_gene_id)


bottomly.dup.mat <- bottomly.mat[dup.ids,]


bottomly.dup.mat.nozero <- bottomly.dup.mat[rowSums(bottomly.dup.mat)>0,]

no.dup.genes.nozero <- subset(bottomly.gene.info, ensembl_gene_id %in% rownames(bottomly.dup.mat.nozero))

nz.dups <- subset(no.dup.genes.nozero,external_gene_id %in% no.dup.genes.nozero$external_gene_id[duplicated(no.dup.genes.nozero$external_gene_id)])

bottomly.dup.mat.nozero <- bottomly.dup.mat[nz.dups$ensembl_gene_id,]

dim(nz.dups)
dim(bottomly.dup.mat.nozero)

cbind(nz.dups$ensembl_gene_id,rownames(bottomly.dup.mat.nozero))

Sizes=MedianNorm(bottomly.dup.mat.nozero)


EBOut <- EBTest(Data=bottomly.dup.mat.nozero,Conditions=bottomly.pheno$strain,sizeFactors=Sizes, maxround=itermax)


PP <- GetPPMat(EBOut)
post.fc <- PostFC(EBOut)


# Fit model with GEE

cluster.names <- nz.dups$external_gene_id
clust.iter <- cluster.names[1]

n.clusters <- length(unique(cluster.names))

results.fc <- matrix(NA,n.clusters,2)
colnames(results.fc) <- c("FC","pvalue")
rownames(results.fc) <- unique(cluster.names)

pheno.factor <- bottomly.pheno$strain

library(reshape2)

for(clust.iter in rownames(results.fc)){

gene.data <- Sortby(data.frame(melt(t(bottomly.dup.mat.nozero[cluster.names==clust.iter,])),pheno.factor),"Var1")

gee.out <- geeglm(value~pheno.factor+Var2,data=gene.data,id=Var1,family=poisson)

gee.sum <- summary(gee.out)$coef

results.fc[clust.iter,] <- c(exp(-gee.sum["pheno.factorDBA/2J","Estimate"]),gee.sum["pheno.factorDBA/2J","Pr(>|W|)"])

}

PP2 <- data.frame(PP,ensembl_gene_id=rownames(PP))

post.fc2 <- data.frame(post.fc$PostFC,ensembl_gene_id=names(post.fc$PostFC))


PP.join <- merge(PP2,subset(nz.dups,select=c("external_gene_id","ensembl_gene_id")),by="ensembl_gene_id")
post.fc2.join <- merge(post.fc2,subset(nz.dups,select=c("external_gene_id","ensembl_gene_id")),by="ensembl_gene_id")

post.fc2.results <- merge(post.fc2.join ,data.frame(external_gene_id=rownames(results.fc),results.fc),by="external_gene_id")

melted.fc <- melt(subset(post.fc2.results,select=c("external_gene_id","post.fc.PostFC")),id.vars="external_gene_id")

fc.lm <- lm(value~external_gene_id,data=melted.fc)
Anova(fc.lm)

summary()

wider.fc <- dcast(melted.fc,external_gene_id ~variable+.)

PP.results <- merge(PP.join ,data.frame(external_gene_id=rownames(results.fc),results.fc),by="external_gene_id")

PP.results$pvalue <- PP.results$pvalue + 0.00000001


with(PP.results,plot(PPDE,-log10(pvalue),main="GEE combined vs. EBSeq individual calls"))
abline(h=-log10(0.05),v=0.5)
with(PP.results,text(PPDE,-log10(pvalue),external_gene_id,col="red"))

with(post.fc2.results,plot(post.fc.PostFC,FC,main="GEE combined vs. EBSeq individual calls"))
abline(0,1)
with(post.fc2.results,text(post.fc.PostFC,FC,external_gene_id,col="red"))
abline(h=1,v=1,lty=2)

with(post.fc2.results,cor.test(post.fc.PostFC,FC))


hist(results.fc[,"pvalue"],main="GEE p-values")
legend("topright",legend=paste(100*round(mean(results.fc[,"pvalue"]<0.05),2),"% P<0.05"),bg="white")



fit.gee.out <- fit.gee.2.array(bottomly.dup.mat.nozero,n.clusters,pheno.factor)
fit.gee.perm.out <- fit.gee.2.array.perm(bottomly.dup.mat.nozero,n.clusters,pheno.factor)





par(mfrow=c(1,2))

hist(fit.gee.out[,"pvalue"],breaks=50)

hist(fit.gee.perm.out[,"pvalue"],breaks=50)

par(mfrow=c(1,2))

plot(log(fit.gee.out[,1]),log(fit.gee.perm.out[,1]))
plot((fit.gee.out[,1]),(fit.gee.perm.out[,1]))

cor.test(log(fit.gee.out[,1]),log(fit.gee.perm.out[,1]))
cor.test((fit.gee.out[,1]),(fit.gee.perm.out[,1]))


cor.test(-log(fit.gee.out[,2]),-log(fit.gee.perm.out[,2]))
plot(-log(fit.gee.out[,2]),-log(fit.gee.perm.out[,2]))

bum.perm <- Bum(fit.gee.perm.out[,2])
par(mfrow=c(1,1))

hist(bum.perm,main="Permuted P-values as Beta-Uniform Mixture")

alpha.hat <- bum.perm@ahat
l.hat <- bum.perm@lhat
pi0.hat <- bum.perm@pihat

abline(h=pi0.hat,lty=2,col=2)

x <- seq(0,1,length=1000)

y <- l.hat + (1-l.hat)*alpha.hat *x^(alpha.hat-1)

empirical.p.function <- function(x,l.hat,a.hat){
	
	return(l.hat*x+(1-l.hat)*x^a.hat)
	
	}

hist(empirical.p.function(fit.gee.perm.out[,2],l.hat,alpha.hat),main="Corrected Permuted P-Values")

bum.perm.cor <- Bum(empirical.p.function(fit.gee.perm.out[,2],l.hat,alpha.hat))

hist(bum.perm.cor,main="Corrected Permuted P-values as Beta-Uniform Mixture")

par(mfrow=c(1,2))

bum.original <- Bum(fit.gee.out[,2])
hist(bum.original,main="P-values as Beta-Uniform Mixture")


bum.original.cor <- Bum(empirical.p.function(fit.gee.out[,2],l.hat,alpha.hat))
hist(bum.original.cor,main="Corrected P-values as Beta-Uniform Mixture")


lines(x,y,col=2,lty=2)

bum.original@pihat
bum.original.cor@pihat

GEEseq.results <- aggregate(list(DE=fit.gee.out[,2]<0.05),list(FC=de.info[,4][seq(1,n.clusters*p,by=p)]),mean)

GEEseq.cor.results <- aggregate(list(DE=empirical.p.function(fit.gee.out[,2],l.hat,alpha.hat)<0.05),list(FC=de.info[,4][seq(1,n.clusters*p,by=p)]),mean)

tex.list <- lcat(tex.list,latex.out.table.2(GEEseq.results,
	file=NULL,
	caption = paste("GEE results alpha=0.05"),
	label = "",
	rownamer = NULL,
	col.label=NULL))
	
tex.list <- lcat(tex.list,latex.out.table.2(GEEseq.cor.results,
	file=NULL,
	caption = paste("GEE Corrected results alpha=0.05"),
	label = "",
	rownamer = NULL,
	col.label=NULL))
	

EBseq.results <- aggregate(list(DE=1*(PP[,"PPDE"]>=0.95)),list(FC=de.info[,4]),mean)


tex.list <- lcat(tex.list,latex.out.table.2(EBseq.results,
	file=NULL,
	caption = paste("EBSeq Results ppde=95"),
	label = "",
	rownamer = NULL,
	col.label=NULL))


table.files <- write.includer.vector(tex.list,path=tex.dir,include.file.base="includemetables")
make.latex.doc.vector(pastes(tex.dir,"summary.tex"),includer=table.files,title="Simulated Correlated Negative binomial Analysis",author="Jon Gelfond",date=as.character(Sys.time()),path=NULL)


print(paste("EOF:",source.file))






