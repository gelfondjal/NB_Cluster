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
##############################################################################



library(geepack)
library(ClassComparison)
source("http://bioinformatics.mdanderson.org/OOMPA/oompaLite.R")
oompaLite()
library(BMA)

set.seed(2013)
rm(list=ls())


source.file <- gsub("\\.R","/","simulate_neg_binom_v0.2.R")

# define the base directory	

basedir <- "/Users/jonathangelfond/Documents/Projects/Joe_Lab/NMR_Methods/"

# Creates necessary directories

analysisdir <- paste(basedir,"Programs/",sep="") # where the programs are
datadir <- paste(basedir,"Data/",sep="")  # where the data are
resultsdir <- paste(basedir,"Documents/Results/",source.file,sep="") # Standard output
tex.dir <- paste(resultsdir,"texdir/",sep="") # Publication quality output


#loads libraries and adds common functions to workspace
source(paste(analysisdir,"support_functions.R",sep="")) 




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


EBOut <- EBTest(Data=tissue.counts,Conditions=pheno.factor,sizeFactors=Sizes, maxround=itermax)


PP <- GetPPMat(EBOut)
post.fc <- PostFC(EBOut)

mean.df <- data.frame(gene=names(EBOut$C1Mean[[1]]),gene2=names(EBOut$C2Mean[[1]]),Drug.mean=EBOut$C1Mean[[1]],Placebo.mean=EBOut$C2Mean[[1]])

if(identical(mean.df$gene,mean.df$gene2)){mean.df <- Drop.variable(mean.df,"gene2")}

post.fc.df <- data.frame(gene=names(post.fc$PostFC),post.fc$PostFC)


tissue.df <- data.frame(gene=rownames(tissue.counts),tissue.counts)


PPsort <- merge(data.frame(gene=rownames(PP),PP),post.fc.df,by="gene")

PPsort <- merge(PPsort,mean.df,by="gene")

PPsort <- merge(PPsort,tissue.df,by="gene")


DEfound=rownames(PP)[which(PP[,"PPDE"]>=.95)]

#write.csv(PPsort,paste(resultsdir,tissue,data.name,".csv",sep="_"))

par(mfrow=c(1,1))
hist(PP[,"PPDE"],xlab="Posterior Probability DE",main="EBSeq")
legend("topright",legend=paste(100*round(prop.table(table(PP[,"PPDE"]>0.5)),2),"%",c("Not DE","DE")),bg="white")

ave.mean <- 0.5*(EBOut$C1Mean[[1]]+EBOut$C2Mean[[1]])


boxplot(post.fc[[1]]~de.info[,4],log="y",main="True vs. Estimated FC")
abline(h=unique(de.info[,4]))

library(ggplot2)
library(reshape2)

qplot(log10(post.fc[[1]]),1*(PP[,"PPDE"]>=.95),geom=c("point","smooth"),formula=y~ns(x,5))

EBseq.results <- aggregate(list(DE=1*(PP[,"PPDE"]>=0.9)),list(FC=de.info[,4]),mean)


tex.list <- lcat(tex.list,latex.out.table.2(EBseq.results,
	file=NULL,
	caption = paste("EBSeq Results"),
	label = "",
	rownamer = NULL,
	col.label=NULL))


# Fit model with GEE

cluster.names <- gsub("\\.[0-9]$","",rownames(tissue.counts))
clust.iter <- cluster.names[1]

results.fc <- matrix(NA,n.clusters,2)
colnames(results.fc) <- c("FC","pvalue")
rownames(results.fc) <- unique(cluster.names)

for(clust.iter in rownames(results.fc)){

gene.data <- Sortby(data.frame(melt(t(tissue.counts[cluster.names==clust.iter,])),pheno.factor),"Var1")

gee.out <- geeglm(value~pheno.factor+Var2,data=gene.data,id=Var1,family=poisson)

gee.sum <- summary(gee.out)$coef

results.fc[clust.iter,] <- c(exp(-gee.sum["pheno.factor2","Estimate"]),gee.sum["pheno.factor2","Pr(>|W|)"])

}

hist(results.fc[,"pvalue"],main="GEE p-values")
legend("topright",legend=paste(100*round(mean(results.fc[,"pvalue"]<0.05),2),"% P<0.05"),bg="white")

boxplot(results.fc[,"FC"]~de.info[,4][seq(1,n.clusters*p,by=p)],log="y",main="True vs. GEE Estimated FC")
abline(h=unique(de.info[,4]))

GEEseq.results <- aggregate(list(DE=results.fc[,"pvalue"]<0.02),list(FC=de.info[,4][seq(1,n.clusters*p,by=p)]),mean)



tex.list <- lcat(tex.list,latex.out.table.2(GEEseq.results,
	file=NULL,
	caption = paste("GEE results"),
	label = "",
	rownamer = NULL,
	col.label=NULL))



fit.gee.out <- fit.gee.2.array(tissue.counts,n.clusters,pheno.factor)


fit.gee.perm.out <- fit.gee.2.array.perm(tissue.counts,n.clusters,pheno.factor)

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






