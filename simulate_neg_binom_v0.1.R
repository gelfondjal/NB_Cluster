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


set.seed(2013)
rm(list=ls())


source.file <- gsub("\\.R","/","simulate_neg_binom_v0.1.R")

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

simulate.2n.arrays <- function(n.arrays,delta,p.delta,p,sigmat.mat,mu,nb.size){
	
	# simulate 2*n.arrays of 2-condiction differential expression n.arrays on right & n.arrays on left
	# n.arrays = # of arrays in 1 condition
	# delta = effect size if differentiall expressed (DE)
	# p.delta = probability of DE
	# sigma.mat = correlation between nb elements of length p
	# n.clusters = number of clusters
	# mu = average gene expression counts
	# nb.size = NB parameter size
	
	de.indicator <- rep(rbinom(n.clusters,1,p.delta),each=1)
	flipper <- rep(sample(c(-1,1),n.clusters,replace=TRUE),each=1)
	flipper <- ifelse(de.indicator==0,0,flipper)

	shifter <- (de.indicator*delta)^flipper
	
	
	left.right <- rep(rbinom(n.clusters,1,0.5),each=1)
	
	shifter.left <- ifelse(left.right==0,shifter,1)
	shifter.right <- ifelse(left.right==1,shifter,1)
	
	#fill out left matrix
	
	mat.out <- matrix(NA,n.clusters*p,n.arrays*2)
	
	rownames(mat.out) <- paste("homolog",rep(1:n.clusters,each=p),rep(1:p,n.clusters),sep=".")
	
	
	# stack correlated variables on top of one another
	
	for(i in 1:n.arrays){
		mat.out[,i] <- as.vector(t(r.cor.nb(n.clusters,p,sigma.mat,mu*shifter.left,nb.size)))
		mat.out[,i+n.arrays] <- as.vector(t(r.cor.nb(n.clusters,p,sigma.mat,mu*shifter.right,nb.size)))
	}
	
	return(list(array=mat.out,shifter=cbind(rep(shifter.left,each=p),rep(shifter.right,each=p),rep(shifter,each=p))))
	
	}


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





qplot(log10(results.fc[,"FC"]),1*(results.fc[,"pvalue"]<0.05),geom=c("point","smooth"),formula=y~ns(x,15))


cor(cor.nb.out)

colMeans(cor.nb.out)

plot((cor.nb.out))

plot(log(cor.nb.out))





#tex.list <- lcat(tex.list,latex.out.table.2(numberOfPeptides,
#	file=NULL,
#	caption = paste("Number Peptides"),
#	label = "",
#	rownamer = NULL,
#	col.label=NULL))


pdf(cor.nb.file <- pastes(resultsdir,"cor_nb_figs.pdf"))

par(mfrow=c(2,1))

plot((cor.nb.out),main="Correlated Negative Binomial")

plot(log(cor.nb.out),main="Log Correlated Negative Binomial")


dev.off()

tex.list <- lcat(tex.list,latex.figure.out(cor.nb.file,figtype="htpb",width="17cm",height="15cm",caption="Plots of 2D Correlated negative binomial",label=""))







table.files <- write.includer.vector(tex.list,path=tex.dir,include.file.base="includemetables")
make.latex.doc.vector(pastes(tex.dir,"summary.tex"),includer=table.files,title="Simulated Correlated Negative binomial Analysis",author="Jon Gelfond",date=as.character(Sys.time()),path=NULL)


print(paste("EOF:",source.file))






