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
#source("http://bioinformatics.mdanderson.org/OOMPA/oompaLite.R")
#oompaLite()
library(BMA)

set.seed(2013)
rm(list=ls())


source.file <- gsub("\\.R","/","simulate_neg_binom_v0.quick.R")

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

n.clusters <- 100
p.delta <- 0.2


array.data.out <- simulate.2n.arrays(n.arrays,delta,p.delta,p,sigmat.mat,mu,nb.size)

tissue.counts <- array.data.out$array
de.info <- array.data.out$shifter


de.info <- cbind(de.info,de.info[,1]/de.info[,2])

pheno.factor <- as.factor(rep(c(1,2),each=n.arrays))

colnames(tissue.counts) <- paste("c",pheno.factor,rep(1:n.arrays,2),sep=".")
#heatmap.2(log(tissue.counts+0.5),col=redblue(256),trace="none",title="Simulated Data")

hist(log10(tissue.counts))

vars <- apply(tissue.counts,1,var)
means <- apply(tissue.counts,1,mean)

plot(log(means),log(vars))

library(gplots)
library(blockmodeling)
library(EBSeq)

Sizes=MedianNorm(tissue.counts)

itermax <- 5



match.info <- data.frame(V1=gsub("\\.[0-9]+$","",rownames(tissue.counts)),V2=rownames(tissue.counts),stringsAsFactors=FALSE)

n.transcripts <- p*n.clusters

families <- list(poisson=poisson,quasipoisson=quasipoisson,nbinom1="nbinom1",MGLM="MGLM",gee="GEE")
system.time({
analysis.out <- lapply(families,function(x){fit.gee.5.array(tissue.counts,match.info,pheno.factor,crash.list=1:n.transcripts,family=x)})
analysis.out.perm <- lapply(families,function(x){fit.gee.5.array(tissue.counts,match.info,pheno.factor,crash.list=1:n.transcripts,family=x,perm=4)})
})

perm.combo <- list()
for(namer in names(analysis.out)){perm.combo[[namer]] <- list(orig=analysis.out[[namer]][[1]],perm=analysis.out.perm[[namer]][[1]])}

corrected.ps <- lapply(perm.combo,function(x){return(p.bum.correct.right.truncate(x$orig[,"pvalue"],x$perm[,"pvalue"]))})


pdf(paste0(resultsdir,"comparing_simulated_count_models_prod.pdf"),height=8,width=11)


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

df.block <- rbind.fill(data.frame(corrected.ps))

corrected.df0 <- melt(data.frame(rowid=1:nrow(df.block),df.block),id.vars="rowid")

corrected.df <- merge(corrected.df0,data.frame(rowid=1:n.transcripts,FC=de.info[1:n.transcripts,4]),by="rowid")

corrected.df <- rename(corrected.df,replace=c("variable"="model","value"="p.value"))
alpha <- 0.05
sum.power <- ddply(corrected.df,c("FC","model"),function(x){
	
	y.out <- data.frame("Mean.sig"=mean(x$p.value<alpha,na.rm=TRUE),"Sum.sig"=sum(x$p.value<alpha,na.rm=TRUE),"out.of.Total"=nrow(x))
	
	return(y.out)
	
	})

cor.table <- cor(as.matrix(rbind.fill(data.frame(corrected.ps))),use="pair")



tex.list <- lcat(tex.list,latex.out.table.2(round(cor.table,2),caption="Correlation of -log10(pvalues) for various methods"))




tex.list <- lcat(tex.list,latex.out.table.2(Norow.mat(round.data(sum.power,2)),caption=paste("Operating characteristics at alpha=",alpha)))




table.files <- write.includer.vector(tex.list,path=tex.dir,include.file.base="includemetables")
make.latex.doc.vector(pastes(tex.dir,"summary.tex"),includer=table.files,title="Simulated Correlated Negative binomial Analysis",author="Jon Gelfond",date=as.character(Sys.time()),path=NULL)

Create.pdflatex("summary.tex",tex.dir)

print(paste("EOF:",source.file))






