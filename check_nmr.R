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


source.file <- gsub("\\.R","/","check_nmr.R")

# define the base directory	

basedir <- "/Users/jonathangelfond/Documents/Projects/Buffenstein/NMR_Genomics/"

# Creates necessary directories

analysisdir <- paste(basedir,"Programs/",sep="") # where the programs are
datadir <- paste(basedir,"Data/Hibbs/countFiles/",sep="")  # where the data are
resultsdir <- paste(basedir,"Documents/Results/",source.file,sep="") # Standard output
tex.dir <- paste(resultsdir,"texdir/",sep="") # Publication quality output


#loads libraries and adds common functions to workspace
source(paste(analysisdir,"support_functions.R",sep="")) 


dir.create(resultsdir)
dir.create(tex.dir)

tex.list <- list()

# Assign inpute filenames


d.files <- list.files(datadir)

all.d <- lapply(d.files,function(x){
	
	out <- read.delim(pastes(datadir,x),as.is=TRUE)
	
	})
names(all.d) <- gsub("\\..*","",d.files)



data.name <- "NMR2nmr"

pdf(pastes(resultsdir,"significance_tests.pdf"))

for(data.name in names(all.d)){


nmr.d <- all.d[[data.name]]

data.mat <- as.matrix(nmr.d[,grep("(Liver)|(Skin)",names(nmr.d))])

#rownames(data.mat) <- nmr.d[,grep("(ymbol)|(^gene)",names(nmr.d))[1]]

lapply(all.d,names)

# Multiply by effective length

data.mat.counts <- data.mat * nmr.d$effective_length
data.mat.counts <- round(data.mat * nmr.d$length,0)

#tissue <- "Liver"

for(tissue in c("Liver","Skin")){

liver.counts <- data.mat.counts[,grep(tissue,colnames(data.mat.counts))]
liver.counts <- log2(liver.counts[rowSums(liver.counts)>0,]+0.5)

pheno.data <- data.frame(treatment=grep.to.label(colnames(liver.counts),"C","Placebo","Drug"))
pheno.data$replicate <- c(1:3,1:2)


pvals <- apply(liver.counts,1,function(x){
	
	p.out <- NA
	
	try(p.out <- t.test(x~as.factor(pheno.data$treatment))$p.value)
	
	return(p.out)
	
	})

pvalueQQplot(pvals,title=paste(data.name,tissue,"p-values"))

} #loop over tissue
} 

dev.off()

#cds <- newCountDataSet(liver.counts,pheno.data$treatment)
#head( counts(cds) )
# cds <- estimateSizeFactors( cds )
# sizeFactors( cds )
#head( counts( cds, normalized=TRUE ) )
#cds <- estimateDispersions( cds )
#res <- nbinomTest( cds, "Placebo", "Drug" )


#table.files <- write.includer.vector(tex.list,path=tex.dir,include.file.base="includemetables")
#make.latex.doc.vector(pastes(tex.dir,"summary.tex"),includer=table.files,title="Analysis",author="Jon Gelfond",date=as.character(Sys.time()),path=NULL)


stop()

## examine nmr to GP

nmr.d <- all.d[["NMR2guineapig"]]

data.mat <- as.matrix(nmr.d[,grep("(Liver)",names(nmr.d))])

#rownames(data.mat) <- nmr.d[,grep("(ymbol)|(^gene)",names(nmr.d))[1]]

lapply(all.d,names)

# Multiply by effective length

data.mat.counts <- data.mat * nmr.d$effective_length
data.mat.counts <- round(data.mat * nmr.d$length,0)

#tissue <- "Liver"

tissue <- "(Liver)|(Skin)"

liver.counts <- data.mat.counts[,grep(tissue,colnames(data.mat.counts))]
liver.counts <- log2(liver.counts[rowSums(liver.counts)>0,]+0.5)

pheno.data <- data.frame(treatment=grep.to.label(colnames(liver.counts),"C","Placebo","Drug"))
pheno.data$replicate <- c(1:3,1:2)

vars <- apply(data.mat,1,var,na.rm=TRUE)
means <- apply(data.mat,1,mean,na.rm=TRUE)

plot(means,vars,log="xy")

plot(nmr.d$X..Identity.with.respect.to.Mouse.gene,vars,log="y")

log.vars <- log10(vars)

log.vars[is.na(log.vars)|(log.vars==-Inf)] <- NA

pct.id <- nmr.d$X..Identity.with.respect.to.Mouse.gene

pct.id <- as.numeric(as.character(pct.id))

pct.id[is.na(pct.id)] <- 0

log.means <- log10(means)

xmat <- na.exclude(data.frame(log.vars,log.means,pct.id,x= poly(pct.id,7)))

xmat <- subset(xmat,pct.id!=0)

xmat <- xmat[order(xmat$pct.id),]

temp <- lm(log.vars ~ x.1 + x.2 + x.3+x.4+x.5+x.6+x.7,data=xmat)

plot(pct.id,log.vars,pch=16,main=)

pred <- predict(temp,newdata=xmat,interval="predict")

lines(xmat$pct.id,pred[,1],col=2,lty=2,lwd=3)

temp <- lm(log.vars ~ pct.id*log.means,data=xmat)

plot(log.means,log.vars)

test.pct <- c(5,95)

for(pct.iter in test.pct){

xmat2 <- xmat
xmat2$pct.id <- pct.iter

pred <- predict(temp,newdata=xmat2,interval="predict")
points(xmat2$log.means,pred[,1],col=which(test.pct==pct.iter)+1,pch=".")

}
legend("topleft",legend=c("5%ID","95%ID"),lty=1,col=c(2,3))


mean.slice <- subset(xmat,(log.means>log10(100))&(log.means<log10(500)))

mean.slice$pct.id.factor <- cut(mean.slice$pct.id,breaks=quantile(mean.slice$pct.id,c(0,.25,.5,.75,1.00)))



boxplot(log.vars~pct.id.factor,data=mean.slice,xlab="% ID Quartile",ylab="log10 Variance",main="Mean counts [100, 500]")


mean.slice$pct.id.factor <- cut(mean.slice$pct.id,breaks=seq(0,100,by=20))

boxplot(log.vars~pct.id.factor,data=mean.slice,xlab="% ID Factor",ylab="log10 Variance",main="Mean counts [100, 500]")

library(ggplot2)

qplot(mean.slice$pct.id,mean.slice$log.vars,geom=c("point","smooth"),formula=y~ns(x,5))


temp <- lm(log.means ~ x.1 + x.2 + x.3+x.4+x.5+x.6+x.7,data=xmat)

plot(pct.id,log.means,pch=16,main=)

pred <- predict(temp,newdata=xmat,interval="predict")

lines(xmat$pct.id,pred[,1],col=2,lty=2,lwd=3)




install.packages("gplots")
install.packages("blockmodeling")
install.packages("/Users/jonathangelfond/Documents/Projects/Buffenstein/NMR_Genomics/Documents/Background/EBSeq_1.1.6.tar.gz", repos=NULL, type="source")
library(gplots)
library(blockmodeling)
library(EBSeq)



#demo(EBSeq)


data(GeneMat)


gene.name.col <- cbind(c("Symbol","MGI.symbol","Symbol","Symbol","gene_id"))
rownames(gene.name.col) <- names(all.d)

ebseq.results <- list()

pdf(pastes(resultsdir,"eb_seq_significance_tests.pdf"))

for(data.name in names(all.d)){


nmr.d <- all.d[[data.name]]

data.mat <- as.matrix(nmr.d[,grep("(Liver)|(Skin)",names(nmr.d))])

rownames(data.mat) <- nmr.d[,gene.name.col[data.name,]]


# Multiply by effective length

data.mat.counts <- data.mat * nmr.d$effective_length
data.mat.counts <- round(data.mat * nmr.d$length,0)

#tissue <- "Liver"

for(tissue in c("Liver","Skin")){

try({

tissue.counts <- data.mat.counts[,grep(tissue,colnames(data.mat.counts))]

tissue.counts <- tissue.counts[!duplicated(rownames(tissue.counts)),]
tissue.counts <- tissue.counts[rownames(tissue.counts)!="",]


Sizes=MedianNorm(tissue.counts)


pheno.data <- data.frame(treatment=as.factor(grep.to.label(colnames(tissue.counts),"C","Placebo","Drug")))
pheno.data$replicate <- c(1:3,1:2)

EBOut <- EBTest(Data=tissue.counts,Conditions=pheno.data$treatment,sizeFactors=Sizes, maxround=5)

PP <- GetPPMat(EBOut)
post.fc <- PostFC(EBOut)

mean.df <- data.frame(gene=names(EBOut$C1Mean[[1]]),gene2=names(EBOut$C2Mean[[1]]),Drug.mean=EBOut$C1Mean[[1]],Placebo.mean=EBOut$C2Mean[[1]])

if(identical(mean.df$gene,mean.df$gene2)){mean.df <- Drop.variable(mean.df,"gene2")}

post.fc.df <- data.frame(gene=names(post.fc$PostFC),post.fc$PostFC)


tissue.df <- data.frame(gene=rownames(tissue.counts),tissue.counts)


PPsort <- merge(data.frame(gene=rownames(PP),PP),post.fc.df,by="gene")

PPsort <- merge(PPsort,mean.df,by="gene")

PPsort <- merge(PPsort,tissue.df,by="gene")


PPsort <- Sortby(PPsort,"PPEE")

DEfound=rownames(PP)[which(PP[,"PPDE"]>=.95)]

ebseq.results[[paste(data.name,tissue)]] <- PPsort

write.csv(PPsort,paste(resultsdir,tissue,data.name,".csv",sep="_"))

hist(PP[,"PPDE"],xlab="Posterior Probability DE",main=paste(data.name,tissue))
legend("topright",legend=paste(100*round(prop.table(table(PP[,"PPDE"]>0.5)),2),"%",c("Not DE","DE")))

})

} #loop over tissue
} #loop over dataset 

dev.off()




# Look at human liver




data.name <- "NMR2guineapig"

nmr.d <- all.d[[data.name]]

data.mat <- as.matrix(nmr.d[,grep("(Liver)|(Skin)",names(nmr.d))])

rownames(data.mat) <- nmr.d[,gene.name.col[data.name,]]


# Multiply by effective length

data.mat.counts <- data.mat * nmr.d$effective_length
data.mat.counts <- round(data.mat * nmr.d$length,0)

#tissue <- "Liver"

tissue <- "Liver"


otholog.types <- unique(nmr.d$Homology.Type[duplicated(nmr.d$Homology.Type)])

pct.ids <- as.numeric(as.character(nmr.d$X..Identity.with.respect.to.Mouse.gene))

tissue.counts <- data.mat.counts[,grep(tissue,colnames(data.mat.counts))]


pct.ids.nodups <- pct.ids[!duplicated(rownames(tissue.counts))]

pct.ids.nodups.factor <- cut(pct.ids.nodups,breaks=quantile(pct.ids.nodups,seq(0,1,length=10),na.rm=TRUE))

tissue.counts.nodups <- tissue.counts[!duplicated(rownames(tissue.counts)),]

pct.id.analysis.df <- data.frame(pct.id=NULL,alpha=NULL,beta=NULL,pct.de=NULL,median.var=NULL)

pheno.data <- data.frame(treatment=as.factor(grep.to.label(colnames(tissue.counts.no.dups.select),"C","Placebo","Drug")))
pheno.data$replicate <- c(1:3,1:2)


for(pct.iter in levels(pct.ids.nodups.factor)){

try({

select.condition <- pct.ids.nodups.factor==pct.iter
select.condition <- ifelse(is.na(select.condition),FALSE,select.condition)

tissue.counts.nodups.select <- tissue.counts.nodups[select.condition,]

Sizes=MedianNorm(tissue.counts.nodups.select)

itermax <- 5

EBOut <- EBTest(Data=tissue.counts.nodups.select,Conditions=pheno.data$treatment,sizeFactors=Sizes, maxround=itermax)


PP <- GetPPMat(EBOut)
post.fc <- PostFC(EBOut)

mean.df <- data.frame(gene=names(EBOut$C1Mean[[1]]),gene2=names(EBOut$C2Mean[[1]]),Drug.mean=EBOut$C1Mean[[1]],Placebo.mean=EBOut$C2Mean[[1]])

if(identical(mean.df$gene,mean.df$gene2)){mean.df <- Drop.variable(mean.df,"gene2")}

post.fc.df <- data.frame(gene=names(post.fc$PostFC),post.fc$PostFC)


tissue.df <- data.frame(gene=rownames(tissue.counts),tissue.counts)


PPsort <- merge(data.frame(gene=rownames(PP),PP),post.fc.df,by="gene")

PPsort <- merge(PPsort,mean.df,by="gene")

PPsort <- merge(PPsort,tissue.df,by="gene")


PPsort <- Sortby(PPsort,"PPEE")

DEfound=rownames(PP)[which(PP[,"PPDE"]>=.95)]

#write.csv(PPsort,paste(resultsdir,tissue,data.name,".csv",sep="_"))

hist(PP[,"PPDE"],xlab="Posterior Probability DE",main=paste(data.name,tissue))
legend("topright",legend=paste(100*round(prop.table(table(PP[,"PPDE"]>0.5)),2),"%",c("Not DE","DE")))

ave.mean <- 0.5*(EBOut$C1Mean[[1]]+EBOut$C2Mean[[1]])



pct.id.analysis.df <- rbind(pct.id.analysis.df,data.frame(pct.id=pct.iter,alpha=EBOut$Alpha[itermax],beta=EBOut$Beta[itermax],pct.de=mean(PP[,"PPDE"]>0.5,na.rm=TRUE),median.var=median(EBOut$PoolVar[[1]]/ave.mean,na.rm=TRUE)))

})

} #loop over pct.iter


pdf(pastes(resultsdir,"pct_id_plots.pdf"),h=8,w=11)


with(pct.id.analysis.df,plot(as.numeric(pct.id),alpha,xlab="Pct ID Quantile",ylab="Alpha"))
with(pct.id.analysis.df,plot(as.numeric(pct.id),beta,xlab="Pct ID Quantile",ylab="Beta"))
with(pct.id.analysis.df,plot(alpha,beta,xlab="Alpha",ylab="Beta"))
with(pct.id.analysis.df,text(alpha,beta,pct.id))

with(pct.id.analysis.df,plot(as.numeric(pct.id),median.var,xlab="Pct ID",ylab="Median Coef of Variance"))

dev.off()



# Look at mapping id



data.name <- "NMR2guineapig"

nmr.d <- all.d[[data.name]]

data.mat <- as.matrix(nmr.d[,grep("(Liver)|(Skin)",names(nmr.d))])

rownames(data.mat) <- nmr.d[,gene.name.col[data.name,]]


# Multiply by effective length

data.mat.counts <- data.mat * nmr.d$effective_length
data.mat.counts <- round(data.mat * nmr.d$length,0)

#tissue <- "Liver"

for(tissue in c("Liver","Skin")){


otholog.types <- unique(nmr.d$Homology.Type[duplicated(nmr.d$Homology.Type)])


tissue.counts <- data.mat.counts[,grep(tissue,colnames(data.mat.counts))]


homology.types.nodups <- nmr.d$Homology.Type[!duplicated(rownames(tissue.counts))]

tissue.counts.nodups <- tissue.counts[!duplicated(rownames(tissue.counts)),]

pct.id.analysis.df <- data.frame(ortho.type=NULL,alpha=NULL,beta=NULL,pct.de=NULL,median.var=NULL)

by.gene.analysis.df <- data.frame(ortho.type=NULL,mean=NULL,variance=NULL)

pheno.data <- data.frame(treatment=as.factor(grep.to.label(colnames(tissue.counts.no.dups.select),"C","Placebo","Drug")))
pheno.data$replicate <- c(1:3,1:2)


for(ortho.iter in otholog.types){

try({

select.condition <- (homology.types.nodups==ortho.iter)&(rownames(tissue.counts.nodups)!="")
select.condition <- ifelse(is.na(select.condition),FALSE,select.condition)

tissue.counts.nodups.select <- tissue.counts.nodups[select.condition,]

Sizes=MedianNorm(tissue.counts.nodups.select)

itermax <- 5


EBOut <- EBTest(Data=tissue.counts.nodups.select,Conditions=pheno.data$treatment,sizeFactors=Sizes, maxround=itermax)


PP <- GetPPMat(EBOut)
post.fc <- PostFC(EBOut)

mean.df <- data.frame(gene=names(EBOut$C1Mean[[1]]),gene2=names(EBOut$C2Mean[[1]]),Drug.mean=EBOut$C1Mean[[1]],Placebo.mean=EBOut$C2Mean[[1]])

if(identical(mean.df$gene,mean.df$gene2)){mean.df <- Drop.variable(mean.df,"gene2")}

post.fc.df <- data.frame(gene=names(post.fc$PostFC),post.fc$PostFC)


tissue.df <- data.frame(gene=rownames(tissue.counts),tissue.counts)


PPsort <- merge(data.frame(gene=rownames(PP),PP),post.fc.df,by="gene")

PPsort <- merge(PPsort,mean.df,by="gene")

PPsort <- merge(PPsort,tissue.df,by="gene")


PPsort <- Sortby(PPsort,"PPEE")

DEfound=rownames(PP)[which(PP[,"PPDE"]>=.95)]

#write.csv(PPsort,paste(resultsdir,tissue,data.name,".csv",sep="_"))

hist(PP[,"PPDE"],xlab="Posterior Probability DE",main=paste(data.name,tissue))
legend("topright",legend=paste(100*round(prop.table(table(PP[,"PPDE"]>0.5)),2),"%",c("Not DE","DE")))

ave.mean <- 0.5*(EBOut$C1Mean[[1]]+EBOut$C2Mean[[1]])



pct.id.analysis.df <- rbind(pct.id.analysis.df,data.frame(ortho.type=ortho.iter,alpha=EBOut$Alpha[itermax],beta=EBOut$Beta[itermax],pct.de=mean(PP[,"PPDE"]>0.5,na.rm=TRUE),median.var=median(sqrt(EBOut$PoolVar[[1]])/ave.mean,na.rm=TRUE)))


by.gene.analysis.df <- rbind(by.gene.analysis.df,data.frame(ortho.type=ortho.iter,mean=ave.mean,variance=EBOut$PoolVar[[1]]))



})

} #loop over pct.iter


pdf(paste(resultsdir,tissue,"_ortho_id_plots.pdf",sep=""),h=8,w=11)


with(pct.id.analysis.df,plot(ortho.type,alpha,xlab="Homology Type",ylab="Alpha",main=paste(tissue,"by Orthology")))
with(pct.id.analysis.df,plot(ortho.type,beta,xlab="Homology Type",ylab="Beta",main=paste(tissue,"by Orthology")))
with(pct.id.analysis.df,plot(alpha,beta,xlab="Alpha",ylab="Beta",main=paste(tissue,"by Orthology")))
with(pct.id.analysis.df,text(alpha,beta,ortho.type))

box.stats <- with(by.gene.analysis.df,boxplot(sqrt(variance)/mean~ortho.type,xlab="Homology Type",ylab="Median Coef of Variation",main=paste(tissue,"by Orthology")))

mtext(pastes("n=",box.stats$n),side=1,line=-1,at=1:length(box.stats$n))

with(by.gene.analysis.df,plot(ortho.type,sqrt(variance)/mean,xlab="Homology Type",ylab="Median Coef of Variation",main=paste(tissue,"by Orthology")))


dev.off()

} #loop over tissue







# Look at mapping pct id



data.name <- "NMR2guineapig"

nmr.d <- all.d[[data.name]]

data.mat <- as.matrix(nmr.d[,grep("(Liver)|(Skin)",names(nmr.d))])

rownames(data.mat) <- nmr.d[,gene.name.col[data.name,]]


# Multiply by effective length

data.mat.counts <- data.mat * nmr.d$effective_length
data.mat.counts <- round(data.mat * nmr.d$length,0)

#tissue <- "Liver"

for(tissue in c("Liver","Skin")){


otholog.types <- unique(nmr.d$Homology.Type[duplicated(nmr.d$Homology.Type)])


tissue.counts <- data.mat.counts[,grep(tissue,colnames(data.mat.counts))]


homology.types.nodups <- nmr.d$Homology.Type[!duplicated(rownames(tissue.counts))]


pct.ids.nodups <- pct.ids[!duplicated(rownames(tissue.counts))]

pct.ids.nodups.factor <- cut(pct.ids.nodups,breaks=quantile(pct.ids.nodups,seq(0,1,length=10),na.rm=TRUE))


tissue.counts.nodups <- tissue.counts[!duplicated(rownames(tissue.counts)),]

pct.id.analysis.df <- data.frame(pct.id=NULL,alpha=NULL,beta=NULL,pct.de=NULL,median.var=NULL)

by.gene.analysis.df <- data.frame(pct.id=NULL,mean=NULL,variance=NULL)

pheno.data <- data.frame(treatment=as.factor(grep.to.label(colnames(tissue.counts.no.dups.select),"C","Placebo","Drug")))
pheno.data$replicate <- c(1:3,1:2)


for(pct.id.iter in levels(pct.ids.nodups.factor)){

try({

select.condition <- (pct.ids.nodups.factor==pct.id.iter)&(rownames(tissue.counts.nodups)!="")
select.condition <- ifelse(is.na(select.condition),FALSE,select.condition)

tissue.counts.nodups.select <- tissue.counts.nodups[select.condition,]

Sizes=MedianNorm(tissue.counts.nodups.select)

itermax <- 5


EBOut <- EBTest(Data=tissue.counts.nodups.select,Conditions=pheno.data$treatment,sizeFactors=Sizes, maxround=itermax)


PP <- GetPPMat(EBOut)
post.fc <- PostFC(EBOut)

mean.df <- data.frame(gene=names(EBOut$C1Mean[[1]]),gene2=names(EBOut$C2Mean[[1]]),Drug.mean=EBOut$C1Mean[[1]],Placebo.mean=EBOut$C2Mean[[1]])

if(identical(mean.df$gene,mean.df$gene2)){mean.df <- Drop.variable(mean.df,"gene2")}

post.fc.df <- data.frame(gene=names(post.fc$PostFC),post.fc$PostFC)


tissue.df <- data.frame(gene=rownames(tissue.counts),tissue.counts)


PPsort <- merge(data.frame(gene=rownames(PP),PP),post.fc.df,by="gene")

PPsort <- merge(PPsort,mean.df,by="gene")

PPsort <- merge(PPsort,tissue.df,by="gene")


PPsort <- Sortby(PPsort,"PPEE")

DEfound=rownames(PP)[which(PP[,"PPDE"]>=.95)]

#write.csv(PPsort,paste(resultsdir,tissue,data.name,".csv",sep="_"))

hist(PP[,"PPDE"],xlab="Posterior Probability DE",main=paste(data.name,tissue))
legend("topright",legend=paste(100*round(prop.table(table(PP[,"PPDE"]>0.5)),2),"%",c("Not DE","DE")))

ave.mean <- 0.5*(EBOut$C1Mean[[1]]+EBOut$C2Mean[[1]])



pct.id.analysis.df <- rbind(pct.id.analysis.df,data.frame(ortho.type=pct.id.iter,alpha=EBOut$Alpha[itermax],beta=EBOut$Beta[itermax],pct.de=mean(PP[,"PPDE"]>0.5,na.rm=TRUE),median.var=median(sqrt(EBOut$PoolVar[[1]])/ave.mean,na.rm=TRUE)))


by.gene.analysis.df <- rbind(by.gene.analysis.df,data.frame(ortho.type=pct.id.iter,mean=ave.mean,variance=EBOut$PoolVar[[1]]))



})

} #loop over pct.iter


pdf(paste(resultsdir,tissue,"_pct_id_plots.pdf",sep=""),h=8,w=11)


with(pct.id.analysis.df,plot(ortho.type,alpha,xlab="Pct ID",ylab="Alpha",main=paste(tissue,"by Pct ID")))
with(pct.id.analysis.df,plot(ortho.type,beta,xlab="Pct ID",ylab="Beta",main=paste(tissue,"by Pct ID")))
with(pct.id.analysis.df,plot(alpha,beta,xlab="Alpha",ylab="Beta",main=paste(tissue,"by Pct ID")))
with(pct.id.analysis.df,text(alpha,beta,ortho.type))

box.stats <- with(by.gene.analysis.df,boxplot(sqrt(variance)/mean~ortho.type,xlab="Pct ID",ylab="Median Coef of Variation",main=paste(tissue,"by Orthology")))

mtext(pastes("n=",box.stats$n),side=1,line=-1,at=1:length(box.stats$n))



dev.off()

} #loop over tissue



temp <- Sortby(subset(nmr.d,Homology.Type=="ortholog_one2many"),"MGI.ID")


temp





print(paste("EOF:",source.file))






