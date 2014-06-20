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
library(ggplot2)
library(DESeq)
library(reshape2)

source.file <- gsub("\\.R","/","check_nmr_homologs.R")

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


nmr2mus <- read.table(nmr2mus.file,as.is=TRUE)
mus2nmr <- read.table(mus2nmr.file,as.is=TRUE)

mus2nmr.match <- subset(mus2nmr,ifelse(V2=="UNALIGNED",V3<1,TRUE))

mus2nmr.matches <- subset(mus2nmr,V2!="UNALIGNED")
match.table <- table(mus2nmr.matches$V1)
temp <- hist(match.table,breaks=seq(0.5,max(match.table)+0.5,by=1),main="MUS2NMR 1 to Many Matches")

matplot(rbind(temp$mids,temp$mids),rbind(0,temp$counts),lty=1,col=1,type="l",xlab="Number of Matches",ylab="Count",main="Frequency vs. Number of Matches")
points(temp$mids,temp$counts,pch=16)
text(temp$mids,temp$counts-500,temp$counts,col=2)

match.count <- data.frame(V1=names(match.table),Count=match.table)

mus2nmr.matches.count <- merge(mus2nmr.matches,match.count,by="V1")

matches.2 <- mus2nmr.matches.count

matches.2 <- matches.2[order(matches.2$V1,matches.2$V3),]

matches.2$Streak.ct <- Streak.count(matches.2$V1)

melted.2 <- subset(matches.2,Count.Freq==2,select=c("V3","Streak.ct","V1"))

id.matrix <- dcast(melted.2, V1 ~ Streak.ct  ,value.var="V3")

with(id.matrix,plot(`1`,`2`,pch=19))

id.matrix$ratio1over2 <- id.matrix$"1"/id.matrix$"2"

hist(-log2(id.matrix$ratio1over2),main="Ratio of Bigger over Smaller Match ID")

with(mus2nmr.matches.count,plot(Count.Freq,V3,pch=19))

qplot(`1`,`2`,data=id.matrix,alpha=I(1/10)) + xlab("Prop ID Match 1") + ylab("Prop ID Match 2") + labs(title="Average ID by Match") 

qplot(log10(`1`),log10(`2`),data=id.matrix,alpha=I(1/5),geom=c("point","smooth")) + xlab("Prop ID Match 1") + ylab("Prop ID Match 2") + labs(title="Average ID by Match (log10)") 



ave.id <- aggregate(subset(mus2nmr.matches.count,select="V3"),subset(mus2nmr.matches.count,select="Count.Freq"),mean,na.rm=TRUE)

matplot(rbind(ave.id$Count.Freq,ave.id$Count.Freq),rbind(0,ave.id$V3),lty=1,col=1,type="l",xlab="Number of Matches",ylab="Prop ID",main="Proportion ID vs. Number of Matches")


qplot(Count.Freq,log10(V3),data=mus2nmr.matches.count,geom=c("point","smooth")) + labs(title="Average ID by Match") + xlab("# Matches") + ylab("Prop ID")



## examine nmr to GP

head(all.d[["NMR2nmr"]])

expression.mat <- as.matrix(all.d[["NMR2nmr"]][,grep("(Liver)|(Skin)",names(all.d[["NMR2nmr"]]),value=TRUE)])
rownames(expression.mat) <- all.d[["NMR2nmr"]]$gene_id

perm.exp <- expression.mat[sample(nrow(expression.mat)),]
rownames(perm.exp) <- rownames(perm.exp)[sample(nrow(expression.mat))]

mean(all.d[["NMR2nmr"]]$gene_id %in% matches.2$V2)

mean(matches.2$V2 %in% all.d[["NMR2nmr"]]$gene_id)

length(unique(matches.2$V2))


cor.calc.out <- ddply(subset(matches.2,Count.Freq>1),"V1",compute.cluster.correlations,e.mat=expression.mat)


cor.calc.out.perm <- ddply(subset(matches.2,Count.Freq>1),"V1",compute.cluster.correlations,e.mat=perm.exp)

	
par(mfrow=c(1,2))

boxplot(mean.cor~Count.Freq,data=cor.calc.out,main="Cluster Correlation vs. Cluster Size")
abline(h=0,lwd=2,col=2)

boxplot(mean.cor~Count.Freq,data=cor.calc.out.perm,main="Cluster Correlation vs. Cluster Size (Permuted)")
abline(h=0,lwd=2,col=2)

cor.data <- rbind(data.frame(cor.calc.out,Perm="NotPerm"),data.frame(cor.calc.out.perm,Perm="Perm"))

qplot(Count.Freq,mean.cor,data=cor.calc.out,geom=c("point","smooth"),alpha=I(0.1)) + labs(title="Cluster Correlation vs. Cluster Size") + xlab("# Matches") + ylab("Correlation")

qplot(Count.Freq,mean.cor,data=cor.data,geom=c("point","smooth")) + labs(title="Cluster Correlation vs. Cluster Size") + xlab("# Matches") + ylab("Correlation") + facet_grid(.~Perm)


cor.merge <- merge(cor.calc.out,id.matrix,by="V1")

qplot(ratio1over2,mean.cor,data=cor.merge,geom=c("point","smooth")) + labs(title="Mean Correlation of Cluster vs. Ratio of Prop ID in Cluster Size = 2")



source("/Users/jonathangelfond/Documents/Projects/Joe_Lab/NMR_Methods/Programs/support_functions.R")


tissue.counts <- expression.mat

pheno.factor <- as.factor(grep.to.label(colnames(tissue.counts),"Liver","Liver","Skin"))

matches.2.1tomany <- subset(matches.2,Count.Freq>1)

tissue.counts <- tissue.counts[matches.2.1tomany$V2,]


fit.gee.out <- fit.gee.3.array(tissue.counts,matches.2.1tomany$V1,pheno.factor)
#fit.gee.perm.out <- fit.gee.2.array.perm(tissue.counts,n.clusters,pheno.factor)





print(paste("EOF:",source.file))






