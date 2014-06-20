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

source.file <- gsub("\\.R","/","check_gee_analysis_Jan_2014_data_update.R")

# define the base directory	

basedir <- "/Users/jonathangelfond/Documents/Projects/Buffenstein/NMR_Genomics/"

# Creates necessary directories

analysisdir <- paste(basedir,"Programs/",sep="") # where the programs are
datadir <- paste(basedir,"Data/Hibbs/countFiles/",sep="")  # where the data are
resultsdir <- paste(basedir,"Documents/Results/",source.file,sep="") # Standard output
tex.dir <- paste(resultsdir,"texdir/",sep="") # Publication quality output


mus2nmr.file <- "/Users/jonathangelfond/Documents/Projects/Buffenstein/NMR_Genomics/Data/Hibbs/Transcript_Mappings/mmu2nmr.hits.txt"
nmr2mus.file <- "/Users/jonathangelfond/Documents/Projects/Buffenstein/NMR_Genomics/Data/Hibbs/Transcript_Mappings/nmr2mmu.hits.txt"

newalignment.file <- "/Users/jonathangelfond/Documents/Projects/Buffenstein/NMR_Genomics/Data/Matt_homolgy/nmr2mmu.25frag.nwalign"

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
		
	return(p.out)
	
	})

#pvalueQQplot(pvals,title=paste(data.name,tissue,"p-values"))

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

nmr2mus.new <- read.delim(newalignment.file,as.is=TRUE)
# Columns in dataset per Matt H
#0 = NMR gene ID
#1 = Systematic mouse gene ID (from UCSC)
#2 = Gene symbol for mouse
#3 = NMR transcript length
#4 = mouse transcript length
#5 = length of best alignment
#6 = no. of mismatches in best alignment
#7 = homology score as percentage (between 0 and 1)
names(nmr2mus.new ) <- c("NMR.gene","mus.UCSCID","mus.Gene.Symbol","NMR.length","mus.length","alignment.length","n.mismatches","proportion.homology")

mus2nmr.matches <- subset(nmr2mus.new,mus.UCSCID!="UNALIGNED")
match.table <- table(mus2nmr.matches$mus.UCSCID)
temp <- hist(match.table,breaks=seq(0.5,max(match.table)+0.5,by=1),main="MUS2NMR 1 to Many Matches")

matplot(rbind(temp$mids,temp$mids),rbind(0,temp$counts),lty=1,col=1,type="l",xlab="Number of Matches",ylab="Count",main="Frequency vs. Number of Matches")
points(temp$mids,temp$counts,pch=16)
text(temp$mids,temp$counts-500,temp$counts,col=2)

match.count <- data.frame(mus.UCSCID=names(match.table),Count=match.table)

mus2nmr.matches.count <- merge(mus2nmr.matches,match.count,by="mus.UCSCID")

matches.2 <- mus2nmr.matches.count

matches.2 <- matches.2[order(matches.2$mus.UCSCID,matches.2$proportion.homology),]

matches.2$Streak.ct <- Streak.count(matches.2$mus.UCSCID)

melted.2 <- subset(matches.2,Count.Freq==2,select=c("proportion.homology","Streak.ct","mus.UCSCID"))

id.matrix <- dcast(melted.2, mus.UCSCID ~ Streak.ct  ,value.var="proportion.homology")

with(id.matrix,plot(`1`,`2`,pch=19))

id.matrix$ratio1over2 <- id.matrix$"1"/id.matrix$"2"

hist(-log2(id.matrix$ratio1over2),main="Ratio of Bigger over Smaller Match ID")

with(mus2nmr.matches.count,plot(Count.Freq,proportion.homology,pch=19))

ggplot2::qplot(`1`,`2`,data=id.matrix,alpha=I(1/10)) + xlab("Prop ID Match 1") + ylab("Prop ID Match 2") + labs(title="Average ID by Match") 

ggplot2::qplot(log10(`1`),log10(`2`),data=id.matrix,alpha=I(1/5),geom=c("point","smooth")) + xlab("Prop ID Match 1") + ylab("Prop ID Match 2") + labs(title="Average ID by Match (log10)") 



ave.id <- aggregate(subset(mus2nmr.matches.count,select="proportion.homology"),subset(mus2nmr.matches.count,select="Count.Freq"),mean,na.rm=TRUE)

matplot(rbind(ave.id$Count.Freq,ave.id$Count.Freq),rbind(0,ave.id$proportion.homology),lty=1,col=1,type="l",xlab="Number of Matches",ylab="Prop ID",log="x",main="Proportion ID vs. Number of Matches")


ggplot2::qplot(Count.Freq,log10(proportion.homology),data=mus2nmr.matches.count,geom=c("point","smooth")) + labs(title="Average ID by Match") + xlab("# Matches") + ylab("Prop ID")



## examine nmr to GP

head(all.d[["NMR2nmr"]])

expression.mat <- as.matrix(all.d[["NMR2nmr"]][,grep("(Liver)|(Skin)",names(all.d[["NMR2nmr"]]),value=TRUE)])
rownames(expression.mat) <- all.d[["NMR2nmr"]]$gene_id



sum(duplicated(rownames(expression.mat) ))

perm.exp <- expression.mat[sample(nrow(expression.mat)),]
rownames(perm.exp) <- rownames(perm.exp)[sample(nrow(expression.mat))]

mean(all.d[["NMR2nmr"]]$gene_id %in% matches.2$V2)

mean(matches.2$NMR.gene %in% all.d[["NMR2nmr"]]$gene_id)

length(unique(matches.2$NMR.gene))


cor.calc.out <- ddply(subset(matches.2,Count.Freq>1),"mus.UCSCID",compute.cluster.correlations.new,e.mat=expression.mat)


cor.calc.out.perm <- ddply(subset(matches.2,Count.Freq>1),"mus.UCSCID",compute.cluster.correlations.new,e.mat=perm.exp)

	
par(mfrow=c(1,2))

boxplot(mean.cor~Count.Freq,data=cor.calc.out,main="Cluster Correlation vs. Cluster Size")
abline(h=0,lwd=2,col=2)

boxplot(mean.cor~Count.Freq,data=cor.calc.out.perm,main="Cluster Correlation vs. Cluster Size (Permuted)")
abline(h=0,lwd=2,col=2)

cor.data <- rbind(data.frame(cor.calc.out,Perm="NotPerm"),data.frame(cor.calc.out.perm,Perm="Perm"))

ggplot2::qplot(Count.Freq,mean.cor,data=cor.calc.out,geom=c("point","smooth"),alpha=I(0.1)) + labs(title="Cluster Correlation vs. Cluster Size") + xlab("# Matches") + ylab("Correlation")

ggplot2::qplot(Count.Freq,mean.cor,data=cor.data,geom=c("point","smooth")) + labs(title="Cluster Correlation vs. Cluster Size") + xlab("# Matches") + ylab("Correlation") + facet_grid(.~Perm)


cor.merge <- merge(cor.calc.out,id.matrix,by="mus.UCSCID")

ggplot2::qplot(ratio1over2,mean.cor,data=cor.merge,geom=c("point","smooth")) + labs(title="Mean Correlation of Cluster vs. Ratio of Prop ID in Cluster Size = 2")


workspace.file <- file.path(resultsdir,"all.R.Obj.Rdata")

save(list=ls(),file=workspace.file)
#load(workspace.file)
source("/Users/jonathangelfond/Documents/Projects/Joe_Lab/NMR_Methods/Programs/support_functions.R")
library(geepack)



install.packages("gplots")
install.packages("blockmodeling")
install.packages("/Users/jonathangelfond/Documents/Projects/Buffenstein/NMR_Genomics/Documents/Background/EBSeq_1.1.6.tar.gz", repos=NULL, type="source")
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

tissue.counts.norm <- QuantileNorm(tissue.counts,0.75)

fit.gee.out <- fit.gee.5.array(tissue.counts,matches.2.1tomany,pheno.factor,crash.list=c(100:nrow(tissue.counts)))





#fit.gee.perm.out <- fit.gee.2.array.perm(tissue.counts,n.clusters,pheno.factor)

fit.merge <- merge(fit.gee.out,cor.calc.out,by="V1")

fit.merge$log10pvalue <- -log10(fit.merge$pvalue+1e-15)

hist(fit.merge$pvalue)

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


par(mfrow=c(1,1))

with(merge.comp,plot(PPDE,-log10(pvalue),main="Cluster mean PPDE (EBSeq) vs. -log10(pvalue) Cluster (GEE)"))

#with(merge.comp,cor.test(PPDE,-log10(pvalue+1e-15)))

hist(log2(pp.df$EB.FC),xlim=c(0,10))

with(merge.comp,plot(log2(FC),log.EB.FC,main="Cluster mean Fold-Change (EBSeq) vs. Fold-Change  Cluster (GEE)"))

with(merge.comp,cor.test(log2(FC),log.EB.FC))


print(paste("EOF:",source.file))






