##############################################################################
# University of Texas Health Science Center at San Antonio
# Department of Epidemiology and Biostatistics                        
##############################################################################
# Filename: analyze_v0.3.R
# Author:              Jonathan Gelfond                                   
# Project Name:          NMR RNA Seq
# Input:                  Codebook Scoresheet.csv    & Code Book Data Sheet.csv
# Output: 
#
# Modification History:
# v 0.1 Creation           
##############################################################################



library(samr)
library(impute)
library(lars)
library(leaps)
library(qvalue)
library(xtable)
library(multtest)
#library(sva)
library(devEMF)
library(psych)
library(gplots)
library(binom)
library(lme4)
library(splines)
library(Hmisc)
#library(Design)
library(BMA)
library(survival)
library(network)
library(glasso)
#library(nlme)


set.seed(2012)

rm(list=ls())

sourceFile <- "analyze_skin_v0.1.R"
sourceFileClip <- gsub("\\.R","/",sourceFile)

base1 <- "/Users/jonathangelfond/Documents/Projects/Buffenstein/NMR_Genomics/"
#base2 <- "/Users/jonathangelfond/Dropbox/Temp_Work/Weiner/"

basedir <- base1

analysisdir <- paste(basedir,"Programs/",sep="") 

datadir <- paste(basedir,"Data/",sep="")

resultsdir <- paste(basedir,"Documents/Results/",sourceFileClip,sep="")
resultsdir.r2sas <- paste(resultsdir,"R2SAS/",sep="")
resultsdir.sas <- paste("//psf/Home/Dropbox/OSCE_Project2/Jons_Template/Documents/Results/",sourceFileClip,"R2SAS/",sep="")


source("/Users/jonathangelfond/Documents/Projects/sudheer/Analysis/table_functions.R")

loudenPath <- "/Users/jonathangelfond/Documents/Projects/Louden/"
sapply(paste(loudenPath,list.files(loudenPath),sep=""),source)
source(pastes(analysisdir,"support_functions.R"))


dir.create(resultsdir)
dir.create(resultsdir.r2sas)


important.genes <- read.xls(pastes(datadir,"RNASeq Genes.xlsx"))
important.genes$SYMBOL <- toupper(important.genes$Official.Symbol)


# Read in Sequence Data

nmr.skin.data <- read.delim(paste(datadir,"NMRskin_DESeq_result/NMRskin_stat_symbol.xls",sep=""),as.is=TRUE,sep="\t")

sum(duplicated(nmr.skin.data$id))

mouse.skin.data <- read.delim(paste(datadir,"MMskin_DESeq_result/MMskin_NOskin2_DESeq/MMskin_stat.xls",sep=""),as.is=TRUE,sep="\t")

sum(duplicated(mouse.skin.data$id))

sum(mouse.skin.data$id %in% nmr.skin.data$id)
#NMR.symbol.txt
nmr.symbols <- read.delim(paste(datadir,"NMR.symbol.txt",sep=""),as.is=TRUE,sep="\t")



sym.table <- table(nmr.symbols$mouse)

nmr.symbols.short <- subset(nmr.symbols,mouse!=names(sym.table[which(sym.table>1000)]))

with(nmr.symbols,mean.x(mouse!=names(sym.table[which(sym.table>1000)])))

dim(nmr.symbols.short)


nmr.skin.data$symbol <- gsub(".*#","",nmr.skin.data$mouse_cDNA)

mouse.ids <- mouse.skin.data$id

write(mouse.ids,pastes(datadir,"mouse.ids.txt"))

mean(nmr.skin.data$id  %in% gsub("_seq1","",nmr.symbols.short$transcriptID))

mouse.id.table <- read.delim(pastes(datadir,"mouse_id_batchReport.txt"),as.is=TRUE)


mouse.skin.data$symbol <- mouse.id.table$Symbol[match(mouse.skin.data$id,mouse.id.table$Symbol)]

outcomes <- c("log2FoldChange","baseMeanA_control","baseMeanB_treatment")

mouse.brief <- subset(mouse.skin.data,select=c(outcomes,"symbol"))
nmr.brief <- subset(nmr.skin.data,select=c(outcomes,"symbol"))

names(mouse.brief) <- pastes("mus.",names(mouse.brief))
names(nmr.brief) <- pastes("nmr.",names(nmr.brief))

names(mouse.brief)[grep("symbol",names(mouse.brief))] <- "symbol"
names(nmr.brief)[grep("symbol",names(nmr.brief))] <- "symbol"


mouse.brief.agg <- aggregate(subset(mouse.brief,select=pastes("mus.",outcomes)),by=subset(mouse.brief,select="symbol"),median,na.rm=TRUE)

nmr.brief.agg <- aggregate(subset(nmr.brief,select=pastes("nmr.",outcomes)),by=subset(nmr.brief,select="symbol"),median,na.rm=TRUE)


for(oc.iter in outcomes){

nmr.brief.agg[[pastes("nmr.",oc.iter)]] <- cleaner(nmr.brief.agg[[pastes("nmr.",oc.iter)]])
mouse.brief.agg[[pastes("mus.",oc.iter)]] <- cleaner(mouse.brief.agg[[pastes("mus.",oc.iter)]])

}


sum(mouse.brief.agg$symbol %in% nmr.brief.agg$symbol)

merge.data <- merge(nmr.brief.agg,mouse.brief.agg,by="symbol")

dim(merge.data)

pdf(pastes(resultsdir,"summary_mouse_nmr.pdf"),h=8,w=11)

for(oc.iter in outcomes[c(2,3,1)]){

oc.mouse <- pastes("mus.",oc.iter)
oc.nmr <- pastes("nmr.",oc.iter)

if(oc.iter != outcomes[1]){
	logger <- log10
	}else{
		logger <-  identity
		}

plot(logger(merge.data[[oc.mouse]]),logger(merge.data[[oc.nmr]]),main=paste("NMR vs. Mouse:",oc.iter),xlab="Mouse",ylab="NMR",logy=logger,logx=logger)
cors <- cor.test(cleaner(logger(merge.data[[oc.mouse]])),cleaner(logger(merge.data[[oc.nmr]])))

legend("bottomright",legend=paste(c("Cor =","p ="),c(round(cors$estimate,2),p.rounder(cors$p.value))),inset=0.05,bg="white")
}


dev.off()





pdf(pastes(resultsdir,"summary_clipped_mouse_nmr.pdf"),h=8,w=11)


exp.threshold <- 3

merge.data2 <- subset(merge.data,(nmr.baseMeanA_control>exp.threshold)&(mus.baseMeanA_control>exp.threshold))

for(oc.iter in outcomes[c(2,3,1)]){

oc.mouse <- pastes("mus.",oc.iter)
oc.nmr <- pastes("nmr.",oc.iter)

if(oc.iter != outcomes[1]){
	logger <- log10
	}else{
		logger <- identity		}

plot(logger(merge.data2[[oc.mouse]]),logger(merge.data2[[oc.nmr]]),main=paste("NMR vs. Mouse:",oc.iter),xlab="Mouse",ylab="NMR",logy=logger,logx=logger)
cors <- cor.test(cleaner(logger(merge.data2[[oc.mouse]])),cleaner(logger(merge.data2[[oc.nmr]])))

legend("bottomright",legend=paste(c("Cor =","p ="),c(round(cors$estimate,2),p.rounder(cors$p.value))),inset=0.05,bg="white")
}


dev.off()


normals.controls <- pastes("Sample_",c("MM_skin_C1.norm.log", "MM_skin_C2.norm.log", "MM_Skin1.norm.log", "MM_skin_3.norm.log"))

mouse.skin.data.agg <- aggregate(subset(mouse.skin.data,select=normals.controls),by=subset(mouse.skin.data,select=id),mean,na.rm=TRUE)

names(mouse.skin.data.agg) <- ifelse(names(mouse.skin.data.agg)=="MM_Skin1.norm.log","MM_skin_1.norm.log",names(mouse.skin.data.agg))

nmr.normals.controls <- c("NMR_skin_C1.norm.log" ,"NMR_skin_C2.norm.log" ,"NMR_skin_1.norm.log" ,"NMR_skin_2.norm.log", "NMR_skin_3.norm.log")

nmr.skin.data.agg <- aggregate(subset(nmr.skin.data,select=nmr.normals.controls),by=subset(nmr.skin.data,select=symbol),function(x){mean(cleaner(x),na.rm=TRUE)})

mouse.skin.data.agg$symbol <- mouse.skin.data.agg$id

mouse.skin.data.agg$Species <- "mus"
nmr.skin.data.agg$Species <- "nmr"

gene.pick <- nmr.skin.data.agg$symbol[10]

mouse.skin.data.agg$SYMBOL <- toupper(mouse.skin.data.agg$symbol)
nmr.skin.data.agg$SYMBOL <- toupper(nmr.skin.data.agg$symbol)

all.genes <- unique(nmr.skin.data.agg$SYMBOL)

pdf(pastes(resultsdir,"important_genes_skin.pdf"),h=8,w=11)


for(gene.pick in unique(important.genes$SYMBOL)){
	
	
	plot.mus.nmr.by.gene(gene.pick,mus.d=mouse.skin.data.agg,nmr.d=nmr.skin.data.agg,plot.label="skin")
	
	}


dev.off()

n <- length(all.genes)

#n <- 10

int.p.vec <- matrix(NA,n,1)
rownames(int.p.vec) <- all.genes[1:n]
step <- 1
for(gene.pick in all.genes[1:n]){
	
	
	int.p.vec[gene.pick,1] <- interactions.mus.nmr.by.gene(gene.pick,mus.d=mouse.skin.data.agg,nmr.d=nmr.skin.data.agg,plot.label="skin")
	
	print.skip(step,itermax=length(int.p.vec),percent=5)
	
	step <- step+ 1
	
	}

int.p.vec <- as.numeric(int.p.vec)

most.sig <- all.genes[order(int.p.vec)]


df.out <- data.frame(all.genes,int.p.vec)

df.out <- df.out[order(df.out$int.p.vec),]

write.csv(df.out ,pastes(resultsdir,"skin_most_significant_difference_in_change.csv"),row.names=FALSE)



pdf(pastes(resultsdir,"interaction_p_genes_skin_most_sig.pdf"),h=8,w=11)


for(gene.pick in most.sig[1:150] ){
	
	
	plot.mus.nmr.by.gene(gene.pick,mus.d=mouse.skin.data.agg,nmr.d=nmr.skin.data.agg,plot.label="skin")
	
	}


dev.off()

write.csv(mouse.skin.data.agg,pastes(resultsdir,"mouse_skin.csv"))
write.csv(nmr.skin.data.agg,pastes(resultsdir,"nmr_skin.csv"))



