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
library(asbio)

set.seed(2012)


sourceFile <- "analyze_liver_v0.1.R"
sourceFileClip <- gsub("\\.R","/",sourceFile)

base1 <- "/Users/jonathangelfond/Documents/Projects/Buffenstein/NMR_Genomics/"
#base2 <- "/Users/jonathangelfond/Dropbox/Temp_Work/Weiner/"

basedir <- base1

analysisdir <- paste(basedir,"Programs/",sep="") 

datadir <- paste(basedir,"Data/",sep="")

trinity.data.file <- "/Users/jonathangelfond/Documents/Projects/Buffenstein/NMR_Genomics/Data/BLAST_Trinity/Triniy_mouse_cDNA_blastn.txt"

resultsdir <- paste(basedir,"Documents/Results/",sourceFileClip,sep="")
resultsdir.r2sas <- paste(resultsdir,"R2SAS/",sep="")
resultsdir.sas <- paste("//psf/Home/Dropbox/OSCE_Project2/Jons_Template/Documents/Results/",sourceFileClip,"R2SAS/",sep="")


source("/Users/jonathangelfond/Documents/Projects/sudheer/Analysis/table_functions.R")

loudenPath <- "/Users/jonathangelfond/Documents/Projects/Louden/"
sapply(paste(loudenPath,list.files(loudenPath),sep=""),source)
source(pastes(analysisdir,"support_functions.R"))


dir.create(resultsdir)
dir.create(resultsdir.r2sas)

# important genes

important.genes <- read.xls(pastes(datadir,"RNASeq Genes.xlsx"))
important.genes$SYMBOL <- toupper(important.genes$Official.Symbol)
# Read in Sequence Data

nmr.liver.data <- read.delim(paste(datadir,"NMRliver_DESeq_result/NMRliver_stat_symbol.xls",sep=""),as.is=TRUE,sep="\t")

sum(duplicated(nmr.liver.data$id))

mouse.liver.data <- read.delim(paste(datadir,"MMliver_DESeq_result/MMliver_stat.xls",sep=""),as.is=TRUE,sep="\t")

sum(duplicated(mouse.liver.data$id))

sum(mouse.liver.data$id %in% nmr.liver.data$id)
#NMR.symbol.txt
nmr.symbols <- read.delim(paste(datadir,"NMR.symbol.txt",sep=""),as.is=TRUE,sep="\t")


sym.table <- table(nmr.symbols$mouse)

nmr.symbols.short <- subset(nmr.symbols,mouse!=names(sym.table[which(sym.table>1000)]))

with(nmr.symbols,mean.x(mouse!=names(sym.table[which(sym.table>1000)])))

dim(nmr.symbols.short)


nmr.liver.data$symbol <- gsub(".*#","",nmr.liver.data$mouse_cDNA)

mouse.ids <- mouse.liver.data$id

write(mouse.ids,pastes(datadir,"mouse.ids_liver.txt"))

mean(nmr.liver.data$id  %in% gsub("_seq1","",nmr.symbols.short$transcriptID))

mouse.id.table <- read.delim(pastes(datadir,"mouse_id_batchReport.txt"),as.is=TRUE)


mouse.liver.data$symbol <- mouse.id.table$Symbol[match(mouse.liver.data$id,mouse.id.table$Symbol)]

outcomes <- c("log2FoldChange","baseMeanA_control","baseMeanB_treatment")

mouse.brief <- subset(mouse.liver.data,select=c(outcomes,"symbol"))
nmr.brief <- subset(nmr.liver.data,select=c(outcomes,"symbol"))

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


nmr.liver.id.normed <- subset(nmr.liver.data,select=c("id",grep("norm\\.log",names(nmr.liver.data),value=TRUE),"symbol"))


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



# Process BLAST Data

blast.data <- read.delim(trinity.data.file)

blast.names <- gsub("^ *","",c(" NMR_ID","Length of NMR transcript","Mouse ID","Length of mouse transcript",
					" E-value","Identities","Strand","Query start position","Query end position",
					"Subject start positio","Subject end position"))
					
blast.header <- make.names(blast.names)
names(blast.data) <- blast.header
blast.data$nmr.id <- gsub("_seq","",gsub(" .*","",blast.data$NMR_ID))

x <- as.character(blast.data$Identities[1])

Divider <- function(x){ 
	# Takes "A/B" returns A/B
	y <- as.numeric(strsplit(as.character(x),split="/")[[1]])
	return(y[1]/y[2])
}

blast.data$Pct.ID <- apply(matrix(as.character(blast.data$Identities)),1,Divider)



# Export Transcripts

write(unique(as.character(blast.data$ens.num)),pastes(resultsdir,"transcript_ids.csv"))

blast.data$ens.num <- gsub(" .*","",gsub("^ ","",as.character(blast.data$Mouse.ID)))# Import mapping to symbol

ens.num.to.symbol <- read.delim(pastes(datadir,"ensemble_id_match.txt"))
# Remove Duplicates 
ens.num.to.symbol <- subset(ens.num.to.symbol,!duplicated(Input))
blast.data.2 <- merge(blast.data,ens.num.to.symbol,by.x="ens.num",by.y="Input")

mean(as.character(blast.data.2$ens.num) %in% ens.num.to.symbol$Input)

#blast.data.2<- subset(ens.num.to.symbol,!duplicated(Input))

blast.data.2$id <- gsub(" .*","",blast.data.2$NMR_ID)

mean(blast.data.2$id %in% nmr.liver.id.normed$id)

blast.data.2$n <- 1

match.quality <- aggregate(subset(blast.data.2,select="Pct.ID"),by=subset(blast.data.2,select="Symbol"),mean,na.rm=TRUE)

match.numbers <- aggregate(subset(blast.data.2,select="n"),by=subset(blast.data.2,select="Symbol"),sum,na.rm=TRUE)

match.stats <- merge(match.numbers,match.quality,by="Symbol")

pdf(pastes(resultsdir,"mapping_plots.pdf"),h=8,w=11)

mean(nmr.liver.data$id  %in% blast.data$nmr.id )
mapping.reps <- table(blast.data$nmr.id)
plot(1:length(mapping.reps),sort(mapping.reps,decreasing=TRUE),main="Maps from NMR to >1 Mouse Sequence",xlab="Number of NMR Seqs",ylab="Number of Mapping",type="l")
abline(h=1)
legend("topright",legend=paste(sum(1-duplicated(blast.data$nmr.id)),"Total NMR Seqs"),inset=0.02)
hist(mapping.reps[mapping.reps<100],breaks=100)



legend("bottomright",legend=paste(sum(1-duplicated(blast.data$ens.num)),"Total Mouse Genes"),inset=0.1)

hist(blast.data$Pct.ID,main="Histogram of Match Quality",xlab="% Match")

x <- with(blast.data, abs(Query.start.position-Query.end.position))
y <- blast.data$Pct.ID
plot(cbind(x,y)[as.logical(rbinom(length(y),1,0.1)),],ylab="% ID",xlab="Seq Length",pch=".",main="Match Quality vs. Sequence Length")
abline(coef=lm(y~x)$coef,lty=2,lwd=2,col=2)

hist(match.quality$Pct.ID,main="% ID per Gene")
hist(match.numbers$n[match.numbers$n<400],main="Number NMR seq per mouse Gene")

n.threshold <- 100

with(subset(match.stats,n<n.threshold),plot(n,Pct.ID,main="Match Quality vs. Number of Matches",xlab="Number of NMR Seq Matches per Mouse Gene"))

loessFit <- loess(Pct.ID~n,data=subset(match.stats,n<n.threshold))
points(1:n.threshold,predict(loessFit,newdata=data.frame(n=1:n.threshold)),pch=19,col=2)


dev.off()


names(match.stats) <- ifelse(names(match.stats)=="Symbol","symbol",names(match.stats))


merge.data.matches <- merge(merge.data,match.stats,by="symbol")


pdf(pastes(resultsdir,"summary_mouse_nmr_match_quality.pdf"),h=8,w=11)

for(oc.iter in outcomes[c(2,3,1)]){

oc.mouse <- pastes("mus.",oc.iter)
oc.nmr <- pastes("nmr.",oc.iter)

if(oc.iter != outcomes[1]){
	logger <- log10
	}else{
		logger <-  identity
		}

plot(logger(merge.data.matches[[oc.mouse]]),logger(merge.data.matches[[oc.nmr]]),main=paste("NMR vs. Mouse:",oc.iter),xlab="Mouse",ylab="NMR",logy=logger,logx=logger)


plot(merge.data.matches$Pct.ID,logger(merge.data.matches[[oc.mouse]])-logger(merge.data.matches[[oc.nmr]]),main=paste("NMR vs. Mouse:",oc.iter,"by % ID"),ylab="Mouse-NMR",xlab="% ID with Mouse",logy=logger,logx=logger)

x <- cleaner(merge.data.matches$Pct.ID)
y <- cleaner(logger(merge.data.matches[[oc.mouse]])-logger(merge.data.matches[[oc.nmr]]))

x1 <- logger(merge.data.matches[[oc.mouse]])
x2 <- logger(merge.data.matches[[oc.nmr]])
x3 <- cleaner(merge.data.matches$Pct.ID) 

loessFit <- loess(y~x,data=data.frame(x=x,y=y))
points(sort(x),predict(loessFit,newdata=data.frame(x=sort(x))),pch=19,col=2)


data.mat <- na.exclude(cbind(cleaner(x1),cleaner(x2),x3))

colnames(data.mat) <- c("Mouse","NMR","Pct.ID")

data.mat.small <- data.mat[rbinom(nrow(data.mat),1,p=0.9)==1,]

loess2d <- loess(data.mat.small[,"Pct.ID"]~data.mat.small[,c("Mouse","NMR")])
preds <- predict(loess2d) 

n.colors <- 256

pred.fit <- round(n.colors*(preds-min(preds,na.rm=TRUE))/abs(diff(range(preds,na.rm=TRUE))),0)
pred.fit <- ifelse(pred.fit==0,1,pred.fit)

#plot(pred.fit,preds)

plot(x1,x2,pch=".",xlab="Mouse",ylab="NMR",main=paste("NMR vs. Mouse:",oc.iter))

colors <- bluered(n.colors)

points(data.mat.small[,1:2],col=colors[pred.fit],pch=19)

legend("topleft",legend=c("Max Pct ID","Min Pct ID"),pch=19,col=c(colors[n.colors],colors[1]))




plot(log10(merge.data.matches$n),logger(merge.data.matches[[oc.mouse]])-logger(merge.data.matches[[oc.nmr]]),main=paste("NMR vs. Mouse:",oc.iter,"by # of matches"),ylab="Mouse-NMR",xlab="log10(# NMR seq) matches with Mouse",logy=logger,logx=logger)

x <- log10(cleaner(merge.data.matches$n))
y <- cleaner(logger(merge.data.matches[[oc.mouse]])-logger(merge.data.matches[[oc.nmr]]))


abline(h=0,lty=2,col=3,lwd=4)

loessFit <- loess(y~x,data=data.frame(x=x,y=y))
points(sort(x),predict(loessFit,newdata=data.frame(x=sort(x))),pch=19,col=2)




}


dev.off()


normals.controls <- c("MM_liver_C1.norm.log", "MM_liver_C2.norm.log", "MM_liver_1.norm.log", "MM_liver_2.norm.log", "MM_liver_3.norm.log")

mouse.liver.data.agg <- aggregate(subset(mouse.liver.data,select=normals.controls),by=subset(mouse.liver.data,select=id),mean,na.rm=TRUE)

nmr.normals.controls <- c("NMR_liver_C1.norm.log" ,"NMR_liver_C2.norm.log" ,"NMR_liver_1.norm.log" ,"NMR_liver_2.norm.log", "NMR_liver_3.norm.log")

nmr.liver.data.agg <- aggregate(subset(nmr.liver.data,select=nmr.normals.controls),by=subset(nmr.liver.data,select=symbol),function(x){mean(cleaner(x),na.rm=TRUE)})

mouse.liver.data.agg$symbol <- mouse.liver.data.agg$id

mouse.liver.data.agg$Species <- "mus"
nmr.liver.data.agg$Species <- "nmr"

gene.pick <- nmr.liver.data.agg$symbol[10]

mouse.liver.data.agg$SYMBOL <- toupper(mouse.liver.data.agg$symbol)
nmr.liver.data.agg$SYMBOL <- toupper(nmr.liver.data.agg$symbol)

pdf(pastes(resultsdir,"important_genes_liver.pdf"),h=8,w=11)


for(gene.pick in unique(important.genes$SYMBOL)){
	
	
	plot.mus.nmr.by.gene(gene.pick,mus.d=mouse.liver.data.agg,nmr.d=nmr.liver.data.agg,plot.label="Liver")
	
	}


dev.off()



write.csv(mouse.liver.data.agg,pastes(resultsdir,"mouse_liver.csv"))
write.csv(nmr.liver.data.agg,pastes(resultsdir,"nmr_liver.csv"))

mouse.all <- merge(mouse.liver.data.agg,mouse.skin.data.agg,by="SYMBOL")
nmr.all <- merge(nmr.liver.data.agg,nmr.skin.data.agg,by="SYMBOL")

mouse.data <- as.matrix(mouse.all[,grep("log",names(mouse.all),value=TRUE)])
rownames(mouse.data) <- mouse.all$SYMBOL

nmr.data <- as.matrix(nmr.all[,grep("log",names(nmr.all),value=TRUE)])
rownames(nmr.data) <- nmr.all$SYMBOL








