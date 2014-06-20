##############################################################################
# University of Texas Health Science Center at San Antonio
# Department of Epidemiology and Biostatistics                        
##############################################################################
# Filename:  day_analysis_v0.1.R                                
# Author:    Jonathan Gelfond                                             
# Project Name:  Shireman Revision         
# Input:        E2F3 Results.csv, miR351 Results.csv                     
# Output: 
#
# Modification History:  
# v 0.1 Creation                              
##############################################################################





library(plyr)
library(gdata)
library(gplots)


set.seed(2011)


sourceFile <- "day_analysis_v0.1.R"


sourceFileClip <- gsub("\\.R","/",sourceFile)

base1 <- "/Users/jonathangelfond/Documents/Projects/Shireman_miRNA/Menton/SOP_Template/"
#base2 <- "/Users/jonathangelfond/Dropbox/Temp_Work/Weiner/"

basedir <- base1

analysisdir <- paste(basedir,"Programs/",sep="") 


datadir <- paste(basedir,"Data/",sep="") 


resultsdir <- paste(basedir,"Documents/Results/",sourceFileClip,sep="")
resultsdir.r2sas <- paste(resultsdir,"R2SAS/",sep="")
resultsdir.sas <- paste("//psf/Home/Documents/Projects/Shireman_miRNA/Menton/SOP_Template/Documents/Results/",sourceFileClip,"R2SAS/",sep="")

source(paste(analysisdir,"support_functions.R",sep=""))
source("/Users/jonathangelfond/Documents/Projects/Weiner/Analysis/support_functions.R")
source("/Users/jonathangelfond/Documents/Projects/sudheer/Analysis/table_functions.R")
source("/Users/jonathangelfond/Documents/Projects/OB_Residents/Brumm/Programs/support_functions.R")
chris.path <- "/Users/jonathangelfond/Dropbox/Work_Shared/R Scripts/Chris/"


chris.path2 <- "/Users/jonathangelfond/Documents/Projects/Louden/"

temp <- lapply(list.files(chris.path,full.names=1,pattern="r$"),source)
temp <- lapply(list.files(chris.path2,full.names=1,pattern="r$"),source)


dir.create(resultsdir)
dir.create(resultsdir.r2sas)

files <- c("E2F3 Results.csv","miR351 Results.csv")

numerators <- c("E2f3","microRNA-351")

denominators <- c("GAPDH","snoRNA55")

skips <- c(0,1)

all.summary <- NULL

for(gene.iter in 1:2){


indata <- read.delim(pastes(datadir,files[gene.iter]),sep=",",as.is=TRUE,header=TRUE,skip=skips[gene.iter])
indata <- indata[,1:10]
indata <- subset(indata,Sample.Name!="")

indata$sample.num <- with(indata,gsub(".*#","",Sample.Name))

indata$day <- with(indata,gsub("d.*","",gsub(" .*","",Sample.Name)))

indata$day.num <- with(indata,ifelse(day=="Base",0,as.numeric(day)))

table(indata$day.num)

ct.threshold <- 30

indata$ct.cut <- with(indata,ifelse(Ct>30,NA,as.numeric(Ct)))


with(indata,Ct[which(is.na(as.numeric(Ct)))])

control.data <- aggregate(list(median=indata$ct.cut),subset(indata,select=c("sample.num","day.num","Detector.Name")),median,na.rm=TRUE)

numerator   <- numerators[gene.iter]
denominator <- denominators[gene.iter]

lograts <- ddply(control.data,c("sample.num","day.num"),function(x){
	num <- x$median[x$Detector.Name==numerator]
	denom <- x$median[x$Detector.Name==denominator]
	ratio <- num-denom
	names(ratio) <- "logratio"	
	return(ratio)
	})
	
lograts$Gene <- numerator 

all.summary <- rbind(all.summary,lograts)

}




	
	
	




pdf(pastes(resultsdir,"boxplots_RNA.pdf"),height=8,width=11)

for(gene.iter in 1:length(numerators)){

indata <- subset(all.summary,Gene==numerators[gene.iter])

pretty.name <- numerators[gene.iter]

boxstats <- boxplot(logratio~as.factor(day.num),data=indata,ylog=TRUE,ylab="Delta-CT",xlab="day",main=paste(pretty.name,"Control",denominators[gene.iter]))
#temp <- plotmeans(logratio~as.factor(day.num),data=indata,add=TRUE,ylab="",xlab="")
mtext(at=1:length(boxstats$n),line=-1,pastes("N=",boxstats$n),side=1)


stripchart(logratio~as.factor(day.num),data=indata,add=TRUE,vertical=TRUE,pch=19,col=2)

logmeans <- aggregate(list(logmean = indata$logratio),by=subset(indata,select="day.num"),mean,na.rm=TRUE)
logCIs <- aggregate(list(tCI = indata$logratio),by=subset(indata,select="day.num"),tCI,conf.level=0.9)
IQRs <- aggregate(list(Percentiles = indata$logratio),by=subset(indata,select="day.num"),quantile,probs=c(0.25,0.5,0.75),na.rm=TRUE)


logmeans$geo.mean <- exp(logmeans$logmean)
logCIs$geo.mean.lwr90CI <- exp(logCIs$tCI[,1])
logCIs$geo.mean.upr90CI <- exp(logCIs$tCI[,2])


all.stats <- merge(logmeans,logCIs,by=c("day.num"))
all.stats <- merge(all.stats,IQRs,by=c("day.num"))
all.stats$x.pos <- as.numeric(as.factor(all.stats$"day.num"))

plotCI(all.stats$x.pos,uiw=0.5,all.stats$logmean,ui=all.stats$tCI[,2],li=all.stats$tCI[,1],add=TRUE,lwd=3,pch=4,col="grey")

legend("topleft",inset=0.01,legend=c("Mean (90% CI)"),pch=c(14),lty=1,lwd=2,col=c("grey"),bg="white")

}


dev.off()












pdf(pastes(resultsdir,"boxplots_RNA_raw.pdf"),height=8,width=11)

for(gene.iter in 1:length(numerators)){

indata <- subset(all.summary,Gene==numerators[gene.iter])

indata$rawratio <- 2^-indata$logratio

pretty.name <- numerators[gene.iter]

plot(indata$rawratio~as.factor(day.num),data=indata,log="y",ylab="2^-(delta CT)",xlab="Day",main=paste(pretty.name,"Control",denominators[gene.iter]))

boxstats <- boxplot(logratio~as.factor(day.num),data=indata,ylog=TRUE,ylab="Delta-CT",xlab="day",plot=FALSE)
#temp <- plotmeans(logratio~as.factor(day.num),data=indata,add=TRUE,ylab="",xlab="")
mtext(at=1:length(boxstats$n),line=-1,pastes("N=",boxstats$n),side=1)


stripchart(rawratio~as.factor(day.num),data=indata,add=TRUE,vertical=TRUE,pch=19,col=1)

logmeans <- aggregate(list(logmean = indata$logratio),by=subset(indata,select="day.num"),mean,na.rm=TRUE)
logCIs <- aggregate(list(tCI = indata$logratio),by=subset(indata,select="day.num"),tCI,conf.level=0.9)
IQRs <- aggregate(list(Percentiles = indata$logratio),by=subset(indata,select="day.num"),quantile,probs=c(0.25,0.5,0.75),na.rm=TRUE)


logmeans$geo.mean <- 2^-(logmeans$logmean)
logCIs$geo.mean.lwr90CI <- 2^-(logCIs$tCI[,2])
logCIs$geo.mean.upr90CI <- 2^-(logCIs$tCI[,1])


all.stats <- merge(logmeans,logCIs,by=c("day.num"))
all.stats <- merge(all.stats,IQRs,by=c("day.num"))
all.stats$x.pos <- as.numeric(as.factor(all.stats$"day.num"))

with(all.stats,matlines(rbind(x.pos,x.pos),rbind(geo.mean.upr90CI,geo.mean.lwr90CI),lwd=3,pch=4,col="grey",lty=1))
with(all.stats,points(x.pos,geo.mean,pch=4,lwd=2,col="grey"))

legend("topleft",inset=0.01,legend=c("Geo. Mean (90% CI)"),pch=c(14),lty=1,lwd=2,col=c("grey"),bg="white")

}


dev.off()
































