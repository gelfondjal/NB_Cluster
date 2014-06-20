##############################################################################
# University of Texas Health Science Center at San Antonio
# Department of Epidemiology and Biostatistics                        
##############################################################################
# Filename:  day_analysis_v0.3.R                                
# Author:    Jonathan Gelfond                                             
# Project Name:  Shireman Revision         
# Input:        E2F3 Results.csv, miR351 Results.csv                     
# Output: 
#
# Modification History:  
# v 0.1 Creation    
#
#  0.3 Plots the CI for the ANOVA effects
#  includes keeping the low baseline values for "snoRNA55"
                          
##############################################################################





library(plyr)
library(gdata)
library(gplots)


set.seed(2011)


sourceFile <- "day_analysis_v0.3.R"


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

indata$ct.cut <- with(indata,ifelse(Ct>ct.threshold,NA,as.numeric(Ct)))
indata$ct.cut <- with(indata,ifelse((day.num==0)&(Detector.Name %in% c("microRNA-351","snoRNA55")),as.numeric(Ct),as.numeric(ct.cut)))

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

pdf(pastes(resultsdir,"boxplots_RNA_delta_delta_2.pdf"),height=8,width=11)


for(pars in 1:1){
	


par(mfrow=c(1,pars))


for(gene.iter in 1:length(numerators)){

indata <- subset(all.summary,Gene==numerators[gene.iter])

indata$rawratio <- 2^-indata$logratio

indata$logratio <- -indata$logratio

indata$day.factor <- with(indata,as.factor(day.num))

pretty.name <- numerators[gene.iter]

anova.fitter <- lm(logratio~day.factor,data=indata)

estimates <- coef(anova.fitter)[-1]

CIs <- confint(anova.fitter)[-1,]

pvals <- summary(anova.fitter)$coef[,4][-1]

plot.data <- data.frame(pvals,CIs,est=as.matrix(estimates),day=as.factor(as.numeric(gsub("day\\.factor","",rownames(CIs)))))

#with(plot.data,plotCI(as.numeric(day),est,ui=plot.data$X97,li=plot.data$X2.5,xlab="Day",xaxt="n",ylab="Log Fold Change, Delta Delta CT",main=paste(pretty.name,"vs. baseline")))
#with(plot.data,axis(1,labels=levels(day),at=as.numeric(day)))
#abline(h=0)

ylimmer <- c( min(c(1,2^min(plot.data$X2.5))  ),max(2^plot.data$X97))
xlimmer <- with(plot.data,range(as.numeric(day)))+c(-0.25,.25)

with(plot.data,plotCI(as.numeric(day),2^est,ui=2^plot.data$X97,li=2^plot.data$X2.5,xlab="Day",xaxt="n",main=paste(pretty.name,"vs. baseline"),log="y",ylim=ylimmer,xlim=xlimmer,ylab=expression(paste(2^{paste("-",Delta,Delta,"CT")})),gap=0))
with(plot.data,axis(1,labels=levels(day),at=as.numeric(day)))
abline(h=1,lty=3)

with(plot.data,text(as.numeric(day),rep(1,nrow(plot.data)),paste(c("p=",rep("",nrow(plot.data)-1)),p.rounder(pvals,digits=2)),pos=3))


legend("topright",inset=0.01,legend=c("Fold Change (95% CI)"),pch=c(1),lty=1,lwd=1,col=c("black"),bg="white")

} # loop over gene.inter

} #loop over pars
dev.off()






plot(c(1,1),c(0,1),ylab=expression("A**2 Delta"))

labNames <- c('xLab','yLabl')
plot(c(1:10),xlab=expression(paste(labName[1], x^Delta,Delta)),ylab=expression(paste(labName[2], y^2)))



labNames <- c('xLab','yLabl')
plot(c(1:10),xlab=expression(paste(labName[1], x^Delta)),ylab=expression(paste(labName[2], y^{(paste(delta,delta))},delta,delta)))

























