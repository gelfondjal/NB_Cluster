##############################################################################
# University of Texas Health Science Center at San Antonio
# Department of Epidemiology and Biostatistics                        
##############################################################################
# Filename: analyze_v0.3.R
# Author:              Jonathan Gelfond                                   
# Project Name:          KLugman OSCE
# Input:                  Codebook Scoresheet.csv    & Code Book Data Sheet.csv
# Output: 
#
# Modification History:
# v 0.1 Creation   
# 0.2 add intervention status & case study scores     
# 0.3 add modification for poster output                      
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


sourceFile <- "analyze_v0.5.R"
sourceFileClip <- gsub("\\.R","/",sourceFile)

base1 <- "/Users/jonathangelfond/Dropbox/OSCE_Project2/Jons_Template/"
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


# Read in Video Analysis Data

osceData <- read.csv(paste(datadir,"Code Book Data Sheet.csv",sep=""),as.is=TRUE)
varData <- read.csv(paste(datadir,"Codebook Scoresheet_Direction.csv",sep=""),as.is=TRUE)

# Read in Standardized Patient Data


patientData <- read.csv(paste(datadir,"Patient_eval_questions.csv",sep=""),as.is=TRUE)

names(patientData) <- gsub(" $","",names(patientData))


patient.scan <- scan(paste(datadir,"Patient_eval_questions.csv",sep=""),what=character(),sep=",")


osceData$Subject <- gsub(", ",",",osceData$Subject..Last..First.)

mean(osceData$Subject %in% patientData$STUDENT.NAME)

missing <- osceData$Subject[!(osceData$Subject %in% patientData$STUDENT.NAME)]

patientData$STUDENT.NAME[pmatch(missing,patientData$STUDENT.NAME)]
osceData$Subject <- ifelse(osceData$Subject!="Zoffato,Teresa",osceData$Subject,"Zoffuto,Teresa")
osceData$Subject <- ifelse(osceData$Subject!="Hernard,Kristin",osceData$Subject,"Hemard,Kristin")

osceData$Subject[!(osceData$Subject %in% patientData$STUDENT.NAME)]


patientData$Subject <- patientData$STUDENT.NAME

merged.osce <- merge(osceData,patientData,by="Subject")


patient.evals <- names(patientData)[grep("^IMC",names(patientData))]





varData2 <- subset(varData,!is.na(QNum))
varData2$QNum <- pastes("Q",varData2$QNum)

names(osceData)

# Clean UP the data:

length(unique(grep("^Q",names(osceData))))
qnames <- names(osceData)[grep("^Q",names(osceData))]

for(namer in qnames){
	
	osceData[[namer]] <- gsub(" ","",	osceData[[namer]])
}


# Q41 is the time variable

osceData2 <- osceData

seconds <- as.numeric(gsub(".*:","",osceData2$Q41))
minutes <- as.numeric(gsub(":.*","",osceData2$Q41))

osceData2$Q41 <- minutes + (seconds/60)

# Rename variables meaningfully

for(name.row in 1:nrow(varData2)) {
	newname <- varData2$ShortName[name.row]
	osceData2[[newname]] <- osceData2[[varData2$QNum[name.row]]]
}


# Read in Case Summary scores data


case.summary.data <- read.csv(pastes(datadir,"IMC AD Case-Category Scores.csv"),as.is=TRUE)

case.summary.data <- case.summary.data[,1:which(names(case.summary.data)=="X.1")] # 6 is the number of meaningful columns

case.summary.data$Subject <- case.summary.data$Student.Name

# Read in Case Question scores data


case.questions.data <- read.csv(pastes(datadir,"IMC AD Case Only.csv"),as.is=TRUE,skip=1)
case.questions.data$Subject <- gsub(",",", ",case.questions.data$Subject)
names(case.questions.data) <- gsub("Interpesonal","Interpersonal",names(case.questions.data))

important.questions <- names(case.questions.data)[grep("(Interpersonal)|(History)",names(case.questions.data))]

for(var in important.questions){case.questions.data[[var]] <- ifelse(case.questions.data[[var]]==-1,NA,case.questions.data[[var]])}


# Read in Patient assignment data

intervention.data <- read.csv(pastes(datadir,"Med Students.csv"),as.is=TRUE)

intervention.data$Subject <- intervention.data$Student

# Read in self evaluation data

self.eval.data <- read.csv(pastes(datadir,"student_self_eval_questions.csv"),as.is=TRUE)

names(self.eval.data) <- gsub("STUDENT.NAME","Subject",names(self.eval.data))

self.questions <- names(table(self.eval.data$Question))

question.table <- matrix(c(c("Self.Review","Self.Comfortable","Self.Legal.Authority"),self.questions),length(self.questions),2)

self.eval.data$Short.Answer <- question.table[match(self.eval.data$Question,question.table[,2]),1]

with(self.eval.data,table(Short.Answer,Question))

with(self.eval.data,table(RESPONSES))

self.eval.wide <- reshape(subset(self.eval.data,select=c("Subject","RESPONSES","Short.Answer")),idvar="Subject",v.names="RESPONSES",direction="wide",timevar="Short.Answer")

setdiff(self.eval.wide$Subject ,case.summary.data$Subject)

setdiff(case.summary.data$Subject,self.eval.wide$Subject )

self.eval.wide$Subject <- ifelse(self.eval.wide$Subject!= "Garcia-Reyes, Lawren",self.eval.wide$Subject,"Garcia-Reyes, Lawrence")
self.eval.wide$Subject <- ifelse(self.eval.wide$Subject!="De La Garza, Elizabe",self.eval.wide$Subject,"De La Garza, Elizabeth")
self.eval.wide$Subject <- ifelse(self.eval.wide$Subject!="Osadebe, Uchechukwuk",self.eval.wide$Subject,"Osadebe, Uchechukwuka")


setdiff(self.eval.wide$Subject ,case.summary.data$Subject)


match.likert <- c("Strongly Agree","Agree","Neutral/Don't Know","Disagree","Strongly Disagree")

self.eval.wide$Self.Comfortable.1 <- self.eval.wide$RESPONSES.Self.Comfortable

self.eval.wide$Self.Comfortable.2 <- match(self.eval.wide$Self.Comfortable.1,match.likert)

self.eval.wide$Self.Review <- self.eval.wide$RESPONSES.Self.Review

with(self.eval.wide,table(Self.Comfortable.1,Self.Comfortable.2))


# Get only the relevant students

table(intervention.data$Group)




intervention.data.yn <- subset(intervention.data,Group %in% c("Control","Intervention"))

# Check the names

mean(intervention.data.yn$Subject %in% case.summary.data$Subject)

intervention.data.yn$Subject[!(intervention.data.yn$Subject %in% case.summary.data$Subject)]
intervention.data.yn$Subject <- ifelse(intervention.data.yn$Subject!="Zoffato, Teresa",intervention.data.yn$Subject,"Zoffuto, Teresa")
intervention.data.yn$Subject <- ifelse(intervention.data.yn$Subject!="Hernard, Kristin",intervention.data.yn$Subject,"Hemard, Kristin")
intervention.data.yn$Subject <- ifelse(intervention.data.yn$Subject!="Agapito,Adrian",intervention.data.yn$Subject,"Agapito, Adrian")

mean(intervention.data.yn$Subject %in% case.summary.data$Subject)

# merge

sum(intervention.data.yn$Subject %in% case.summary.data$Subject)

intervention.data.all <- merge(intervention.data.yn,case.summary.data,by="Subject")

mean(intervention.data.yn$Subject %in% case.questions.data$Subject)

setdiff(intervention.data.yn$Subject,case.questions.data$Subject)

intervention.data.all <- merge(intervention.data.all,case.questions.data,by="Subject")




osce.data.3 <- osceData2

osce.data.3$Subject <- gsub(",",", ",osce.data.3$Subject)

mean(intervention.data.all$Subject %in% osce.data.3$Subject)

osce.data.all <- merge(intervention.data.all,osce.data.3,by="Subject")

osce.data.all <- merge(osce.data.all,self.eval.wide,by="Subject")

osce.data.all$livingWill <- ifelse(osce.data.all$livingWill=="C","A",osce.data.all$livingWill)

osce.data.all$handGest <- ifelse(osce.data.all$handGest=="",NA,osce.data.all$handGest)



summary.vars <- c("Case.Score" ,"History" ,"Interpersonal.Skills")


summary(aov(Case.Score ~ Self.Comfortable.1,data=osce.data.all))


self.comfort.case.score <- lm(Case.Score ~ Self.Comfortable.2,data=osce.data.all)
self.review.case.score <- lm(Case.Score ~ Self.Review,data=osce.data.all)

self.comfort.history.score <- lm(History ~ Self.Comfortable.2,data=osce.data.all)
self.review.history.score <- lm(History ~ Self.Review,data=osce.data.all)


self.comfort.interpersonal.score <- lm(`Interpersonal.Skills` ~ Self.Comfortable.2,data=osce.data.all)
self.review.interpersonal.score <- lm(`Interpersonal.Skills` ~ Self.Review,data=osce.data.all)


regressions <- list(self.comfort.case.score,self.review.case.score,self.comfort.history.score,self.review.history.score,self.comfort.interpersonal.score,self.review.interpersonal.score)
names(regressions) <- c("self.comfort.case.score","self.review.case.score","self.comfort.history.score","self.review.history.score","self.comfort.interpersonal.score","self.review.interpersonal.score") 



lm(Case.Score ~ Self.Review,data=osce.data.all)


with(osce.data.all,fisher.test(Self.Review,Group))

with(osce.data.all,fisher.test(Self.Comfortable.1,Group))

self.vars <- c("Self.Comfortable.1","Self.Comfortable.2","Self.Review")

orderedVars <- which(names(osce.data.all) %in% self.vars)

bySelfTable <- summarize.data.frame(data.=osce.data.all,outcome.=as.factor(osce.data.all$Group),columns.=orderedVars)$table



bySelfTable <- makePrettyTests(subset(bySelfTable,select=setdiff(names(bySelfTable),c("row","panel"))))

byself.table <- simple.table.n(resultsdir.r2sas,resultsdir.sas,"selfdatasummary",as.matrix(bySelfTable),
				title=pastes("Table 1. Group effect on Self Evaluation.",ifelse(FALSE,"","cont'd")))

regression.tables <- list()

for(i in 1:length(regressions)){
				
regression.tables[[i]] <- simple.table.n(resultsdir.r2sas,resultsdir.sas,names(regressions)[i],as.matrix(linearRegression.table(regressions[[i]])),
				title=paste("Table 2. Case Scores regressed onto Self Evaluation",names(regressions)[i],ifelse(i>1,"","cont'd")))
}



output.tables.0 <- c(byself.table,unlist(regression.tables))

summaryTablesTable4SAS <- simple.multi.table(resultsdir.r2sas,resultsdir.sas,"self_comparison","self_comparison.rtf",output.tables.0,close.ods=TRUE)




importantVars <- c(varData2$ShortName)

importantVars[!(importantVars %in% names(osce.data.all))]
#summary(subset(allData,select=importantVars))

for(i in importantVars){
if(!is.numeric(osceData2[[i]])){
	print(i)
	print(table(osceData2[i]))
}
}


set.seed(2011)

# Do All the subjects:

orderedVars <- which(names(osce.data.all) %in% importantVars)

byGroupTable <- summarize.data.frame(data.=osce.data.all,outcome.=as.factor(osce.data.all$Group),columns.=orderedVars)$table


byGroupTable.all <- summarize.data.frame(data.=osce.data.all,outcome.=as.factor(osce.data.all$Group),columns.=orderedVars)

pdf(pastes(resultsdir,"pvalues_by_group.pdf"))
hist(byGroupTable.all$pvalues,breaks=20)
dev.off()

alpha01 <- .1
byGroupTable.sig <- subset(byGroupTable.all$table,panel %in% which(byGroupTable.all$pvalues<alpha01))

orderedVars.scores <- which(names(osce.data.all) %in% summary.vars)

byScoreTable <- as.data.frame(summarize.data.frame(data.=osce.data.all,outcome.=as.factor(osce.data.all$Group),columns.=orderedVars.scores)$table)

byScoreTablePretty <- makePrettyTests(subset(byScoreTable,select=setdiff(names(byScoreTable),c("row","panel"))))

byGroupTable.sig.pretty <- makePrettyTests(subset(byGroupTable.sig,select=setdiff(names(byGroupTable.sig),c("row","panel"))))

byGroupTable.sig.pretty$test <- ifelse(byGroupTable.sig.pretty$Pval=="","",byGroupTable.sig.pretty$test)

score.table <- simple.table.n(resultsdir.r2sas,resultsdir.sas,"scoredatasummary",as.matrix(byScoreTablePretty),
				title=pastes("Table 1. Summarized patient data by group.",ifelse(FALSE,"","cont'd")))


#orderedVars2 <- orderedVars[-c(3)]

#bySPTable <- summarize.data.frame(data.=osceData2,outcome.=as.factor(osceData2$SP..),columns.=orderedVars2)$table

byGroupTablePretty <- makePrettyTests(subset(byGroupTable,select=setdiff(names(byGroupTable),c("row","panel"))))

pages <- Split.runs(byGroupTable$panel,21)

group.table <- list()
for(page in 1:max(pages)){
group.table[[page]] <- simple.table.n(resultsdir.r2sas,resultsdir.sas,pastes("alldatasummary",page),as.matrix(byGroupTablePretty[pages==page,]),
				title=pastes("Table 2. Summarized video analysis data by group.",ifelse(page==1,"","cont'd")))
}


byQuestionsTable <- as.data.frame(summarize.data.frame(data.=osce.data.all,outcome.=as.factor(osce.data.all$Group),columns.=which(names(osce.data.all) %in% important.questions))$table)
byQuestionsTablePretty <- makePrettyTests(subset(byQuestionsTable,select=setdiff(names(byQuestionsTable),c("row","panel"))))


byQuestionsTable.all <- summarize.data.frame(data.=osce.data.all,outcome.=as.factor(osce.data.all$Group),columns.=which(names(osce.data.all) %in% important.questions))
byQuestionsTable.sig <- subset(as.data.frame(byQuestionsTable.all$table),panel %in% which(byQuestionsTable.all$pvalues<alpha01))

byQuestionsTable.sig.pretty <- makePrettyTests(subset(byQuestionsTable.sig,select=setdiff(names(byQuestionsTable.sig),c("row","panel"))))



pages <- Split.runs(byQuestionsTable$panel,21)


group.table2 <- list()
for(page in 1:max(pages)){
group.table2[[page]] <- simple.table.n(resultsdir.r2sas,resultsdir.sas,pastes("alldatasummary_sp_",page),as.matrix(byQuestionsTablePretty[pages==page,]),
				title=pastes("Table 3. Summarized standardized patient question data by group.",ifelse(page==1,"","cont'd")))
}

output.tables <- c(score.table,unlist(group.table),unlist(group.table2))

summaryTablesTable4SAS <- simple.multi.table(resultsdir.r2sas,resultsdir.sas,"group_comparison","group_comparison.rtf",output.tables,close.ods=TRUE)

write(summaryTablesTable4SAS,pastes(resultsdir.r2sas,"sasComm_1.txt"))






sig.table.video <- simple.table.n(resultsdir.r2sas,resultsdir.sas,pastes("alldatasummary_sp_sig_",page),as.matrix(byGroupTable.sig.pretty),
				title=pastes("Table 1. Suggestive results (p<0.1) from video analysis.",ifelse(page==1,"","cont'd")))

sig.table.spquestions <- simple.table.n(resultsdir.r2sas,resultsdir.sas,pastes("questions_sp_sig_",page),as.matrix(byQuestionsTable.sig.pretty),
				title=pastes("Table 2. Suggestive results (p<0.1) from standardized patient questions.",ifelse(page==1,"","cont'd")))


sig.table.output <- c(sig.table.video,sig.table.spquestions)
summaryTablesTable4SAS <- simple.multi.table(resultsdir.r2sas,resultsdir.sas,"group_comparison_sigtables","group_comparison_sigtables.rtf",sig.table.output,close.ods=TRUE)




numberMatrix0 <- as.matrix(osce.data.all[,c(importantVars,summary.vars,self.vars[-1],"Group")])

numberMatrix <- matrix(0,nrow(numberMatrix0),ncol(numberMatrix0))





for(i in 1:ncol(numberMatrix)){
	
	y <- numberMatrix0[,i]
	
	if(sum(is.na(as.numeric(y)))==length(y)){
	
		y <- as.factor(y)
	
	}
	
	numberMatrix[,i] <- as.numeric(y)
	
	}

dimnames(numberMatrix) <- dimnames(numberMatrix0)

numberMatrixFlip <- numberMatrix

for(meta.var.iter in unique(varData2$Meta.variable)){
	
	numberMatrixFlip.temp <- numberMatrix*0
	
	varData.temp <- subset(varData2,Meta.variable==meta.var.iter)
	
	print(meta.var.iter)

	for(j in intersect(colnames(numberMatrix),varData.temp$ShortName)){
			
			flip.code <- varData.temp[which(varData.temp$ShortName==j),]$A.is.good
			
			
			flip.factor <- ifelse(flip.code=="Yes",1,ifelse(flip.code=="No",-1,0))

			numberMatrixFlip.temp[,j] <- flip.factor					
			
			print(j)
			print(flip.code)
			print(flip.factor)

		
		}
	
	
	summary.stat <- as.matrix(rowSums(numberMatrix*numberMatrixFlip.temp,na.rm=TRUE))
	
	colnames(summary.stat) <- meta.var.iter
	
	numberMatrixFlip <- cbind(numberMatrixFlip,summary.stat)



} #loop over meta variables



meta.var.data <- data.frame(numberMatrixFlip[,unique(varData2$Meta.variable)],Group=osce.data.all[,"Group"],Case.Score=osce.data.all[,"Case.Score"])


by.meta.var.table <- as.data.frame(summarize.data.frame(data.=meta.var.data,outcome.=as.factor(meta.var.data$Group),columns.=which(names(meta.var.data) %in% make.names(unique(varData2$Meta.variable))))$table)
by.meta.var.table.pretty <- makePrettyTests(subset(by.meta.var.table,select=setdiff(names(by.meta.var.table),c("row","panel"))))

plot(subset(meta.var.data,select=setdiff(names(meta.var.data),"Group")))

par(mfrow=c(1,3))

for(meta.var.iter in make.names(unique(varData2$Meta.variable))){
	
	temp <- with(meta.var.data,cor.test(meta.var.data[[meta.var.iter]],Case.Score))
	
	with(meta.var.data,plot(meta.var.data[[meta.var.iter]],Case.Score))

	
	print(temp)
	
	}

meta.variable.table <- list()

meta.variable.table[[1]] <- simple.table.n(resultsdir.r2sas,resultsdir.sas,"meta_variable_vs_group",as.matrix(by.meta.var.table.pretty),
				title="Table 1. Metavariable summary vs. group assignment.")

summaryTablesTable4SAS <- simple.multi.table(resultsdir.r2sas,resultsdir.sas,"group_comparison_metavariables","group_comparison_metavariables.rtf",meta.variable.table,close.ods=TRUE)



colnames(numberMatrix)

sds <- apply(numberMatrix,2,sd,na.rm=TRUE)

numberMat2 <- numberMatrix[,sds>0]

corMat <- cor(numberMat2,use="pair")

all.cors <- NULL
ncols <- ncol(numberMat2)

cor.pvalue <- corMat*0+1

for(i in 1:ncols){
  for(j in 1:ncols){
    cor.pvalue[i,j] <- cor.test(numberMat2[,i],numberMat2[,j])$p.value
    	}
    	}




pdf(pastes(resultsdir,"key_correlations_between_variables.pdf"),height=11,width=11)


par(mfrow=c(1,1),mar=c(15, 4, 4, 2) + 0.1)
temp <- heatmap.2(corMat,col=bluered(256),margins=c(7,7),main="Positive Correlation Clusters")

carpet <- temp$carpet
par(mfrow=c(1,1),mar=c(12, 12, 4, 2) + 0.1)



cor.pvalue <- cor.pvalue[,match(colnames(carpet),colnames(cor.pvalue))]
cor.pvalue <- cor.pvalue[match(rownames(carpet),rownames(cor.pvalue)),]

label.matrix <- ifelse(cor.pvalue<0.05,round(carpet,2),"")

heatText(carpet,labelMatrix=label.matrix,col=bluered(256),main="Between Variable Positive Correlations",textSize=0.5)



par(mfrow=c(1,1),mar=c(15, 4, 4, 2) + 0.1)
temp <- heatmap.2(corMat,col=bluered(256),margins=c(7,7),main="Absolute Correlation Clusters",dist=corDist)

carpet <- temp$carpet

video.colorder <- colnames(carpet)

par(mfrow=c(1,1),mar=c(12, 12, 4, 2) + 0.1)



heatText(carpet,roundDigits=2,col=bluered(256),main="Variables Clustered by Correlations",textSize=0.5)

cor.pvalue <- cor.pvalue[,match(colnames(carpet),colnames(cor.pvalue))]
cor.pvalue <- cor.pvalue[match(rownames(carpet),rownames(cor.pvalue)),]

label.matrix <- ifelse(cor.pvalue<0.05,round(carpet,2),"")

importantVars
summary.vars
self.vars
group.vars <- "Group"

colors <- (video.colorder %in% summary.vars)*4+ (video.colorder %in% self.vars)*6 + (video.colorder %in% importantVars)*1+ (video.colorder %in% group.vars)*2


carpet2 <- carpet
rownames(carpet2) <- rep("",nrow(carpet2))
colnames(carpet2) <- rep("",ncol(carpet2))



heatText(carpet2,labelMatrix=label.matrix,col=bluered(256),main="Variables Clustered by Correlations",textSize=0.5)

hticks <- seq(0.,ncol(corMat)-1,by=1)/(ncol(corMat)-1)
vticks <- seq(0,ncol(corMat)-1,by=1)/(ncol(corMat)-1)

abline(h=hticks[seq(2,length(hticks),by=2)],lty=2,lwd=0.5,col="gray")
abline(v=hticks[seq(2,length(vticks),by=2)],lty=2,lwd=0.5,col="gray")

mtext(video.colorder,side=1,at=vticks,line=1.0,col=colors,las=2)
mtext(video.colorder,side=2,at=vticks,line=1.0,col=colors,las=2)

colNum <- ncol(carpet2)
rowNum <- nrow(carpet2)

y <- as.vector(matrix(hticks,colNum,rowNum))
x <- as.vector(matrix(vticks,colNum,rowNum,byrow=TRUE))
		
texter <- as.vector(t(label.matrix))
text(x,y,texter,cex=0.5)




label.matrix.p <- ifelse(cor.pvalue<0.05,p.rounder(cor.pvalue),"")

heatText(carpet2,labelMatrix=label.matrix.p ,col=bluered(256),main="pvalues Between Variable Absolute Correlations",textSize=0.25)


abline(h=hticks[seq(2,length(hticks),by=2)],lty=2,lwd=0.5,col="gray")
abline(v=hticks[seq(2,length(vticks),by=2)],lty=2,lwd=0.5,col="gray")

mtext(video.colorder,side=1,at=vticks,line=1.0,col=colors,las=2)
mtext(video.colorder,side=2,at=vticks,line=1.0,col=colors,las=2)

		
texter <- as.vector(t(label.matrix.p))
text(x,y,texter,cex=0.25)



dev.off()










par(mfrow=c(1,1),mar=c(15, 4, 4, 2) + 0.1)
temp <- heatmap.2(corMat,col=bluered(256),margins=c(7,7),main="Positive Correlation Clusters")

carpet <- temp$carpet
par(mfrow=c(1,1),mar=c(12, 12, 4, 2) + 0.1)



cor.pvalue <- cor.pvalue[,match(colnames(carpet),colnames(cor.pvalue))]
cor.pvalue <- cor.pvalue[match(rownames(carpet),rownames(cor.pvalue)),]

label.matrix <- ifelse(cor.pvalue<0.05,round(carpet,2),"")

heatText(carpet,labelMatrix=label.matrix,col=bluered(256),main="Between Variable Positive Correlations",textSize=0.5)



par(mfrow=c(1,1),mar=c(15, 4, 4, 2) + 0.1)
temp <- heatmap.2(corMat,col=bluered(256),margins=c(7,7),main="Absolute Correlation Clusters",dist=corDist)

carpet <- temp$carpet

video.colorder <- colnames(carpet)

par(mfrow=c(1,1),mar=c(12, 12, 4, 2) + 0.1)



heatText(carpet,roundDigits=2,col=bluered(256),main="Variables Clustered by Correlations",textSize=0.5)

cor.pvalue <- cor.pvalue[,match(colnames(carpet),colnames(cor.pvalue))]
cor.pvalue <- cor.pvalue[match(rownames(carpet),rownames(cor.pvalue)),]

label.matrix <- ifelse(cor.pvalue<0.05,round(carpet,2),"")

importantVars
summary.vars
self.vars
group.vars <- "Group"

colors <- (video.colorder %in% summary.vars)*4+ (video.colorder %in% self.vars)*6 + (video.colorder %in% importantVars)*1+ (video.colorder %in% group.vars)*2


carpet2 <- carpet
rownames(carpet2) <- rep("",nrow(carpet2))
colnames(carpet2) <- rep("",ncol(carpet2))


pdf(pastes(resultsdir,"key_correlations_between_variables_2.pdf"),height=11,width=11)


par(mfrow=c(1,1),mar=c(12, 12, 4, 2) + 0.1)

heatText(carpet2,labelMatrix=label.matrix,col=bluered(256),main="Variables Clustered by Correlations",textSize=0.43)

hticks <- seq(0.,ncol(corMat)-1,by=1)/(ncol(corMat)-1)
vticks <- seq(0,ncol(corMat)-1,by=1)/(ncol(corMat)-1)

abline(h=hticks[seq(2,length(hticks),by=2)],lty=2,lwd=0.5,col="gray")
abline(v=hticks[seq(2,length(vticks),by=2)],lty=2,lwd=0.5,col="gray")

mtext(video.colorder,side=1,at=vticks,line=1.0,col=colors,las=2)
mtext(video.colorder,side=2,at=vticks,line=1.0,col=colors,las=2)

colNum <- ncol(carpet2)
rowNum <- nrow(carpet2)

y <- as.vector(matrix(hticks,colNum,rowNum))
x <- as.vector(matrix(vticks,colNum,rowNum,byrow=TRUE))
		
texter <- as.vector(t(label.matrix))
text(x,y,texter,cex=0.43,col=ifelse(carpet>(-0.1),"black","white"))




label.matrix.p <- ifelse(cor.pvalue<0.05,p.rounder(cor.pvalue),"")

heatText(carpet2,labelMatrix=label.matrix.p ,col=bluered(256),main="pvalues Between Variable Absolute Correlations",textSize=0.25)


abline(h=hticks[seq(2,length(hticks),by=2)],lty=2,lwd=0.5,col="gray")
abline(v=hticks[seq(2,length(vticks),by=2)],lty=2,lwd=0.5,col="gray")

mtext(video.colorder,side=1,at=vticks,line=1.0,col=colors,las=2)
mtext(video.colorder,side=2,at=vticks,line=1.0,col=colors,las=2)

		
texter <- as.vector(t(label.matrix.p))
text(x,y,texter,cex=0.25)



dev.off()











patient.number.matrix <- as.matrix(osce.data.all[,important.questions])


sds <- apply(patient.number.matrix,2,sd,na.rm=TRUE)


patient.number.matrix <- patient.number.matrix[,sds>0]



corMat <- cor(patient.number.matrix,use="pair")

all.cors <- NULL
ncols <- ncol(patient.number.matrix )

cor.pvalue <- corMat*0+1

for(i in 1:ncols){
  for(j in 1:ncols){
    cor.pvalue[i,j] <- cor.test(patient.number.matrix[,i],patient.number.matrix[,j])$p.value
    	}
    	}






pdf(pastes(resultsdir,"correlations_between_patient.pdf"),height=8,width=11)



par(mfrow=c(1,1),mar=c(15, 4, 4, 2) + 0.1)
temp <- heatmap.2(corMat,col=bluered(256),margins=c(7,7),main="Patient Data Absolute Correlation Clusters",trace="none",dist=corDist,cexCol=.4,cexRow=0.35)

carpet <- temp$carpet
par(mfrow=c(1,1),mar=c(12, 12, 4, 2) + 0.1)


patient.colorder <- colnames(carpet)

cor.pvalue <- cor.pvalue[,match(colnames(carpet),colnames(cor.pvalue))]
cor.pvalue <- cor.pvalue[match(rownames(carpet),rownames(cor.pvalue)),]

label.matrix <- ifelse(cor.pvalue<0.05,round(carpet,2),"")


heatText(carpet,labelMatrix=label.matrix,col=bluered(256),main="Patient Data Between Variable Absolute Correlations",textSize=0.5)





dev.off()










combined.number.matrix0 <- osce.data.all[,c(video.colorder,patient.colorder)]

combined.number.matrix <- matrix(NA,nrow(combined.number.matrix0),ncol(combined.number.matrix0))

for(col.iter in 1:ncol(combined.number.matrix)){
	
	
	y <- combined.number.matrix0[,col.iter]
	
	if(sum(is.na(as.numeric(y)))==length(y)){
	
		y <- as.factor(y)
	
	}
	
	combined.number.matrix[,col.iter] <- as.numeric(y)

	
	
	
	
	}




dimnames(combined.number.matrix) <- dimnames(combined.number.matrix0)

sds <- apply(combined.number.matrix,2,sd,na.rm=TRUE)

combined.number.matrix<- combined.number.matrix[,sds>0]

corMat <- cor(combined.number.matrix,use="pair")

colnames(corMat)

all.cors <- NULL
ncols <- ncol(combined.number.matrix )

cor.pvalue <- corMat*0+1

for(i in 1:ncols){
  for(j in 1:ncols){
  	
  	try({
    cor.pvalue[i,j] <- cor.test(combined.number.matrix[,i],combined.number.matrix[,j])$p.value
    })
    	}
    	}







pdf(pastes(resultsdir,"correlations_between_patient_video.pdf"),height=8,width=11)



par(mfrow=c(1,1),mar=c(15, 4, 4, 2) + 0.1)


label.matrix <- ifelse(cor.pvalue<0.05,round(corMat,2),"")


heatText(corMat,axisTextSize=c(0.7,0.4),labelMatrix=label.matrix,col=bluered(256),main="Patient Scores vs. Video Variables Absolute Correlations",textSize=0.4)


abline(h=seq(0.,ncol(corMat)-1,by=2)/(ncol(corMat)-1),lty=2,lwd=0.2)
abline(v=seq(0,ncol(corMat)-1,by=2)/(ncol(corMat)-1),lty=2,lwd=0.2)


dev.off()



