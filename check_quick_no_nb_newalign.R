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
# v 0.2 add weighted glmm
# v 0.3 add permuations
##############################################################################


set.seed(2013)
rm(list=ls())
require(ggplot2)
require(DESeq)
library(reshape2)
library(lme4)
library(ClassComparison)
library(glmmADMB)
library(biomaRt)

source.file <- gsub("\\.R","/","check_gee_analysis_convergence_quick_no_nb_newalign.R")
source0 <- source.file
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
workspace.file <- "/Users/jonathangelfond/Documents/Projects/Buffenstein/NMR_Genomics/Documents/Results/check_gee_analysis_1/all.R.Obj.Rdata"

#save(list=ls(),file=workspace.file)
load(file.path(resultsdir,"results_workspace.RData"))


source.file <- source0
basedir <- "/Users/jonathangelfond/Documents/Projects/Buffenstein/NMR_Genomics/"

# Creates necessary directories

analysisdir <- paste(basedir,"Programs/",sep="") # where the programs are
datadir <- paste(basedir,"Data/Hibbs/countFiles/",sep="")  # where the data are
resultsdir <- paste(basedir,"Documents/Results/",source.file,sep="") # Standard output
tex.dir <- paste(resultsdir,"texdir/",sep="") # Publication quality output




source("/Users/jonathangelfond/Documents/Projects/Joe_Lab/NMR_Methods/Programs/support_functions.R")
source("/Users/jonathangelfond/Documents/Projects/sudheer/Analysis/table_functions.R")

library(geepack)
library(gee)
library(MGLM)




#install.packages("gplots")
#install.packages("blockmodeling")
#install.packages("/Users/jonathangelfond/Documents/Projects/Buffenstein/NMR_Genomics/Documents/Background/EBSeq_1.1.6.tar.gz", repos=NULL, type="source")
library(gplots)
library(blockmodeling)
library(EBSeq)





# merge the datasets for analysis

single.data  <- data.frame(type="single",single.analysis.out,stringsAsFactors=FALSE)

single.data$NMR.gene <- single.data$V1

single.data <- Drop.variable(single.data,"V1")

merged.single <- merge(single.data,matches.2,by="NMR.gene")

merged.single$total.matches <- merged.single$alignment.length * merged.single$proportion.homology

merged.single <- Sortby(merged.single,"proportion.homology",decreasing=TRUE)

merged.single.nodups <- subset(merged.single,!duplicated(mus.UCSCID))

merged.single.nodups$pvalue.perm <- ecdf(single.analysis.out.perm[,"pvalue"])(merged.single.nodups$pvalue)

analyses <- lapply(analysis.out,function(x){x[[1]]})

bigger.list <- list()

for(namer in names(analyses)){

	bigger.list[[namer]] <- data.frame(type=namer,analyses[[namer]],stringsAsFactors=FALSE)	
	
	bigger.list[[namer]]$mus.UCSCID <- bigger.list[[namer]]$V1
	
	bigger.list[[namer]] <- merge(bigger.list[[namer]],subset(matches.2,!duplicated(mus.UCSCID)),by="mus.UCSCID")	
	
	bigger.list[[namer]]$pvalue.perm <- ecdf(perm.combo[[namer]]$perm[,"pvalue"])(bigger.list[[namer]]$pvalue)
	
}

#save(list=ls(),file=)

all.results <- rbind.fill(rbind.fill(bigger.list),merged.single.nodups)

all.types.list <- unique(c(as.character(subset(bigger.list[[1]],!is.na(pvalue))$mus.UCSCID),as.character(subset(bigger.list[[2]],!is.na(pvalue))$mus.UCSCID)))

all.results <- subset(all.results,mus.UCSCID %in% all.types.list)

alpha <-0.05

summary.results <- ddply(all.results,c("type","Count.Freq"),function(x){
	
			x$mean.power <- mean(x$pvalue.perm<alpha,na.rm=TRUE)
			x$sum.sig <- sum(x$pvalue.perm<alpha,na.rm=TRUE)
			x$total <- sum(!is.na(x$pvalue.perm))
	
			return(subset(x[1,],select=c("mean.power","sum.sig","total")))
	
	})


tex.list <- list()

summary.results <- subset(summary.results,Count.Freq<=5)

summary.results[["Cluster.size"]] <- summary.results[["Count.Freq"]]

summary.results.2 <- Drop.variable(summary.results,"Count.Freq")


sig.w.mvariate <- subset(all.results,(pvalue.perm<alpha)&(type=="GEE"))$V1
sig.w.singleton <- subset(all.results,(pvalue.perm<alpha)&(type=="single"))$V1

multivariate.only.not.singleton <- setdiff(sig.w.mvariate,sig.w.singleton)

multivariate.only.not.singleton.results <- subset(all.results,(V1 %in% multivariate.only.not.singleton)&(type=="GEE"))

multivariate.only.not.singleton.results <- Sortby(multivariate.only.not.singleton.results ,"pvalue")

multivariate.only.not.singleton.results <- rename(multivariate.only.not.singleton.results,replace=c("V1" ="ucsc","NMR.gene" ="NMR_Gene"))


mart <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
attributes <- c( "chromosome_name",
                "start_position","end_position", "ensembl_gene_id", 
                "external_gene_id", "description","ucsc")
                
genes <- getBM(attributes=attributes, filters="ucsc", 
               values=multivariate.only.not.singleton, mart=mart, uniqueRows=T)


merged.list <- merge(multivariate.only.not.singleton.results,genes,by="ucsc")

merged.list <-Sortby(merged.list,"pvalue")

write.csv(merged.list,file.path(resultsdir,"liver_to_skin_differences_unique_to_new_model.csv"),row.names=FALSE)




tex.list <- lcat(tex.list,latex.out.table.2(Norow.mat(round.data(summary.results.2,2)),
		caption=paste("Liver vs Skin in NMR alpha=",alpha,"Number significant after permutation")))
		
		
summary.results.2$se <- sqrt(summary.results.2$mean.power*(1-summary.results.2$mean.power)/summary.results.2$total)
		
library(ggplot2)	

summary.results.2$Type <- ifelse(summary.results.2$type=="single","Singelton",summary.results.2$type)
summary.results.2$Type <- ifelse(summary.results.2$type=="quasipossion","Quasipoisson",summary.results.2$Type)


temp <- ggplot(summary.results.2, aes(x=as.factor(summary.results.2$Cluster.size), y=summary.results.2$mean.power,colour=summary.results.2$Type)) + 
		    geom_errorbar(aes(ymin=summary.results.2$mean.power-2*summary.results.2$se, ymax=summary.results.2$mean.power+2*summary.results.2$se), width=.1,position=position_dodge(.1)) +
		    geom_line(position=position_dodge(.1)) +
		    geom_point(position=position_dodge(.1)) + xlab("Cluster Size") + ylab("Proportion Significant")+
			 ggtitle(paste("Proportion Significant","vs","Cluster Size"))+ scale_colour_hue(summary.results.2$Type,name="Model") + theme_bw()
	
print(temp)	
		

table.files <- write.includer.vector(tex.list,path=tex.dir,include.file.base="includemetables")
make.latex.doc.vector(pastes(tex.dir,"summary_realdata_newalign.tex"),includer=table.files,title="Liver vs. Skin NMR to mouse homologs",author="Jon Gelfond",date=as.character(Sys.time()),path=NULL)

Create.pdflatex("summary_realdata_newalign.tex",tex.dir)




print(paste("EOF:",source.file))






