##############################################################################
# University of Texas Health Science Center at San Antonio
# Department of Epidemiology and Biostatistics                        
##############################################################################
# Filename:   script_here                            
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


source.file <- gsub("\\.R","/","script_here.R")

# define the base directory	

basedir <- "/Users/jonathangelfond/Documents/Projects/"

# Creates necessary directories

analysisdir <- paste(basedir,"Programs/",sep="") # where the programs are
datadir <- paste(basedir,"Data/",sep="")  # where the data are
resultsdir <- paste(basedir,"Documents/Results/",source.file,sep="") # Standard output
tex.dir <- paste(resultsdir,"texdir/",sep="") # Publication quality output


#loads libraries and adds common functions to workspace
source(paste(analysisdir,"support_functions.R",sep="")) 


dir.create(resultsdir)
dir.create(tex.dir)

# Assign inpute filenames




table.files <- write.includer.vector(tex.list,path=tex.dir,include.file.base="includemetables")
make.latex.doc.vector(pastes(tex.dir,"summary.tex"),includer=table.files,title="Analysis",author="Jon Gelfond",date=as.character(Sys.time()),path=NULL)


print(paste("EOF:",source.file))






