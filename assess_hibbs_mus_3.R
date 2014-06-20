##############################################################################
# University of Texas Health Science Center at San Antonio
# Department of Epidemiology and Biostatistics                        
##############################################################################
# Filename:   assess_hibbs_mus_3.R                            
# Author:    Jon Gelfond                                             
# Project Name:  
# Input:       
# Output:   Source Filename directory
#
# Modification History:
# v 0.1 Creation 
# v 3 compares 3 mappings mus, hs, cp
##############################################################################


set.seed(2013)
#rm(list=ls())


source.file <- gsub("\\.R","/","assess_hibbs_mus_3.R")

# define the base directory	

basedir <- "/Users/jonathangelfond/Documents/Projects/Buffenstein/NMR_Genomics/"

# Creates necessary directories

analysisdir <- paste(basedir,"Programs/",sep="") # where the programs are
datadir <- paste(basedir,"Data/Hibbs/",sep="")  # where the data are
resultsdir <- paste(basedir,"Documents/Results/",source.file,sep="") # Standard output
tex.dir <- paste(resultsdir,"texdir/",sep="") # Publication quality output


#loads libraries and adds common functions to workspace
source(paste(analysisdir,"support_functions.R",sep="")) 


dir.create(resultsdir)
dir.create(tex.dir)

# Assign inpute filenames

mus.all <- read.csv(pastes(datadir,"MMU.genes.results.csv"),as.is=TRUE)

mus.all$Symbol <- toupper(mus.all$Symbol)

mus.data.matrix <- as.matrix(subset(mus.all,select=grep("(Skin.*[0-9])|(Liver.*[0-9])",names(mus.all))))
rownames(mus.data.matrix) <- as.character(mus.all$Symbol)

nmr.to.mus <- read.csv(pastes(datadir,"NMR2mouse.genes.results.csv"),as.is=TRUE)
nmr.to.mus$Symbol <- toupper(nmr.to.mus$Symbol)
nmr.to.mus.data.matrix <- as.matrix(subset(nmr.to.mus,select=grep("(Skin.*[0-9])|(Liver.*[0-9])",names(nmr.to.mus))))
rownames(nmr.to.mus.data.matrix) <- as.character(nmr.to.mus$Symbol)


nmr.to.hs <- read.csv(pastes(datadir,"NMR2human.genes.results.csv"),as.is=TRUE)
nmr.to.hs$Symbol <- toupper(nmr.to.hs$Symbol)
nmr.to.hs.data.matrix <- as.matrix(subset(nmr.to.hs,select=grep("(Skin.*[0-9])|(Liver.*[0-9])",names(nmr.to.hs))))
rownames(nmr.to.hs.data.matrix) <- as.character(nmr.to.hs$Symbol)

nmr.to.cp <- read.csv(pastes(datadir,"NMR2guineapig.genes.results.csv"),as.is=TRUE)
nmr.to.cp$Symbol <- toupper(nmr.to.cp$MGI.symbol)
nmr.to.cp.data.matrix <- as.matrix(subset(nmr.to.cp,select=grep("(Skin.*[0-9])|(Liver.*[0-9])",names(nmr.to.cp))))
rownames(nmr.to.cp.data.matrix) <- as.character(nmr.to.cp$Symbol)


class(mus.data.matrix) <- "numeric"
class(nmr.to.mus.data.matrix) <- "numeric"
class(nmr.to.hs.data.matrix) <- "numeric"
class(nmr.to.cp.data.matrix) <- "numeric"

# Remove duplicates



common.genes <- unique(as.vector(as.matrix(read.csv("/Users/jonathangelfond/Documents/Projects/Thompson/Results/common_genes.csv"))[,2]))

length(unique(common.genes))

mean(common.genes %in% nmr.to.mus$Symbol)

mean(common.genes %in% nmr.to.cp$Symbol)


B.network.read <- as.matrix(read.csv("/Users/jonathangelfond/Documents/Projects/Buffenstein/NMR_Genomics/Documents/Results/analyze_liver_v0.1/mouse_network_nmr_seq.csv"))
B.network <- B.network.read[,-1]

rownames(B.network) <- sort(common.genes)
colnames(B.network) <- sort(common.genes)
B.base.network <- B.network


common.genes.2 <- intersect(common.genes,mus.all$Symbol)
common.genes.2 <- intersect(common.genes.2,nmr.to.mus$Symbol)

B.network <- B.network[,colnames(B.network) %in% common.genes.2]
B.network <- B.network[rownames(B.network) %in% common.genes.2,]


mouse.data.pick <- mus.data.matrix[rownames(mus.data.matrix) %in% common.genes.2,]
nmr.data.pick <- nmr.to.mus.data.matrix[rownames(nmr.to.mus.data.matrix) %in% common.genes.2,]


mouse.data.pick <- mouse.data.pick[!duplicated(rownames(mouse.data.pick)),]
nmr.data.pick <- nmr.data.pick[!duplicated(rownames(nmr.data.pick)),]


mouse.data.pick <- mouse.data.pick[order(rownames(mouse.data.pick)),]
nmr.data.pick <- nmr.data.pick[order(rownames(nmr.data.pick)),]

nmr.t <- apply(nmr.data.pick,1,scale,scale=FALSE)
mouse.t <- apply(mouse.data.pick,1,scale,scale=FALSE)


nmr.t[Inf==nmr.t] <- NA
mouse.t[Inf==mouse.t] <- NA

nmr.t <- t(impute.knn(t(nmr.t))$data)
mouse.t <- t(impute.knn(t(mouse.t))$data)


cors <- matrix(NA,length(common.genes.2),2)
colnames(cors) <- c("mouse","nmr")
rownames(cors) <- common.genes.2


n.lambda <- 1


for(l in 1:n.lambda){
		
yhat.mouse <- 	B.network %*% t(mouse.t)

yhat.nmr <- 	B.network %*% t(nmr.t)

	
for(g in 1:length(common.genes.2)){
			
		
cors[g,"mouse"] <- cor(yhat.mouse[g,],mouse.t[,g],method="spearman")
cors[g,"nmr"] <- cor(yhat.nmr[g,],nmr.t[,g],method="spearman")
			
}
}


yhat.nmr.mus <- yhat.nmr

nmr.t.mus <- nmr.t




osrr.cors.mouse <- data.frame(cors)


medians <- colMedians(cors,na.rm=TRUE)
names(medians) <- colnames(cors)

pdf(pastes(resultsdir,"mouse_network_comparison.pdf"),h=8,w=8)

hist(cors[,"nmr"])
legend("topleft",legend=paste("Median Prediction Cor",round(medians["nmr"],2)))

hist(cors[,"mouse"])
legend("topleft",legend=paste("Median Prediction Cor",round(medians["mouse"],2)))


plot(cors[,"nmr"],cors[,"mouse"],main=paste("Mouse vs. NMR: Correlation with Prediction in",length(common.genes.2)," genes"))
abline(v=0,h=0,lty=1,lwd=2,col=2)


text(0.5,1,"Good Match",col=2,cex=2)
text(-0.5,1,"NMR Different",col=2,cex=2)




dev.off()

common.genes.3 <- intersect(common.genes,mus.all$Symbol)
common.genes.3 <- intersect(common.genes.3,nmr.to.hs$Symbol)

B.network <- B.base.network[,colnames(B.base.network) %in% common.genes.3]
B.network <-B.network[rownames(B.base.network) %in% common.genes.3,]


mouse.data.pick <- mus.data.matrix[rownames(mus.data.matrix) %in% common.genes.3,]

nmr.data.pick <- nmr.to.hs.data.matrix[rownames(nmr.to.hs.data.matrix) %in% common.genes.3,]


mouse.data.pick <- mouse.data.pick[!duplicated(rownames(mouse.data.pick)),]
nmr.data.pick <- nmr.data.pick[!duplicated(rownames(nmr.data.pick)),]


mouse.data.pick <- mouse.data.pick[order(rownames(mouse.data.pick)),]
nmr.data.pick <- nmr.data.pick[order(rownames(nmr.data.pick)),]

nmr.t <- apply(nmr.data.pick,1,scale,scale=FALSE)
mouse.t <- apply(mouse.data.pick,1,scale,scale=FALSE)


nmr.t[Inf==nmr.t] <- NA
mouse.t[Inf==mouse.t] <- NA

nmr.t <- t(impute.knn(t(nmr.t))$data)
mouse.t <- t(impute.knn(t(mouse.t))$data)

nmr.t.hs <- nmr.t

n.test.genes <- length(common.genes.3)

cors <- matrix(NA,n.test.genes,2)
colnames(cors) <- c("mouse","nmr")
rownames(cors) <- common.genes.3

pcors <- cors

resids <- matrix(NA,nrow(mouse.t),n.test.genes)


	
lambdas <- matrix(NA,n.test.genes,n.lambda)

colnames(lambdas) <- lambda

pheno.types <- data.frame(row=colnames(nmr.data.pick))
pheno.types$Tissue <- "Skin"
pheno.types$Tissue[grep("^Liver",pheno.types$row)] <- "Liver"
pheno.types$Tx <- "Drug"
pheno.types$Tx[grep("C",pheno.types$row)] <- "Control"


g=1
lmfit <- lm(nmr.t[,g]~yhat.nmr[g,]+Tissue*Tx,data=pheno.types)			
names <- rownames(Anova(lmfit))
pvals <-  matrix(NA,n.test.genes,2)
pvals2 <- matrix(NA,g,length(names))


for(l in 1:n.lambda){
		
yhat.mouse <- 	B.network %*% t(mouse.t)

yhat.nmr <- 	B.network %*% t(nmr.t)

	
for(g in 1:length(common.genes.3)){

outp <- NA
try({			
lmfit <- lm(nmr.t[,g]~Tissue*Tx,data=pheno.types)			
outp <- Anova(lmfit)["Tx","Pr(>F)"]
})

pvals[g,1] <- outp


outp <- NA
try({			
lmfit <- lm(nmr.t[,g]~yhat.nmr[g,]+Tissue*Tx,data=pheno.types)			
outp <- Anova(lmfit)["Tx","Pr(>F)"]
})

pvals[g,2] <- outp
		
cors[g,"mouse"] <- cor(yhat.mouse[g,],mouse.t[,g],method="spearman")
cors[g,"nmr"] <- cor(yhat.nmr[g,],nmr.t[,g],method="spearman")


resids[,g] <- lm(nmr.t[,g]~yhat.nmr[g,])$resid

pcors[g,"mouse"] <- cor(yhat.mouse[g,],mouse.t[,g],method="pearson")
pcors[g,"nmr"] <- cor(yhat.nmr[g,],nmr.t[,g],method="pearson")


			
}
}



vars <- apply(resids,2,var,na.rm=TRUE)


osrr.cors.nmr <- data.frame(cors)


medians <- colMedians(cors,na.rm=TRUE)
names(medians) <- colnames(cors)

pdf(pastes(resultsdir,"hs_mapping_network_comparison.pdf"),h=8,w=8)

hist(cors[,"nmr"])
legend("topleft",legend=paste("Median Prediction Cor",round(medians["nmr"],2)))

hist(cors[,"mouse"])
legend("topleft",legend=paste("Median Prediction Cor",round(medians["mouse"],2)))


plot(cors[,"nmr"],cors[,"mouse"],main=paste("Mouse vs. NMR: Correlation with Prediction in",length(common.genes.3)," genes"))
abline(v=0,h=0,lty=1,lwd=2,col=2)


text(0.5,1,"Good Match",col=2,cex=2)
text(-0.5,1,"NMR Different",col=2,cex=2)

plot(osrr.cors.nmr[,"nmr"],log(vars),main="Resid Variance vs. Correlation")



dev.off()


temp1 <- osrr.cors.mouse
temp2 <- osrr.cors.nmr

names(temp1) <- pastes("mus.map.",names(temp1))
names(temp2) <- pastes("hs.map.",names(temp2))

temp1$Symbol <- rownames(temp1)
temp2$Symbol <- rownames(temp2)

merged.mappings <- merge(temp1,temp2,by="Symbol")


with(merged.mappings,t.test(mus.map.nmr-hs.map.nmr))

pdf(pastes(resultsdir,"mus_to_hs_mapping_comparison.pdf"))

with(merged.mappings,plot(mus.map.nmr,hs.map.nmr,main="Human vs. Mouse Network Predictions"))

text(c(-0.75,0.75,0.75),c(1,1,-1),c("HS Better","Both","Mus Better"),col=2,cex=2)


with(merged.mappings,hist(mus.map.nmr-hs.map.nmr,breaks=20))

plot(pcors[,"nmr"],cors[,"nmr"],xlab="Pearson",ylab="Spearman",main="NMR to Human Spearman vs. Pearson")


dev.off()



common.genes.4 <- intersect(common.genes,unique(nmr.to.cp$Symbol))

B.network <- B.base.network[,colnames(B.base.network) %in% common.genes.4]
B.network <-B.network[rownames(B.base.network) %in% common.genes.4,]


mean(nmr.to.cp$Symbol %in% rownames(nmr.to.cp.data.matrix))

nmr.data.pick <- nmr.to.cp.data.matrix[rownames(nmr.to.cp.data.matrix) %in% common.genes.4,]


nmr.data.pick <- nmr.data.pick[!duplicated(rownames(nmr.data.pick)),]


nmr.data.pick <- nmr.data.pick[order(rownames(nmr.data.pick)),]

nmr.t <- apply(nmr.data.pick,1,scale,scale=FALSE)


nmr.t[Inf==nmr.t] <- NA

nmr.t <- t(impute.knn(t(nmr.t))$data)


n.test.genes <- length(common.genes.4)

cors <- matrix(NA,n.test.genes,2)
colnames(cors) <- c("mouse","nmr")
rownames(cors) <- common.genes.4

pcors <- cors

resids <- matrix(NA,nrow(mouse.t),n.test.genes)



pheno.types <- data.frame(row=colnames(nmr.data.pick))
pheno.types$Tissue <- "Skin"
pheno.types$Tissue[grep("^Liver",pheno.types$row)] <- "Liver"
pheno.types$Tx <- "Drug"
pheno.types$Tx[grep("C",pheno.types$row)] <- "Control"


g=1
lmfit <- lm(nmr.t[,g]~yhat.nmr[g,]+Tissue*Tx,data=pheno.types)			
names <- rownames(Anova(lmfit))
pvals <-  matrix(NA,n.test.genes,2)
pvals2 <- matrix(NA,n.test.genes,length(names))


for(l in 1:n.lambda){
		

yhat.nmr <- 	B.network %*% t(nmr.t)

	
for(g in 1:length(common.genes.3)){

outp <- NA
try({			
lmfit <- lm(nmr.t[,g]~Tissue*Tx,data=pheno.types)			
outp <- Anova(lmfit)["Tx","Pr(>F)"]
})

pvals[g,1] <- outp


outp <- NA
try({			
lmfit <- lm(nmr.t[,g]~yhat.nmr[g,]+Tissue*Tx,data=pheno.types)			
outp <- Anova(lmfit)["Tx","Pr(>F)"]
pvals2[g,] <- Anova(lmfit)[,"Pr(>F)"]
})

pvals[g,2] <- outp
		
cors[g,"nmr"] <- cor(yhat.nmr[g,],nmr.t[,g],method="spearman")


resids[,g] <- lm(nmr.t[,g]~yhat.nmr[g,])$resid

pcors[g,"nmr"] <- cor(yhat.nmr[g,],nmr.t[,g],method="pearson")


			
}
}


colnames(pvals2) <- names

plot(pcors[,"nmr"],-log10(pvals2[,3]),main="DE p-value vs. Prediction")
plot(pcors[,"nmr"],-log10(pvals2[,3]),main="DE p-value vs. Prediction")

boxplot(-log10(pvals2[,2])~as.numeric(cut(pcors[,"nmr"],10)),main="DE p-value vs. Prediction",ylab="-log10(p)",xlab="Correlation Quantile")

pdf(pastes(resultsdir,"network_vs_de.pdf"),h=8,w=8)
par(mfrow=c(2,2))
for(i in setdiff(names,"Residuals")){

pvalueQQplot(pvals2[,i],title=paste(i,"p-value"))

}
dev.off()


pdf(pastes(resultsdir,"network_vs_de_cors.pdf"),h=8,w=8)
par(mfrow=c(2,2))
for(i in setdiff(names,"Residuals")){

boxplot(-log10(pvals2[,i])~as.numeric(cut(pcors[,"nmr"],10)),main=paste(i,"DE p-value vs. Prediction"),ylab="-log10(p)",xlab="Correlation Quantile")

}
dev.off()




hist(pcors[,"nmr"])

pdf(pastes(resultsdir,"cp_mapping_network_comparison.pdf"),h=8,w=8)

hist(pcors[,"nmr"])
legend("topleft",legend=paste("Median Prediction Cor",round(median(pcors[,"nmr"],na.rm=TRUE),2)))

dev.off()


print(paste("EOF:",source.file))






