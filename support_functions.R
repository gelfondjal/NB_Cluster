#
#  Support_functions.R
#
library(car)

lapply(read.table(paste(analysisdir,"common_libs.txt",sep=""),header=TRUE)[,1],library,character.only=TRUE)

source("/Users/jonathangelfond/Documents/Projects/sudheer/Analysis/table_functions.R")
source("/Users/jonathangelfond/Documents/Projects/OB_Residents/Brumm/Programs/support_functions.R")
chris.path <- "/Users/jonathangelfond/Dropbox/Work_Shared/R Scripts/Chris/"
chris.path2 <- "/Users/jonathangelfond/Documents/Projects/Louden/"

temp <- lapply(list.files(chris.path,full.names=1,pattern="r$"),source)
temp <- lapply(list.files(chris.path2,full.names=1,pattern="r$"),source)

cleaner <- function(x){
	
	# substitute +/- 5 for +/-Inf
	# substitute NA of nan
	
	x <- ifelse(abs(x)==Inf,5*sign(x),x)
	x <- ifelse(is.nan(x),NA,x)
	
	return(x)
	
	}




plot.mus.nmr.by.gene <- function(gene.pick,mus.d=mouse.liver.data.agg,nmr.d=nmr.liver.data.agg,plot.label=""){
	
	tall.data.nmr <- NULL
	
	try({
	
	nmr.gene <- subset(nmr.d,SYMBOL==gene.pick,select=names(nmr.d)[grep("\\.log",names(nmr.d))])
	tall.data <- data.frame(expression=t(as.matrix(nmr.gene)))
	tall.data <- data.frame(id=rownames(tall.data),expression=tall.data[,1])
	tall.data$Species <- "nmr"	
	tall.data.nmr <- tall.data
	
	
	
	tall.data <- NULL
	
	
	
	mus.gene <- subset(mus.d,SYMBOL==gene.pick,select=names(mus.d)[grep("\\.log",names(mus.d))])
	tall.data <- data.frame(expression=t(as.matrix(mus.gene)))
	tall.data <- data.frame(id=rownames(tall.data),expression=tall.data[,1])
	tall.data$Species <- "mus"
	

	
	if(is.null(tall.data)){return(NULL)}
	
	tall.data <- rbind(tall.data,tall.data.nmr)



	tall.data$Tx <- "DBMA"
	tall.data$Tx[grep("_C",tall.data$id)] <- "Placebo"
	
	tall.data$Tx <- reLeveler(as.factor(tall.data$Tx),c("Placebo","DBMA"))
	tall.data$Species <- as.factor(tall.data$Species)
	

	
	factor2.boxplot2("Tx","Species","expression",data=tall.data,main=paste(plot.label,gene.pick))
	
	
	
	mouse.p <- t.test(`expression`~`Tx`,data=subset(tall.data,Species=="mus"))
	nmr.p <- t.test(`expression`~`Tx`,data=subset(tall.data,Species=="nmr"))
	
	aov.result <- Anova(aov(`expression`~`Tx`*`Species`,data=tall.data))
	
	legend("topright",bg="white",legend=c(paste(c("mus","nmr","diff"),"P=",c(p.rounder(mouse.p$p.value),p.rounder(nmr.p$p.value),p.rounder(aov.result[["Pr(>F)"]][3])))))
	
	})
	
	}
	
	
	

interactions.mus.nmr.by.gene <- function(gene.pick,mus.d=mouse.liver.data.agg,nmr.d=nmr.liver.data.agg,plot.label=""){
	
	int.p.value <- NA
	
	tall.data.nmr <- NULL
	
	try({
	
	nmr.gene <- subset(nmr.d,SYMBOL==gene.pick,select=names(nmr.d)[grep("\\.log",names(nmr.d))])
	tall.data <- data.frame(expression=t(as.matrix(nmr.gene)))
	tall.data <- data.frame(id=rownames(tall.data),expression=tall.data[,1])
	tall.data$Species <- "nmr"	
	tall.data.nmr <- tall.data
	
	
	
	tall.data <- NULL
	
	
	
	mus.gene <- subset(mus.d,SYMBOL==gene.pick,select=names(mus.d)[grep("\\.log",names(mus.d))])
	tall.data <- data.frame(expression=t(as.matrix(mus.gene)))
	tall.data <- data.frame(id=rownames(tall.data),expression=tall.data[,1])
	tall.data$Species <- "mus"
	

	
	if(is.null(tall.data)){return(NULL)}
	
	tall.data <- rbind(tall.data,tall.data.nmr)



	tall.data$Tx <- "DBMA"
	tall.data$Tx[grep("_C",tall.data$id)] <- "Placebo"
	
	tall.data$Tx <- reLeveler(as.factor(tall.data$Tx),c("Placebo","DBMA"))
	tall.data$Species <- as.factor(tall.data$Species)
	

	
#	factor2.boxplot2("Tx","Species","expression",data=tall.data,main=paste(plot.label,gene.pick))
	
	
	
	mouse.p <- t.test(`expression`~`Tx`,data=subset(tall.data,Species=="mus"))
	nmr.p <- t.test(`expression`~`Tx`,data=subset(tall.data,Species=="nmr"))
	
	aov.result <- Anova(aov(`expression`~`Tx`*`Species`,data=tall.data))
	
	int.p.value <- aov.result[["Pr(>F)"]][3]
	
	
	})
	
	return(int.p.value)
	
	
	}	
	
	
	
	
	
	
	loess.surf <- function (Y, X, span = 0.75, degree = 1, family = "gaussian", 
    phi = 20, theta = 50, xlab = "X", ylab = "Y", zlab = "Fit", 
    line.col = 1, line.type = 1, scale = TRUE, duplicate = "error", 
    expand = 0.5, ...) 
{
    require(akima)
    lom <- loess(Y ~ X, span = span, degree = degree, family = family)
    r <- ncol(lom$x)
    yhat <- lom$fitted
    if (r > 2) 
        stop("Number of predictors must be < 2")
    if (r == 1) {
        plot(X, Y, xlab = xlab, ylab = ylab, ...)
        d <- 10^round(log(min(X), 10), 0)
        xv <- seq(min(X), max(X), d/100)
        yv <- predict(lom, data.frame(X = xv))
        lines(xv, yv, col = line.col, lty = line.type)
    }
    if (r == 2) {
        iout <- c(1:length(yhat))
        nm1 <- length(yhat) - 1
        for (i in 1:nm1) {
            ip1 <- i + 1
            for (k in ip1:length(yhat)) if (sum(X[i, ] == X[k, 
                ]) == 2) 
                iout[k] <- 0
        }
        yhat1 <- yhat[iout >= 1]
        X1 <- X[iout >= 1, ]
        fit <- interp(X1[, 1], X1[, 2], yhat1, duplicate = duplicate)
        persp(fit, theta = theta, phi = phi, xlab = xlab, ylab = ylab, 
            zlab = zlab, scale = scale, expand = expand)
    }
}





# Support functions

grep.to.label <- function(char.vec,pattern,withlabel,withoutlabel){
	
	y <- char.vec
	y[grep(pattern,char.vec)] <- withlabel
	y[-grep(pattern,char.vec)] <- withoutlabel
	
	return(y)
	
	}



compute.cluster.correlations <- function(x,e.mat=expression.mat){
	
	# use within ddply to compute the correation within clusters defined by x
	
	match.genes <- unique(x$V2)
	
#	print(match.genes)
	
	e.mat.2 <- e.mat[match.genes,]
	
#	print(e.mat.2)
	
	little.mat <- t(e.mat.2)
	
	cor.mat <- cor(little.mat,use="pair")
	
	diag(cor.mat) <- NA
	
	mean.cor <- mean(cor.mat,na.rm=TRUE)
	
	match.count <- min(x$Count.Freq)
	
	return(data.frame(Count.Freq=match.count,mean.cor=mean.cor))
	
	} # END: compute.cluster.correlations
	
	

compute.cluster.correlations.new <- function(x,e.mat=expression.mat){
	
	# use within ddply to compute the correation within clusters defined by x
	
	match.genes <- unique(x$NMR.gene)
	
#	print(match.genes)
	
	
	match.rows <- (rownames(e.mat) %in% match.genes)& !(duplicated(rownames(e.mat)))
	
	
	e.mat.2 <- e.mat[match.rows,]
	
#	print(e.mat.2)
	
	little.mat <- t(e.mat.2)
	
	cor.mat <- cor(little.mat,use="pair")
	
	diag(cor.mat) <- NA
	
	mean.cor <- mean(cor.mat,na.rm=TRUE)
	
	match.count <- min(x$Count.Freq)
	
	return(data.frame(Count.Freq=match.count,mean.cor=mean.cor))
	
	} # END: compute.cluster.correlations
	




empirical.p.function <- function(x,l.hat,a.hat){
	
		return(l.hat*x+(1-l.hat)*x^a.hat)
	
}


p.bum.correct <- function(p.original,p.permuted){
	
		# Fits Beta-Uniform Mixture to p-permuted p-values 
		# adjusts the p.orginal
		
		
		p.permuted <- as.numeric(na.exclude(p.permuted))
		
		delta <- 1e-128
		
		p.permuted <- ifelse(p.permuted>delta,p.permuted,delta)
	
		bum.perm <- Bum(p.permuted)
	
		
		alpha.hat <- bum.perm@ahat
		l.hat <- bum.perm@lhat
		pi0.hat <- bum.perm@pihat

		

		return(empirical.p.function(p.original,l.hat,alpha.hat))
	
	
	} #END p.bum.correct
	
to.char.matrix <- function(x,rounder=2){
	
	# adds rownames and column names to make a character matrix of x
	# rounder is the number of digits
	# x is a numeric matrix
	
	if(rounder>=0){x <- round(x,rounder)}
	
	x.out <- rbind(colnames(x),x)
	
	x.out <- cbind(c("",rownames(x)),x.out)
	
	return(x.out)
	
	}	#end to.char.matrix
	

p.bum.correct.right.truncate <- function(p.original,p.permuted,censor=0.3){
	
		# Fits Beta-Uniform Mixture to p-permuted p-values 
		# adjusts the p.orginal
		
		
		p.permuted <- as.numeric(na.exclude(p.permuted))
				
		delta <- 1e-128
		
		p.permuted <- ifelse(p.permuted>delta,p.permuted,delta)
		
		p.permuted <- p.permuted[p.permuted<censor]
		
		p.permuted <- c(p.permuted,seq(censor,1.0,length=ceiling(length(p.permuted)*(1-censor))))
		
		hist(p.permuted)
	
		bum.perm <- Bum(p.permuted)
	
		
		alpha.hat <- bum.perm@ahat
		l.hat <- bum.perm@lhat
		pi0.hat <- bum.perm@pihat

		

		return(empirical.p.function(p.original,l.hat,alpha.hat))
	
	
	} #END p.bum.correct	
