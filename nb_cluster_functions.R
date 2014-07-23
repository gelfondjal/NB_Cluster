#
#  Support_functions.R
#

library(MASS)
#library(glmmADMB)
library(geepack)


fit.homology.cluster <- function(tissue.counts,match.info,pheno.factor,crash.list=NULL,weight=FALSE,family=poisson,permnumber=4){
	
	#calls fit.gee.5.array
	
	#tissue.counts = matrix with rows are transcripts and columns are samples
	#match.info = cluster identities columns are V1[cluster identitiy] V2[target transcript matches rows in tissue.counts] V3[homology weights]
	#pheno.factor = factor representing differential expression analysis grouping (e.g., treated or controls)
	#crash.list = the clusters that will not be considered in analysis
	#weight = TRUE?FALSE whether to perform homology weighted analysis
	#family = poisson, negbinom1, quasipoisson what type of cluster model to fit
	#permnumber = how many permutations of full transcriptome (nrows(tissue.counts)*permnumber)
	
	orginal <- fit.gee.5.array(tissue.counts,match.info,pheno.factor,crash.list,weight,perm=0,family)
	permuted <- fit.gee.5.array(tissue.counts,match.info,pheno.factor,crash.list,weight,perm=permnumber,family)

	oringal$corrected.pvalue <- ecdf(permuted[[1]][,"pvalue"])(original[[1]][,"pvalue"])

	
	#returns list of length 2
	# 1-data.fram with columns V1[cluster identity] FC[fold-change] pvalue convergence [of cluster model] corrected.pvalue [permutation based pvalue]
	# 2 list of fit objects corresponding to each cluster
	
	return(original)

	
	}



fit.gee.5.array <- function(tissue.counts,match.info,pheno.factor,crash.list=NULL,weight=FALSE,perm=0,family=poisson){

	
	#tissue.counts = matrix with rows are transcripts and columns are samples
	#match.info = cluster identities columns are V1[cluster identitiy] V2[target transcript matches rows in tissue.counts] V3[homology weights]
	#pheno.factor = factor representing differential expression analysis grouping (e.g., treated or controls)
	#crash.list = the clusters that will not be considered in analysis
	#weight = TRUE?FALSE whether to perform homology weighted analysis
	#family = poisson, negbinom1, quasipoisson what type of cluster model to fit
	#perm = if nonzero how many permutations of full transcriptome (nrows(tissue.counts)*permnumber)
	
	
	
	cluster.ids <- match.info$V1

	n.clusters <- length(unique(cluster.ids))

	cluster.names <- unique(cluster.ids)
	clust.iter <- cluster.names[1]
	
	if(!identical(rownames(tissue.counts),match.info$V2)){stop("Cluster ids don't match expression matrix")}

	results.fc <- matrix(NA,n.clusters,3)
	colnames(results.fc) <- c("FC","pvalue","conv")
	rownames(results.fc) <- unique(cluster.names)

	base.time <- timer.out()

	#crash.list.vals <- rownames(results.fc)[crash.list]
	
	cluster.row.fits <- rownames(results.fc)
	
	if(!is.null(crash.list)){cluster.row.fits <- rownames(results.fc)[crash.list]}

	#clust.iter <- hits[57]
	
	pheno.factor.0 <- pheno.factor
	
	gee.list <- list()
	
	
	perm.results <- NULL
	
	for(perm.iter in 1:ifelse(perm==0,1,perm)){
	
	for(clust.iter in cluster.row.fits){


		if(perm>0){pheno.factor <- pheno.factor.0[sample(length(pheno.factor.0))]}
		
		pheno.data <- data.frame(Var1=colnames(tissue.counts), pheno.factor,stringsAsFactors=0)

#		gene.data <- Sortby(data.frame(melt(t(tissue.counts[cluster.ids==clust.iter,])),pheno.factor),"Var1")

		gene.data <- melt(t(tissue.counts[cluster.ids==clust.iter,]))
		
		gene.data$Var1 <- as.character(gene.data$Var1)
				
		gene.data <-Sortby( merge(gene.data,pheno.data,by="Var1"),"Var1")
		
		if(weight){
		
			weight.info <- subset(match.info,V1==clust.iter)
		
			weight.info$w <- weight.info$V3
			weight.info$Var2 <- weight.info$V2
			
			gene.data <- merge(gene.data,weight.info,by="Var2")
	
			
			}#END: weight
		
		

		gene.data$round.value <- floor(gene.data$value)

		zero.count.data <- subset(aggregate(subset(gene.data,select="round.value"),subset(gene.data,select="Var2"),sum,na.rm=TRUE),round.value==0)
		
		print(clust.iter)
		
		gene.data <- subset(gene.data,!(Var2 %in% zero.count.data$Var2))
	  
	  if((mean(gene.data$value==0)<1)&(length(unique(gene.data$Var2))>1)){
	
		print(gene.data)
	
		try({	
	
		
		if(weight){
			gene.data$w <- gene.data$w*sum(!is.na(gene.data$value))/sum(gene.data$w)
			}else{
			gene.data$w <- 1}

	
	
		gene.data$Var2 <- as.factor(as.character(gene.data$Var2))
		
		gene.data <<- aggregate(subset(gene.data,select=c("round.value","value")),subset(gene.data,select=c("Var1","Var2","pheno.factor")),mean,na.rm=TRUE)

	
		#gee.out <- geeglm(value~pheno.factor+Var2,data=gene.data,id=Var1,family=poisson)
		
		#gee.out <- gee(value~pheno.factor+Var2,data=gene.data,id=Var1,family=poisson)
		
		
		#gene.data <- data.frame(Var2=rep(c("A","B"),each=10),Var1=as.factor(rbinom(20,1,0.5)),round.value=rpois(20,2),weights=1,pheno.factor=as.factor(sample(c(1,2),20,replace=TRUE)))
		

		if(identical(family,"GEE")){
			gene.data$Var0 <- as.numeric(as.factor(gene.data$Var1))
			gene.data <- Sortby(gene.data,"Var0")
			gee.out <- geeglm(value~pheno.factor+Var2,data=gene.data,id=Var0,family=poisson("log"))
			gee.sum <- coef(gee.out)["pheno.factor2"]
			pval <- summary(gee.out)$coef["pheno.factor2","Pr(>|W|)"]
			}
		
		
		
		if(identical(family,poisson)){
			gee.out <- glmer(round.value~pheno.factor+Var2+(1|Var1),data=gene.data,weights=gene.data$w,family=family)
			gee.sum <- lme4::fixef(gee.out)["pheno.factor2"]
			pval <- Anova(gee.out)["pheno.factor","Pr(>Chisq)"]
			}

		if(identical(family,quasipoisson)){
			print("quasipoisson")
			gee.out <- glmmPQL(round.value~pheno.factor+Var2,random=~1|Var1,family=family,data=gene.data,weights=gene.data$w)
			gee.sum <- fixef(gee.out)["pheno.factor2"]
			pval <- Anova(gee.out)["pheno.factor","Pr(>Chisq)"]
		
			}

		if(identical(family,"nbinom1")){
			gene.data$Var1 <- as.factor(gene.data$Var1)
			gene.data$Var2 <- as.factor(gene.data$Var2)
			gee.out <- glmmadmb(round.value~pheno.factor+Var2+(1|Var1),data=gene.data,family=family,zeroInflation=TRUE)
			gee.sum <- coef(gee.out)["pheno.factor2"]
			pval <- Anova(gee.out)["pheno.factor","Pr(>Chisq)"]
			}


#		if(identical(family,"MGLM")){
#			temp <- dcast(gene.data,Var1 + pheno.factor ~ Var2,value.var="value")
#			y <- as.matrix(Drop.variable(temp,c("Var1","pheno.factor")))
#			x <- as.matrix(subset(temp,select="pheno.factor"))
#			gee.out <- MGLMreg(y ~ x,dist="NegMN")
#		    gee.sum <- mean(coef(gee.out)$alpha["x2",])
#			pval <- gee.out$test["x2","Pr(>wald)"]
#			}


		#results.fc[clust.iter,] <- 1

		print(clust.iter)

		gee.list[[clust.iter]] <- gee.out

		#results.fc[clust.iter,] <- c(exp(-gee.sum["pheno.factor2","Estimate"]),gee.sum["pheno.factor2","Pr(>|W|)"])
		
		#results.fc[clust.iter,] <- c(exp(-gee.sum["pheno.factor2","Estimate"]),2*pnorm(-abs(gee.sum["pheno.factor2","Robust z"])))
		
		results.fc[clust.iter,] <- c(exp(-gee.sum),pval,0)#gee.out@optinfo$conv

		
		})
		
		}
		
		#results.fc <- results.fc[cluster.row.fits,]
		
		base.time <- timer.out(base.time,step=1,itermax=length(cluster.row.fits)*ifelse(perm==0,1,perm))

	} # loop over clust.iter
	
	if(perm>0){perm.results <- rbind(perm.results,data.frame(V1=rownames(results.fc),results.fc))}
	
	} # loop over perm

	hits <- rownames(results.fc[!is.na(results.fc[,1]),])
	
	if(perm==0){
		out.data <- data.frame(V1=rownames(results.fc),results.fc)
		}else{
			out.data <- perm.results
	}

	return(list(out.data=out.data,gee.list=gee.list))

} # END: fit.gee.5.array


timer.out <-  function(iter.list=NULL,itermax=NA,step=NA,percent=NA){
	
# prints the iteration, % remaining, Estimated elapsed time to completion
# iter = iteration
# itermax = iterations until complete
# step = the # of steps to skip until output
# percent = the % to skip until output
# start = TRUE/FALSE to reset timer
# new.time.print.skip timing variable
	
	elapsed.element <- 3
	
	if(is.null(iter.list)){
		iter <- 1
		base.time <- proc.time()
		return(list(time=base.time,iter=iter))
		
	}else{
		iter <- iter.list$iter + 1
		base.time <- iter.list$time
	}
	
	proportion.completed <- iter/itermax
	

	new.time <- proc.time()
	
	elapsed.time <- new.time[elapsed.element] - base.time[elapsed.element]
	
	time.to.completion <- elapsed.time/proportion.completed - elapsed.time
	
	
	elapsed.time <- round(elapsed.time,1)
	time.to.completion <- round(time.to.completion,1)
	
	if(!(is.na(step))){
		if( !(iter %% step)){
			print(paste("Iter = ",iter,"&",floor(100*iter/itermax),"% ","Elapsed ",elapsed.time,"(secs) Time to Completion",time.to.completion,"(secs)"))
		}
	}else{
		if(!(is.na(percent))){
			pct.completes <- floor(100*c(iter-1,iter)/itermax)
			steps <- floor(pct.completes/percent)
			if( steps[1]!=steps[2]){
				print(paste(pct.completes[2],"% Complete","Elapsed Time",elapsed.time,"(secs) Time to Completion",time.to.completion,"(secs)"))
			}
		}
	}
	
	
	
	return(list(time=base.time,iter=iter))
	
}

