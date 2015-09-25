


# cluster.homology.de fits the cluster homology model for 2 conditions specified by pheno.factor
# tissue.counts.norm is the RNA seq count matrix with named columns
# match.info contains the homology cluster mapping 


cluster.homology.de <- function(tissue.counts.norm,match.info,family.fit="GEE",pheno.factor,permnumber=4,weighted=TRUE,k=20,ncores=1){

	#tissue.counts.norm = matrix with rows are transcripts and columns are samples
	#match.info = cluster identities columns are cluster=[cluster identitiy], transcript=[target transcript matches rows in tissue.counts], weight=[homology weights]
	#pheno.factor = factor representing differential expression analysis grouping (e.g., treated or controls)
	#weighted = TRUE?FALSE whether to perform homology weighted analysis
	#family.fit = "GEE", quasipoisson what type of cluster model to fit
	#permnumber = how many permutations of full transcriptome (nrows(tissue.counts)*permnumber)

	# k is integer that divides the total compute work to reduce memory use for parallel computations	
	# ncores if > 1 then will use parallel processing on multicore CPUs


# Required R pacakges
	
	require(DESeq2)
	require(plyr)
	require(dplyr)
	require(geepack)
	require(MASS)
	require(reshape2)
	require(doMC)
	
	
	familyiter <- family.fit


	analysis.out<- list()
	analysis.out.perm <-list()

	clusters <- unique(match.info$cluster)

	k <- ifelse(ncores==1,1,k)

	kdf <- data.frame(cluster=clusters,subset=sample(1:k,size=length(clusters),replace=TRUE))


	match.info$V1 <- match.info$cluster
	match.info$V2 <- match.info$transcript

# Breaks total compute work into subsets to use less memory
	for(p in 1:k){
		pchar <- as.character(p)
		system.time({

		setter <- as.character(match.info$cluster) %in% as.character(subset(kdf,subset==p)$cluster)
	
		match.info2 <- subset(match.info,setter )	
		tissue.counts <- tissue.counts.norm[setter, ]	
	
		analysis.out[[pchar]] <- fit.gee.parallel.array.test(tissue.counts,match.info2,pheno.factor,crash.list=NULL,weight=weighted,perm=0,family=familyiter,cores=ncores)

		analysis.out.perm[[pchar]] <- fit.gee.parallel.array.test(tissue.counts,match.info2,pheno.factor,crash.list=NULL,
										weight=weighted,perm=permnumber,family=familyiter,cores=ncores)

		})

		}


# Combine permuted and orignal into 1 dataframe

wtddf.perm <- ldply(analysis.out.perm,function(x){
	
	df1 <- data.frame(x[[1]],model=familyiter)
	
	return(df1)

})

wtddf <- ldply(analysis.out,function(x){
	
	df1 <- data.frame(x[[1]],model=familyiter)
	
	return(df1)

})

# Compute permutation p-value

wtd.results <- ddply(wtddf,"model",function(x){
	
	
	model1 <- x$model[1]
	
	permed <- subset(wtddf.perm,model==model1)
	
	x$pvalue.perm <- ecdf(permed$pvalue)(x$pvalue)
	
	return(x)
	
})


#Use DESeq to compute singleton differential expression heterogeneity report

	dds <- DESeqDataSetFromMatrix(round(tissue.counts.norm),data.frame(condition=pheno.factor),design=~condition)

	dds <- DESeq(dds)

	results <- data.frame(results(dds))
	results$V2 <- rownames(results)
	results2 <- merge(match.info,results)

heterogeneity.report <- ddply(results2,c("V1"),function(x){
	
	alpha <- 0.05
	
	bonferroni <- alpha/nrow(x)
	
	upregulated <- sum((x$log2FoldChange>0) & (x$pvalue<bonferroni),na.rm=TRUE)
	downregulated <- sum((x$log2FoldChange<0) & (x$pvalue<bonferroni),na.rm=TRUE)
	
	heterogeneity <- (upregulated>0) & (downregulated>0)
	
	data.frame(total.transcripts=nrow(x),transcripts.upregulated=upregulated,transcripts.downregulated=downregulated,heterogeneity)	
	
})


heterogeneity.report$cluster <- heterogeneity.report$V1
heterogeneity.report$V1 <- NULL

resultsWithHeterogeneity <- merge(wtd.results, heterogeneity.report)

dim(resultsWithHeterogeneity)

resultsWithHeterogeneity$V1 <- NULL
resultsWithHeterogeneity$.id <- NULL
resultsWithHeterogeneity$conv<- NULL

resultsWithHeterogeneity$FoldChange <- resultsWithHeterogeneity$FC

#Output summary by cluter

results <- subset(resultsWithHeterogeneity,select=c("cluster","model","FoldChange", "pvalue","pvalue.perm", "total.transcripts",
						 "transcripts.upregulated", "transcripts.downregulated", "heterogeneity"))
						 
						 
return(results)						 

}






fit.gee.parallel.array.test <- function(tissue.counts,match.info,pheno.factor,crash.list=NULL,weight=FALSE,perm=0,family=poisson,cores=8){

	require(geepack)
	require(MASS)
	require(reshape2)
	#weights 
	# better 0 handling
	
	#match.info <- match.info2
	#crash.list=1:1000
	#perm <- 0
	#perm.iter <- 1
	#weight=TRUE
	#family="GEE"
	#family=quasipoisson
	#cores <- 1
	require(doMC)
	require(MASS)

	gee.out <- NULL

	if(cores>1)	{registerDoMC(cores)}
	
	cluster.ids <- match.info$V1

	n.clusters <- length(unique(cluster.ids))

	cluster.names <- unique(cluster.ids)
	clust.iter <- cluster.names[1]
	
	if(!identical(rownames(tissue.counts),match.info$V2)){stop("Cluster ids don't match expression matrix")}
	

	if(!weight){match.info$weight <- 1}

	all.data <- data.frame(cluster=match.info$V1,gene=rownames(tissue.counts),tissue.counts,weight=match.info$weight)
		
	pheno.factor.0 <- pheno.factor
	
	gee.list <- list()
	

	cluster.order <- sort(all.data$cluster)

	x <- subset(all.data,cluster==all.data$cluster[1])
	system.time({
		
	fitout <- dlply(all.data,"cluster",function(x){

		
		results.fc.list <- list()
		gee.out.list <- list()
		
		#print(which(cluster.order==x$cluster[1])/length(cluster.order))
		
		for(perm.iter in 1:ifelse(perm==0,1,perm)){

			if(perm>0){pheno.factor <- pheno.factor.0[sample(length(pheno.factor.0))]}
		
			pheno.data <- data.frame(Var1=colnames(tissue.counts), pheno.factor,stringsAsFactors=0)

			gene.data <- melt(x,id.vars=c("cluster","gene","weight"))
		
			gene.data$Var1 <- as.character(gene.data$variable)
				
			gene.data <-Sortby( merge(gene.data,pheno.data,by="Var1"),"Var1")
		

			gene.data$round.value <- floor(gene.data$value)

			zero.count.data <- subset(aggregate(subset(gene.data,select="round.value"),subset(gene.data,select="gene"),sum,na.rm=TRUE),round.value==0)
		
			gene.data <- subset(gene.data,!(gene %in% zero.count.data$gene))
	  		
	  		
	  
	  		#print(gene.data)
	  
	  		if((mean(gene.data$value==0)<1)&(length(unique(gene.data$gene))>1)){
	
			# try model fit
			
			try({	
	
		
			if(weight){
				gene.data$w <- gene.data$weight*sum(!is.na(gene.data$value))/sum(gene.data$w)
			}else{
			gene.data$w <- 1}

			gene.data$gene <- as.factor(as.character(gene.data$gene))
		
			gene.data <<- aggregate(subset(gene.data,select=c("round.value","value")),subset(gene.data,select=c("Var1","gene","pheno.factor")),mean,na.rm=TRUE)

		gee.out <- NULL

		if(identical(family,"GEE")){
			gene.data$Var0 <- as.numeric(as.factor(gene.data$Var1))
			gene.data <- Sortby(gene.data,"Var0")
			#print("into GEE")
			gee.out <- geeglm(value~pheno.factor+gene,data=gene.data,weights=gene.data$w,id=Var0,family=poisson("log"))
			gee.sum <- coef(gee.out)["pheno.factor2"]
			pval <- summary(gee.out)$coef["pheno.factor2","Pr(>|W|)"]
			}
		
		if(identical(family,poisson)){
			gee.out <- glmer(round.value~pheno.factor+gene+(1|Var1),data=gene.data,weights=gene.data$w,family=poisson)
			gee.sum <- lme4::fixef(gee.out)["pheno.factor2"]
			pval <- Anova(gee.out)["pheno.factor","Pr(>Chisq)"]
			}


		if(identical(family,quasipoisson)){
			#print("quasipoisson")
			gee.out <- glmmPQL(round.value~pheno.factor+gene,random=~1|Var1,family=quasipoisson,data=gene.data,weights=gene.data$w,verbose=FALSE)
			
			gee.sum <-  gee.out$coef$fixed["pheno.factor2"]
			
		#	print("gee.sum")
		#	print(gee.out)
			
		#	print(Anova(gee.out))
			
			pval <- Anova(gee.out)["pheno.factor","Pr(>Chisq)"]
		#	print("quasidone")
			}

			gee.out.list[[perm.iter]] <- gee.out
			results.fc.list[[perm.iter]] <- data.frame(V1=x$cluster[1],FC=exp(-gee.sum),pvalue=pval,conv=0)

		
		})# end try model fit
		}else{gee.out <- NULL} # if something nonzero
		} #loop over perm
	

	if(perm==0){
		gee.out.list <- gee.out
		if(length(results.fc.list)>0){results.fc.list <- results.fc.list[[perm.iter]]}
		}else{
			results.fc.list <- rbind.fill(results.fc.list)
			}

	return(list(gee.out=gee.out.list,results.fc=results.fc.list))
	},.parallel=(cores>1),.progress="text")
	})	

	#print("finish ddply")


out.data <- ldply(fitout,function(x){
	return(data.frame(x$results.fc))
})

	return(list(out.data=out.data,gee.list=fitout$gee.out.list))

} # END: fit.gee.parallel.array.test 


Sortby <- function(df,sorting.variable,...){

# sort data.frame df by variable named "sorting variable"	
	
return(df[order(df[[sorting.variable]],...),])
	
}



