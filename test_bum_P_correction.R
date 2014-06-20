

library(ClassComparison)
source("http://bioinformatics.mdanderson.org/OOMPA/oompaLite.R")
oompaLite()

n <- 10000

p.skewed <- 1-pchisq(rchisq(n,20),19)

hist(p.skewed)

bum.perm <- Bum(p.skewed)
par(mfrow=c(1,1))

hist(bum.perm,main="Permuted P-values as Beta-Uniform Mixture")

alpha.hat <- bum.perm@ahat
l.hat <- bum.perm@lhat
pi0.hat <- bum.perm@pihat

abline(h=pi0.hat,lty=2,col=2)

x <- seq(0,1,length=1000)

y <- l.hat + (1-l.hat)*alpha.hat *x^(alpha.hat-1)

empirical.p.function <- function(x,l.hat,a.hat){
	
	return(l.hat*x+(1-l.hat)*x^a.hat)
	
	}

hist(empirical.p.function(p.skewed,l.hat,alpha.hat),main="Corrected Permuted P-Values")

bum.perm.cor <- Bum(empirical.p.function(p.skewed,l.hat,alpha.hat))

hist(bum.perm.cor,main="Corrected Permuted P-values as Beta-Uniform Mixture")
