
library(MASS)


simgp = function(NIND, betas, f, S){
	g = rbinom(NIND, 2, f)
	N0 = sum(g == 0)
	N1 = sum(g == 1)
	N2 = sum(g == 2)
	p0 = matrix(nrow = N0, ncol = 2)
	p1 = matrix(nrow = N1, ncol = 2)
	p2 = matrix(nrow = N2, ncol = 2)
	if (N0 >0){ p0 = mvrnorm(N0, mu = c(0,0), Sigma = S)}
	if (N1 >0){ p1 = mvrnorm(N1, mu = c(betas[1],betas[2]), Sigma = S)}
	if (N2 > 0){ p2 = mvrnorm(N2, mu = c(2*betas[1],2*betas[2]), Sigma = S)}
	allp = matrix(nrow = length(g), ncol = 2)
	allp[which(g == 0),] = p0
	allp[which(g == 1),] = p1
	allp[which(g == 2),] = p2
	return(list("pheno" = allp, "geno" = g))
}


#for(i in 1:NSIM){
#	p = mvrnorm(NTOTAL, mu = c(0, 0), Sigma = S)
#	p[1:N1,1] = NA
#	p[(NTOTAL-N2):NTOTAL,2] = NA
#	
#	f = 0.3
#	g = rbinom(NTOTAL, 2, f)
#	lm1 = summary(lm(p[,1] ~g))$coef
#	lm2 = summary(lm(p[,2]~ g))$coef
#	Z1 = lm1[2,1]/lm1[2,2]
#	Z2 = lm2[2,1]/lm2[2,2]
#	toplot[i,1] = Z1
#	toplot[i,2] = Z2
#	toplot[i,3] = lm1[2,1]
#	toplot[i,4] = lm2[2,1]
#}

