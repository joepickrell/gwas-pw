
computeprior = function(z,pi){
	dvec = tabulate(z+1,nbin=3)
	d = length(z)
	return(ifelse(dvec[2]>0,(1-pi)*(1/d) * (1/(dvec[2]+dvec[3])) *
	1/choose(d, dvec[1]) * 1/choose(d-dvec[1],dvec[2]), 0))
}

#picks out the partition gamma with all 1s (ie all in D)
allones = function(gamma){return(prod(gamma==1))}
#note the "drop=FALSE" commands below stop R collapsing matrices into vectors inappropriately
#VYX \approx (1/n) Y'X is d by p
#VYY \approx (1/n) Y'Y is d by d
#VXX is a p-vector of the estimated variances of the SNP
logBF.fromVSummaries = function(VYX,VYY,VXX,U,D,n,m,d,sigmaa){
	dd = sum(D)
	du= sum(U)
	p = dim(VYX)[2]
	if(du>0){
		LUU = chol(VYY[U,U,drop=FALSE]) # a du by du matrix
		VUD = VYY[U,D,drop=FALSE] #a du by dd matrix of the correlations of Yu with Yd
		c = cbind(forwardsolve(t(LUU),VYX[U,,drop=FALSE]))#c solves LUU'c = phiU, c is a du by p matrix
		b = cbind(forwardsolve(t(LUU), VUD)) # b is du by dd, and solves LUU' b = VUD,
		#so b'b = VUD' LUU^-1 LUU'^-1 VUD = VUD' (LUU'LUU)^-1 VUD = VUD'VYYU^-1 VUD
	} else{c=matrix(0,nrow=1,ncol=p); b=matrix(0,nrow=1,ncol=dd);}
	C = VXX - colSums(c*c)
	u = VYX[D,,drop=FALSE] - crossprod(b,c)
	V0 = VYY[D,D,drop=FALSE] - crossprod(b)
	L0 = chol(V0)
	a = forwardsolve(t(L0),u)
	lambda = sigmaa^(-2) / (n*C)
	k = as.numeric(1/(1+lambda))
	return((dd/2) * log(1-k) - 0.5*(n+m-(d-sum(D)-sum(U)))*log(1-(k/C) *colSums(a*a)))
}

#G is an n by p matrix of SNP genotypes
#Y is an n by d matrix of phenotypes
#single-SNP BFs are computed for each SNP (column of G)
logBF.rankone.matrix = function(G,Y,sigmaa,pi0=0.5,m=0){
	if(is.null(dim(G))){G=cbind(G)} # turn a vector of genotypes into a matrix
	subset = complete.cases(Y) & complete.cases(G)
	Y=Y[subset,,drop=FALSE]
	G=G[subset,,drop=FALSE]
	n = dim(Y)[1]
	d = dim(Y)[2]
	p = dim(G)[2] #number of SNPs
	if(m==0){m = d-1}
	Y =scale(Y,center=T,scale=F) #center Y and G to avoid needing intercept in regression
	G = scale(G,center=T,scale=F)
	VYX = (1/n)*crossprod(Y,G) # this is (1/n) t(Y) %*% G, a d by p matrix
	VYY = (1/n)*crossprod(Y) # (1/n) t(Y) %*% Y, a d by d matrix
	VXX = (1/n)*colSums(G*G) # a p vector of (1/n) ||g|| values
	print(VYX)
	print(VYY)
	print(VXX)
	prior = rep(0,3^d)
	gamma=matrix(0,nrow=3^d,ncol=d)
	lbf = matrix(0,nrow=3^d, ncol=p)
	for(i in 0:(3^d-1)){
		for(j in 1:d){
			gamma[i+1,j]= (i %% 3^j) %/% 3^{j-1}
		}
		prior[i+1] = computeprior(gamma[i+1,],pi0)
		U = (gamma[i+1,]==0)
		D = (gamma[i+1,]==1)
		if(prior[i+1]>0){
			BF = 0
			for(ss in 1:length(sigmaa)){
				BF = BF+exp(logBF.fromVSummaries(VYX,VYY,VXX,U,D,n,m,d,sigmaa[ss]))
				#note we just don't bother computing for models with prior = 0
			}
			lbf[i+1,]=log(BF/length(sigmaa))
		} else {lbf[i+1,] = 0}
	}
	prior[1] = pi0
	BF=exp(lbf)
	posterior = prior[-1]*BF[-1,,drop=FALSE]
	normalize=function(x){return(x/sum(x))}
	posterior = apply(posterior,2,normalize)
	p0=t(gamma[-1,]==0) %*% posterior
	p1=t(gamma[-1,]==1) %*% posterior
	p2=t(gamma[-1,]==2) %*% posterior
	#divide by log(10) to convert everything to log base 10
	lbfav = log(colSums(prior[-1]*exp(lbf[-1,,drop=FALSE]))/sum(prior[-1]))/log(10)
	lbfuni.comps = apply(lbf[rowSums(gamma)==(2*d-1),,drop=FALSE],2,rev)/log(10)
	lbfuni = log10(apply(10^lbfuni.comps,2,mean))
	lbfall = lbf[which.max(apply(gamma,1,allones)),,drop=FALSE]/log(10)
	lbf = lbf/log(10)
	return(list(prior=prior,gamma=gamma,lbf=lbf,lbfav = lbfav,lbfuni=lbfuni,
	lbfall=lbfall,lbfuni.comps=lbfuni.comps,p0=p0,p1=p1,p2=p2))
}

