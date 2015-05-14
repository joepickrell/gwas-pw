# approximate BFs using multivarite Wakefield approximation
library(mvtnorm)
approxBF = function(g, p, cc = 0){
	lm1 = summary(lm(p[,1] ~g))$coef
	lm2 = summary(lm(p[,2]~ g))$coef
	Z1 = lm1[2,1]/lm1[2,2]
	Z2 = lm2[2,1]/lm2[2,2]
	cc = cor(p[,1], p[,2])
	V1 = lm1[2,2]^2
	V2 = lm2[2,2]^2
	abf1 = ABF1(Z1, Z2, V1, V2, 0.5, cc)
	abf2 = ABF2(Z1, Z2, V1, V2, 0.5, cc)
	abf3 = ABF3(Z1, Z2, V1, V2, 0.5, cc)
	WABF1 = WABF(Z1, V1, 0.5)
	WABF2 = WABF(Z2, V2, 0.5)
	return(c(abf1, abf2, abf3, Z1, Z2, V1, V2, WABF1, WABF2))	
}

WABF = function(Z, V, W){
        r = W/(V+W)
        toreturn = log(sqrt(1-r))
        tmp = Z*Z*r
        toreturn = toreturn + tmp/ 2
        return(toreturn)
}

ABF1 = function(Z1, Z2, V1, V2, W, C){
	r = W/(V1+W)
	toreturn = log(sqrt(1-r))
	tmp = Z1*Z1*r- 2*C*Z1*Z2*(1-sqrt(1-r))
	toreturn = toreturn + tmp/ (2*(1-C*C))
	return(toreturn)
}

ABF2 = function(Z1, Z2, V1, V2, W, C){

	r = W/ (V2+W)
	toreturn = log ( sqrt(1-r) )

	tmp = Z2*Z2*r- 2*C*Z1*Z2*(1-sqrt(1-r))
	toreturn = toreturn+ tmp/ (2*(1-C*C))
	return(toreturn)
}

ABF3 = function(Z1, Z2, V1, V2, W, C){

	r1 = W/ (V1+W)
	r2 = W/ (V2+W)
	toreturn = log ( sqrt(1-r1) ) + log(sqrt(1-r2));

	tmp = Z1*Z1*r1+Z2*Z2*r2- 2*C*Z1*Z2*(1-sqrt(1-r1)*sqrt(1-r2))
	toreturn = toreturn+ tmp/ (2*(1-C*C))

	return(toreturn)
}

ABF3_sep = function(Z1, Z2, V1, V2, W, C){

	B1 = Z1*sqrt(V1)
	B2 = Z2*sqrt(V2)
	C1 = matrix(nrow = 2, ncol = 2)
	C1[1,1] = V1+W
	#C1[1,2] = C*W
	#C1[2,1] =  C*W
	C1[1,2] = C*sqrt( (V1+W)*(V2+W))
	C1[2,1] =  C*sqrt( (V1+W)*(V2+W))
	C1[2,2] = V2+W
	d1 = dmvnorm(c(B1, B2), mean = c(0, 0), sigma = C1, log = T)
	
	C1[1,2] = 0
	C1[2,1] = 0
	C1[1,1] = V1
	C1[2,2] = V2
	d2 = dmvnorm(c(B1, B2), mean = c(0, 0), sigma = C1, log = T)
	return(d1-d2)
        #r1 = W/ (V1+W)
        #r2 = W/ (V2+W)
        #toreturn = log ( sqrt(1-r1) ) + log(sqrt(1-r2));

        #tmp = Z1*Z1*r1+Z2*Z2*r2- 2*C*Z1*Z2*(1-sqrt(1-r1)*sqrt(1-r2))
        #toreturn = toreturn+ tmp/ (2*(1-C*C))

        #return(toreturn)
}

llk_sep = function(Z1, Z2, V1, V2, W, C){

        B1 = Z1*sqrt(V1)
        B2 = Z2*sqrt(V2)
        C1 = matrix(nrow = 2, ncol = 2)
        C1[1,1] = V1+W
        #C1[1,2] = C*W
        #C1[2,1] =  C*W
        C1[1,2] = C*sqrt( (V1+W)*(V2+W))
        C1[2,1] =  C*sqrt( (V1+W)*(V2+W))
        C1[2,2] = V2+W
        print(C1)
	print(paste(B1, B2))
	d1 = dmvnorm(c(B1, B2), mean = c(0, 0), sigma = C1, log = T)
        return(d1)
        #r1 = W/ (V1+W)
        #r2 = W/ (V2+W)
        #toreturn = log ( sqrt(1-r1) ) + log(sqrt(1-r2));

        #tmp = Z1*Z1*r1+Z2*Z2*r2- 2*C*Z1*Z2*(1-sqrt(1-r1)*sqrt(1-r2))
        #toreturn = toreturn+ tmp/ (2*(1-C*C))

        #return(toreturn)
}
