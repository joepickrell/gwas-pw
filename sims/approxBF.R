

approxBF = function(g, p){
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
	return(c(abf1, abf2, abf3))	
}

ABF1 = function(Z1, Z2, V1, V2, W, C){
	print(paste(Z1, Z2, V1, V2, W, C))
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
