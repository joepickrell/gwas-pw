library(MethComp)


addbetas = function(d){
	d$B1 = d$Z_1*sqrt(d$V_1)
	d$B2 = d$Z_2*sqrt(d$V_2)
	return(d)
}

rankmodels = function(d){
	l12 = llk_1_causes_2(d)
	l21 = llk_2_causes_1(d)
	lindep = llk_indep_effects(d)
	lsame = llk_samepheno(d)
	toreturn = data.frame(matrix(nrow = 4, ncol= 3))
	toreturn[1,1] = "1->2"
	toreturn[2,1] = "2->1"
	toreturn[3,1] = "indep"
	toreturn[4,1] = "same"
	toreturn[1,2] = l12$llk
	toreturn[2,2] = l21$llk
	toreturn[3,2] = lindep$llk
	toreturn[4,2] = lsame$llk
	toreturn[1,3] = l12$beta
	toreturn[2,3] = l21$beta
	toreturn[3,3] = lindep$beta
	toreturn[4,3] = lsame$beta
	toreturn = toreturn[order(-toreturn[,2]),]
	return(toreturn)

}

llk_1_causes_2 = function(d){
	d = addbetas(d)
	tmp = d[d$MODEL == 3 | d$MODEL == 4 | d$MODEL == 1,]
	mv1 = mean(d$V_1)
	mv2 = mean(d$V_2)
	ratio = mv2/mv1
	dm = Deming(tmp$B1, tmp$B2, ratio)
	alpha = dm[1]
	beta = dm[2]
	d$EXPB1 = (ratio * d$B1 + beta *(d$B2 - alpha))/ (ratio + beta*beta)
	d$EXPB2 = alpha + beta * d$EXPB1
	d[d$MODEL == 2 | d$MODEL == "4_2",]$EXPB1 = 0
	d[d$MODEL == 2 | d$MODEL == "4_2",]$EXPB2 = d[d$MODEL == 2|  d$MODEL == "4_2",]$B2
	d$LLK1 = apply(d[c("B1", "EXPB1")], 1, FUN = function(x){ return(dnorm(x[1], x[2], sd = sqrt(mv1), log = T))})
	d$LLK2 = apply(d[c("B2", "EXPB2")], 1, FUN = function(x){ return(dnorm(x[1], x[2], sd = sqrt(mv2), log = T))})
	s1 = sum(d$LLK1)
	s2 = sum(d$LLK2)
	llk = s1+s2
	return(list("llk" = llk, d = d, alpha = alpha, beta = beta))
}


llk_2_causes_1 = function(d){
	d = addbetas(d)
	tmp = d[d$MODEL == 3 | d$MODEL == "4_2" | d$MODEL == 2,]
	mv1 = mean(d$V_1)
	mv2 = mean(d$V_2)
	ratio = mv2/mv1
	dm = Deming(tmp$B1, tmp$B2, ratio)
	alpha = dm[1]
	beta = dm[2]
	d$EXPB1 = (ratio * d$B1 + beta *(d$B2 - alpha))/ (ratio + beta*beta)
	d$EXPB2 = alpha + beta * d$EXPB1
	d[d$MODEL == 1 | d$MODEL == "4",]$EXPB2 = 0
	d[d$MODEL == 1 | d$MODEL == "4",]$EXPB1 = d[d$MODEL == 1|  d$MODEL == "4",]$B1
	d$LLK1 = apply(d[c("B1", "EXPB1")], 1, FUN = function(x){ return(dnorm(x[1], x[2], sd = sqrt(mv1), log = T))})
	d$LLK2 = apply(d[c("B2", "EXPB2")], 1, FUN = function(x){ return(dnorm(x[1], x[2], sd = sqrt(mv2), log = T))})
	s1 = sum(d$LLK1)
	s2 = sum(d$LLK2)
	llk = s1+s2
	return(list("llk" = llk, d = d, alpha = alpha, beta = beta))
}


llk_indep_effects = function(d){
	d = addbetas(d)
	tmp = d[d$MODEL == 3,]
	mv1 = mean(d$V_1)
	mv2 = mean(d$V_2)
	ratio = mv2/mv1
	dm = Deming(tmp$B1, tmp$B2, ratio)
	alpha = dm[1]
	beta = dm[2]
	d$EXPB1 = (ratio * d$B1 + beta *(d$B2 - alpha))/ (ratio + beta*beta)
	d$EXPB2 = alpha + beta * d$EXPB1
	d[d$MODEL == 1 | d$MODEL == "4",]$EXPB2 = 0
	d[d$MODEL == 2 | d$MODEL == "4_2",]$EXPB1 = 0
	d[d$MODEL == 1 | d$MODEL == "4",]$EXPB1 = d[d$MODEL == 1|  d$MODEL == "4",]$B1
	d[d$MODEL == 2 | d$MODEL == "4_2",]$EXPB2 = d[d$MODEL == 2|  d$MODEL == "4_2",]$B2
	d$LLK1 = apply(d[c("B1", "EXPB1")], 1, FUN = function(x){ return(dnorm(x[1], x[2], sd = sqrt(mv1), log = T))})
	d$LLK2 = apply(d[c("B2", "EXPB2")], 1, FUN = function(x){ return(dnorm(x[1], x[2], sd = sqrt(mv2), log = T))})
	s1 = sum(d$LLK1)
	s2 = sum(d$LLK2)
	llk = s1+s2
	return(list("llk" = llk, d = d, alpha = alpha, beta = beta))
}


llk_samepheno = function(d){
	d = addbetas(d)
	mv1 = mean(d$V_1)
	mv2 = mean(d$V_2)
	ratio = mv2/mv1
	dm = Deming(d$B1, d$B2, ratio)
	alpha = dm[1]
	beta = dm[2]
	d$EXPB1 = (ratio * d$B1 + beta *(d$B2 - alpha))/ (ratio + beta*beta)
	d$EXPB2 = alpha + beta * d$EXPB1
	d$LLK1 = apply(d[c("B1", "EXPB1")], 1, FUN = function(x){ return(dnorm(x[1], x[2], sd = sqrt(mv1), log = T))})
	d$LLK2 = apply(d[c("B2", "EXPB2")], 1, FUN = function(x){ return(dnorm(x[1], x[2], sd = sqrt(mv2), log = T))})
	s1 = sum(d$LLK1)
	s2 = sum(d$LLK2)
	llk = s1+s2
	return(list("llk" = llk, d = d, alpha = alpha, beta = beta))
}

