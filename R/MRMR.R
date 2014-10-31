source("approxBF.R")

sumlog = function(logx, logy){
       if (logx > logy) {
		return(logx + log(1 + exp(logy-logx)))
	}
        else {
		return (logy + log(1 + exp(logx-logy)))
	}
}


addBFs = function(d, C = 0, W){
	# assume data has names Z_1 V_1 Z_2 V_2 for Z-scores and variances of each
	# C is the prior covariance of the effect sizes on the 2 phenotypes
	# W is a vector of prior variance (average over all)
	wz1 = which(names(d) == "Z_1")
	wz2 = which(names(d) == "Z_2")
	wv1 = which(names(d) == "V_1")
	wv2 = which(names(d) == "V_2")
	d$BF1 = apply(d[,c(wz1, wz2, wv1, wv2)], 1, FUN = function(x){ return( avgBF1(x[1], x[2], x[3], x[4], C, W))})
	d$BF2 = apply(d[,c(wz1, wz2, wv1, wv2)], 1, FUN = function(x){ return( avgBF2(x[1], x[2], x[3], x[4], C, W))})
	d$BF3 = apply(d[,c(wz1, wz2, wv1, wv2)], 1, FUN = function(x){ return( avgBF3(x[1], x[2], x[3], x[4], C, W))})
	return(d)
		
	

}

llk_no_overlap = function(theta, d){
	alpha_1 = theta[1]
	pi_1 = exp(alpha_1)/( exp(alpha_1)+ 1)
	pi_2 = 1-pi_1

	wbf1 = which(names(d) == "BF1")
	wbf2 = which(names(d) == "BF2")
	sllk = apply( d[,c(wbf1, wbf2)], 1, FUN = function(x){ return(sumlog(log(pi_1)+x[1], log(pi_2)+x[2]))})
	return(sum(sllk))

}


llk_1_causes_2 = function(theta, d){
        alpha_2 = theta[1]
	lambda = theta[2]
        pi_2 = exp(alpha_2)/( exp(alpha_2)+ 1)
        pi_3 = 1-pi_2

       	wz1 = which(names(d) == "Z_1")
        wz2 = which(names(d) == "Z_2")
        wv1 = which(names(d) == "V_1")
        wv2 = which(names(d) == "V_2")
	d$BF3 = apply(d[,c(wz1, wz2, wv1, wv2)], 1, FUN = function(x){ return( avgBF3(x[1], x[2], x[3], x[4], lambda, W))})
        wbf2 = which(names(d) == "BF2")
        wbf3 = which(names(d) == "BF3")
        sllk = apply( d[,c(wbf12, wbf3)], 1, FUN = function(x){ return(sumlog(log(pi_2)+x[1], log(pi_3)+x[2]))})
        return(sum(sllk))

}



llk_2_causes_1 = function(theta, d){
        alpha_1 = theta[1]
        lambda = theta[2]
        pi_1 = exp(alpha_1)/( exp(alpha_1)+ 1)
        pi_3 = 1-pi_1

        wz1 = which(names(d) == "Z_1")
        wz2 = which(names(d) == "Z_2")
        wv1 = which(names(d) == "V_1")
        wv2 = which(names(d) == "V_2")
        d$BF3 = apply(d[,c(wz1, wz2, wv1, wv2)], 1, FUN = function(x){ return( avgBF3(x[1], x[2], x[3], x[4], lambda, W))})
        wbf1 = which(names(d) == "BF1")
        wbf3 = which(names(d) == "BF3")
        sllk = apply( d[,c(wbf1, wbf3)], 1, FUN = function(x){ return(sumlog(log(pi_1)+x[1], log(pi_3)+x[2]))})
        return(sum(sllk))

}


llk_samepheno = function(theta, d){
        lambda = theta[1]

        wz1 = which(names(d) == "Z_1")
        wz2 = which(names(d) == "Z_2")
        wv1 = which(names(d) == "V_1")
        wv2 = which(names(d) == "V_2")
        d$BF3 = apply(d[,c(wz1, wz2, wv1, wv2)], 1, FUN = function(x){ return( avgBF3(x[1], x[2], x[3], x[4], lambda, W))})
        sllk = d$BF3
        return(sum(sllk))

}

llk_indep_effect = function(theta, d){
        alpha_1 = theta[1]
        alpha_2 = theta[2]
        pi_1 = exp(alpha_1)/( exp(alpha_1)+ exp(alpha_2)+1)
	pi_2 = exp(alpha_2)/( exp(alpha_1)+ exp(alpha_2)+1)
        pi_3 = 1-pi_1-pi_2

        wz1 = which(names(d) == "Z_1")
        wz2 = which(names(d) == "Z_2")
        wv1 = which(names(d) == "V_1")
        wv2 = which(names(d) == "V_2")
        d$BF3 = apply(d[,c(wz1, wz2, wv1, wv2)], 1, FUN = function(x){ return( avgBF3(x[1], x[2], x[3], x[4], 0, W))})
        wbf1 = which(names(d) == "BF1")
        wbf2 = which(names(d) == "BF2")
        wbf3 = which(names(d) == "BF3")
        sllk = apply( d[,c(wbf1, wbf2, wbf3)], 1, FUN = function(x){ return(sumlog(log(pi_1)+x[1], sumlog(log(pi_2)+x[2], log(pi_3)+x[2]))})
        return(sum(sllk))

}



avgBF1 = function(Z1, Z2, V1, V2, C = 0, W){
	toreturn = ABF1(Z1, Z2, V1, V2, W[1], C)
	
	if (length(W)>1){
		for (i in 2:length(W)){
			toreturn = sumlog(toreturn,  ABF1(Z1, Z2, V1, V2, W[i], C)) 
		}
	}
	toreturn = toreturn - log(length(W))
	return(toreturn)
}

avgBF2 = function(Z1, Z2, V1, V2, C = 0, W){
        toreturn = ABF2(Z1, Z2, V1, V2, W[1], C)
        
        if (length(W)>1){
                for (i in 2:length(W)){
                        toreturn = sumlog(toreturn,  ABF2(Z1, Z2, V1, V2, W[i], C)) 
                }
        }
        toreturn = toreturn - log(length(W))
        return(toreturn)
}

avgBF3 = function(Z1, Z2, V1, V2, C = 0, W){
        toreturn = ABF3(Z1, Z2, V1, V2, W[1], C)

        if (length(W)>1){
                for (i in 2:length(W)){
                        toreturn = sumlog(toreturn,  ABF3(Z1, Z2, V1, V2, W[i], C))
                }
        }
        toreturn = toreturn - log(length(W))
        return(toreturn)
}



