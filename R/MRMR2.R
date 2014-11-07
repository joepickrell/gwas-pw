source("/Users/jkpickrell/Documents/workspace/gwas-pw/R/Deming_noint.R")

rev = function(d){
	tmpz = d$Z_1
	tmpv = d$V_1
	d$Z_1 = d$Z_2
	d$V_1 = d$V_2
	d$Z_2 = tmpz
	d$V_2 = tmpv
	
	if(sum(d$MODEL == 4) >0){
		d[d$MODEL == 4,]$MODEL = 5
	}
	if (sum(d$MODEL == "4_2") > 0){
		d[d$MODEL == "4_2",]$MODEL = 6
	}
	d[d$MODEL == 1,]$MODEL = 7
	d[d$MODEL == 2,]$MODEL = 1
	d[d$MODEL == 7,]$MODEL = 2
	if (sum(d$MODEL == 6)> 0){
		d[d$MODEL == 6,]$MODEL = 4
	}
	if (sum(d$MODEL == 5)>0){
		d[d$MODEL == 5,]$MODEL = "4_2"
	}
	
	return(d)
}
addbetas = function(d){
	d$B1 = d$Z_1*sqrt(d$V_1)
	d$B2 = d$Z_2*sqrt(d$V_2)
	return(d)
}

setexpbetas_1_causes_2 = function(theta, d, lambda = 1){
	beta = theta[1]
	print(paste(beta, lambda))
	d = addbetas(d)
	d$EXPB1 = d$B1
	d$EXPB2 = d$EXPB1*beta
	
	d[d$MODEL == 2 | d$MODEL == "4_2",]$EXPB2 = 0	
	return(list(d = d, theta = theta))
}

setexpbetas_2_causes_1 = function(theta, d, lambda = 1){
        beta = theta[1]
        d = addbetas(d)
        d$EXPB1 = d$B1
        d$EXPB2 = d$EXPB1*beta

        d[d$MODEL == 1 | d$MODEL == "4",]$EXPB2 = 0
        return(list(d = d, theta = theta))
}



setexpbetas_samepheno = function(theta, d, lambda = 1){
        beta = theta[1]
        d = addbetas(d)
        d$EXPB1 = d$B1
        d$EXPB2 = d$EXPB1*beta
        return(list(d = d, theta = theta))
}



setexpbetas_noeffect = function(theta, d, lambda = 1){
        d = addbetas(d)
        d$EXPB1 = 0
        d$EXPB2 = 0
        return(list(d = d, theta = theta))
}



setexpbetas_indep_effects = function(theta, d, lambda = 1){
        beta = theta[1]
        d = addbetas(d)
        d$EXPB1 = d$B1
        d$EXPB2 = d$EXPB1*beta

       	d[d$MODEL == 1 | d$MODEL == "4",]$EXPB2 = 0
        d[d$MODEL == 2 | d$MODEL == "4_2",]$EXPB2 = 0
        return(list(d = d, theta = theta))
}



rankmodels = function(d, lambda = 1){
	o12 = optim_1_causes_2(d, lambda = lambda)
	o21 = optim_2_causes_1(d, lambda = lambda)
	oindep = optim_indep_effects(d, lambda = lambda)
	osame = optim_samepheno(d, lambda = lambda)
	toreturn = data.frame(matrix(nrow = 4, ncol= 7))
	names(toreturn) = c("MODEL", "LLK", "NPARAM", "BETA", "SIGMA2_2", "SIGMA2_E", "AIC")
	toreturn[1,1] = "1->2"
	toreturn[2,1] = "2->1"
	toreturn[3,1] = "indep"
	toreturn[4,1] = "same"
	
	toreturn[1,2] = -o12$o$value
	toreturn[2,2] = -o21$o$value
	toreturn[3,2] = -oindep$o$value
	toreturn[4,2] = -osame$o$value
	
	toreturn[1,3] = 3
	toreturn[2,3] = 2
	toreturn[3,3] = 3
	toreturn[4,3] = 2
	
	toreturn[1,4] = o12$o$par[1]
	toreturn[2,4] = o21$o$par[1]
	toreturn[3,4] = oindep$o$par[1]
	toreturn[4,4] = osame$o$par[1]
        
	toreturn[1,5] = exp(o12$o$par[2])
        toreturn[2,5] = NA
        toreturn[3,5] = exp(oindep$o$par[2])
        toreturn[4,5] = NA


        toreturn[1,6] = exp(o12$o$par[3])
        toreturn[2,6] = exp(o21$o$par[2])
        toreturn[3,6] = exp(oindep$o$par[3])
        toreturn[4,6] = exp(osame$o$par[2])
	toreturn$AIC = 2*toreturn$NPARAM- 2*toreturn$LLK

	toreturn = toreturn[order(toreturn$AIC),]
	return(toreturn)

}

llk_1_causes_2 = function(theta, d, lambda = 1){
       	#
        # class 1 loci:
        #       beta2 ~ N(\beta \mu, \sigma^2_e)
        #
        # class 2 loci:
        #       beta2 ~ N(0, \sigma^2_2 + \sigma^2_e)
        #
        # class 3 loci:
        #       beta2 ~ N(\beta\mu, \sigma^2_e)
	print(theta)
	d = setexpbetas_1_causes_2(theta, d, lambda  = lambda)$d
        sigma22 = exp(theta[2])
        sigma2e = exp(theta[3])
        c1 =d[d$MODEL == 1 | d$MODEL == 4,]
        c2 =d[d$MODEL == 2 | d$MODEL == "4_2",]
        c3 =d[d$MODEL == 3,]

        c1$LLK2 = apply(c1[c("B2", "EXPB2")], 1, FUN = function(x){ return(dnorm(x[1], x[2], sd = sqrt(sigma2e), log = T))})

        c2$LLK2 = apply(c2[c("B2", "EXPB2")], 1, FUN = function(x){ return(dnorm(x[1], x[2], sd = sqrt(sigma22+sigma2e), log = T))})

        c3$LLK2 = apply(c3[c("B2", "EXPB2")], 1, FUN = function(x){ return(dnorm(x[1], x[2], sd = sqrt(sigma2e), log = T))})
        s2 = sum(c1$LLK2) + sum(c2$LLK2)+ sum(c3$LLK2)
	
	llk = s2
	return(-llk)
}



llk_2_causes_1 = function(theta, d, lambda = 1){
        #
        # class 1 loci:
        #       beta2 ~ N(0, \sigma^2_e)
        #
        # class 2 loci:
        #       beta2 ~ N(\beta \mu,  \sigma^2_e)
        #
        # class 3 loci:
        #       beta2 ~ N(\beta\mu, \sigma^2_e)
        print(theta)
        d = setexpbetas_2_causes_1(theta, d, lambda  = lambda)$d
        sigma2e = exp(theta[2])
        c1 =d[d$MODEL == 1 | d$MODEL == 4,]
        c2 =d[d$MODEL == 2 | d$MODEL == "4_2",]
        c3 =d[d$MODEL == 3,]

        c1$LLK2 = apply(c1[c("B2", "EXPB2")], 1, FUN = function(x){ return(dnorm(x[1], x[2], sd = sqrt(sigma2e), log = T))})

        c2$LLK2 = apply(c2[c("B2", "EXPB2")], 1, FUN = function(x){ return(dnorm(x[1], x[2], sd = sqrt(sigma2e), log = T))})

        c3$LLK2 = apply(c3[c("B2", "EXPB2")], 1, FUN = function(x){ return(dnorm(x[1], x[2], sd = sqrt(sigma2e), log = T))})
        s2 = sum(c1$LLK2) + sum(c2$LLK2)+ sum(c3$LLK2)

        llk = s2
        return(-llk)
}





llk_indep_effects = function(theta, d, lambda = 1){
        #
        # class 1 loci:
        #       beta2 ~ N(0, \sigma^2_e)
        #
        # class 2 loci:
        #       beta2 ~ N(0, \sigma^2_2+ \sigma^2_e)
        #
        # class 3 loci:
        #       beta2 ~ N(\mu, \sigma^2_e)
        print(theta)
        d = setexpbetas_indep_effects(theta, d, lambda  = lambda)$d
        sigma22 = exp(theta[2])
        sigma2e = exp(theta[3])
        c1 =d[d$MODEL == 1 | d$MODEL == 4,]
        c2 =d[d$MODEL == 2 | d$MODEL == "4_2",]
        c3 =d[d$MODEL == 3,]

        c1$LLK2 = apply(c1[c("B2", "EXPB2")], 1, FUN = function(x){ return(dnorm(x[1], x[2], sd = sqrt(sigma2e), log = T))})

        c2$LLK2 = apply(c2[c("B2", "EXPB2")], 1, FUN = function(x){ return(dnorm(x[1], x[2], sd = sqrt(sigma22+ sigma2e), log = T))})

        c3$LLK2 = apply(c3[c("B2", "EXPB2")], 1, FUN = function(x){ return(dnorm(x[1], x[2], sd = sqrt(sigma2e), log = T))})
        s2 = sum(c1$LLK2) + sum(c2$LLK2)+ sum(c3$LLK2)

        llk = s2
        return(-llk)
}



llk_samepheno = function(theta, d, lambda = 1){
        #
        # class 1 loci:
        #       beta2 ~ N(\mu, \sigma^2_e)
        #
        # class 2 loci:
        #       beta2 ~ N(\mu, \sigma^2_e)
        #
        # class 3 loci:
        #       beta2 ~ N(\mu, \sigma^2_e)
        print(theta)
        d = setexpbetas_samepheno(theta, d, lambda  = lambda)$d
        sigma2e = exp(theta[2])
        c1 =d[d$MODEL == 1 | d$MODEL == 4,]
        c2 =d[d$MODEL == 2 | d$MODEL == "4_2",]
        c3 =d[d$MODEL == 3,]

        c1$LLK2 = apply(c1[c("B2", "EXPB2")], 1, FUN = function(x){ return(dnorm(x[1], x[2], sd = sqrt(sigma2e), log = T))})

        c2$LLK2 = apply(c2[c("B2", "EXPB2")], 1, FUN = function(x){ return(dnorm(x[1], x[2], sd = sqrt(sigma2e), log = T))})

        c3$LLK2 = apply(c3[c("B2", "EXPB2")], 1, FUN = function(x){ return(dnorm(x[1], x[2], sd = sqrt(sigma2e), log = T))})
        s2 = sum(c1$LLK2) + sum(c2$LLK2)+ sum(c3$LLK2)

        llk = s2
        return(-llk)
}




llk_noeffect = function(theta, d, lambda = 1){
	# 
	# class 1 loci:
	# 	beta2 ~ N(0, \sigma^2_e)
	#
	# class 2 loci:
	#	beta2 ~ N(0, \sigma^2_2 + \sigma^2_e)
	#
	# class 3 loci:
	#	beta2 ~ N(0, \sigma^2_2 + \sigma^2_e)

	
        d = setexpbetas_noeffect(theta, d, lambda = lambda)$d
        sigma22 = exp(theta[1])
        sigma2e = exp(theta[2])
	c1 =d[d$MODEL == 1 | d$MODEL == 4,]
	c2 =d[d$MODEL == 2 | d$MODEL == "4_2",]
	c3 =d[d$MODEL == 3,]
	
        c1$LLK2 = apply(c1[c("B2", "EXPB2")], 1, FUN = function(x){ return(dnorm(x[1], x[2], sd = sqrt(sigma2e), log = T))})
        
        c2$LLK2 = apply(c2[c("B2", "EXPB2")], 1, FUN = function(x){ return(dnorm(x[1], x[2], sd = sqrt(sigma22+sigma2e), log = T))})
        
        c3$LLK2 = apply(c3[c("B2", "EXPB2")], 1, FUN = function(x){ return(dnorm(x[1], x[2], sd = sqrt(sigma22+sigma2e), log = T))})
        s2 = sum(c1$LLK2) + sum(c2$LLK2)+ sum(c3$LLK2)
        llk = s2
        return(-llk)
}


optim_noeffect = function(d, lambda = 1){
	d = addbetas(d)
	init = c(0, 0)
        o = optim(init, llk_noeffect, d = d, method = "Nelder")
        return(list(o = o, d = d))
}


optim_samepheno = function(d, lambda = 1){
        d = addbetas(d)
        init = c(0, 0)
        o = optim(init, llk_samepheno, d = d, method = "Nelder")
        return(list(o = o, d = d))
}



optim_1_causes_2 = function(d, lambda = 1){
        d = addbetas(d)
        init = c(0, 0, 0)
        o = optim(init, llk_1_causes_2, d = d, method = "Nelder")
        return(list(o = o, d = d))
}

optim_2_causes_1 = function(d, lambda = 1){
        d = addbetas(d)
        init = c(0, 0)
        o = optim(init, llk_2_causes_1, d = d, method = "Nelder")
        return(list(o = o, d = d))
}

optim_indep_effects = function(d, lambda = 1){
        d = addbetas(d)
        init = c(0, 0, 0)
        o = optim(init, llk_indep_effects, d = d, method = "Nelder")
        return(list(o = o, d = d))
}






