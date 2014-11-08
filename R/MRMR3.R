library(mvtnorm)



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



addlk = function(d, C = 0, W = 0.1){
	wz1 = which(names(d) == "Z_1")
        wz2 = which(names(d) == "Z_2")
        wv1 = which(names(d) == "V_1")
        wv2 = which(names(d) == "V_2")
        d$LK1 = apply(d[,c(wz1, wz2, wv1, wv2)], 1, FUN = function(x){ return( lk_M1(x[1], x[2], x[3], x[4], W))})
        d$LK2 = apply(d[,c(wz1, wz2, wv1, wv2)], 1, FUN = function(x){ return( lk_M2(x[1], x[2], x[3], x[4], W))})
        d$LK3 = apply(d[,c(wz1, wz2, wv1, wv2)], 1, FUN = function(x){ return( lk_M3(x[1], x[2], x[3], x[4], W, C))})
        return(d)



}

sumlog = function(logx, logy){
       if (logx > logy) {
                return(logx + log(1 + exp(logy-logx)))
        }
        else {
                return (logy + log(1 + exp(logx-logy)))
        }
}




rankmodels = function(d){
        o12 = optim_1_causes_2(d)
        o21 = optim_2_causes_1(d)
        oindep = optim_indep_effects(d)
        osame = optim_samepheno(d)
	luncor = llk_indep_effects_uncor(d)
        toreturn = data.frame(matrix(nrow = 5, ncol= 6))
	nobs = nrow(d)	
        names(toreturn) = c("MODEL", "LLK", "NPARAM", "COR", "AIC", "AICc")
        toreturn[1,1] = "1->2"
        toreturn[2,1] = "2->1"
        toreturn[3,1] = "indep"
        toreturn[5,1] = "indep(uncor)"
        toreturn[4,1] = "same"

        toreturn[1,2] = -o12$o$objective
        toreturn[2,2] = -o21$o$objective
        toreturn[3,2] = -oindep$o$objective
        toreturn[4,2] = -osame$o$objective
	toreturn[5,2] = -luncor

        toreturn[1,3] = 7
        toreturn[2,3] = 7
        toreturn[3,3] = 7
        toreturn[4,3] = 7
	toreturn[5,3] = 6
	
        toreturn[1,4] = o12$o$minimum
        toreturn[2,4] = o21$o$minimum
        toreturn[3,4] = oindep$o$minimum
        toreturn[4,4] = osame$o$minimum
	toreturn[5,4] = NA
	
	toreturn$AIC = 2*toreturn$NPARAM - 2*toreturn$LLK
	toreturn$AICc = toreturn$AIC + (2*toreturn$NPARAM *(toreturn$NPARAM+1)) / (nobs- toreturn$NPARAM-1)
	toreturn = toreturn[order(toreturn$AIC),]
        return(toreturn)

}



optim_2_causes_1 = function(d){
        o = optimize(llk_2_causes_1, d = d, interval = c(-0.999, 0.999))
        #o = optim(init, llk_2_causes_1, d = d, method = "L-BFGS-B", lower = c(-0.999, -Inf, -Inf, -Inf), upper = c(0.999, Inf, Inf, Inf))
        return(list(o = o, d = d))
}



optim_1_causes_2 = function(d){
 	o = optimize(llk_1_causes_2, d = d, interval = c(-0.999, 0.999))
        #o = optim(init, llk_1_causes_2, d = d, method = "L-BFGS-B", lower = c(-0.999, -Inf, -Inf, -Inf), upper = c(0.999, Inf, Inf, Inf))
        return(list(o = o, d = d))
}


optim_indep_effects = function(d){
        #init = c(0, 0, 0, 0)
	o = optimize(llk_indep_effects, d = d, interval = c(-0.999, 0.999))
        #o = optim(init, llk_indep_effects, d = d, method = "L-BFGS-B", lower = c(-0.999, -Inf, -Inf, -Inf), upper = c(0.999, Inf, Inf, Inf))
        return(list(o = o, d = d))
}




optim_samepheno = function(d){
        #init = c(0, 0, 0)
	o = optimize(llk_samepheno, d = d, interval = c(-0.999, 0.999))
        #o = optim(init, llk_samepheno, d = d, method = "L-BFGS-B", lower = c(-0.999, -Inf, -Inf), upper = c(0.999, Inf, Inf))
        return(list(o = o, d = d))
}




llk_2_causes_1 = function(theta, d){
      	#
        # class 1 loci:
        #       beta ~ MVN(c(0, 0), V1 = \sigma^2_1 +\sigma^2_e, V2 = \sigma^2_e, C = 0)
        #
        # class 2 loci:
        #       beta ~ MVN( c(0, 0), V1 = \sigma^2_1 +\sigma^2_e, V2 = \sigma^2_2+\sigma^2_e, C = C*sqrt(V1*V2))
        #
        # class 3 loci:
        #       beta ~ MVN( c(0, 0), V1 = \sigma^2_1 +\sigma^2_e, V2 = \sigma^2_2+\sigma^2_e, C = c*sqrt(V1*V2))
        d = addbetas(d)
	beta = theta[1]
	S = matrix(nrow = 2, ncol = 2)
        c1 =d[d$MODEL == 1 | d$MODEL == 4,]
        c2 =d[d$MODEL == 2 | d$MODEL == "4_2",]
	c3 = d[d$MODEL ==3,]

        V1 = var(c1$B1)
        V2 = var(c1$B2)
	COV = 0
        S[1,1] = V1
        S[1,2] = COV
        S[2,1] = COV
        S[2,2] = V2
        c1$LLK = apply(c1[c("B1", "B2")], 1, FUN = function(x){ return(dmvnorm( c(x[1], x[2]), mean = c(0,0), sigma = S, log = T))})

        V1 = var(c2$B1)
        V2 = var(c2$B2)
	COV = beta*sqrt(V1*V2)
        S[1,1] = V1
        S[1,2] = COV
        S[2,1] = COV
        S[2,2] = V2
        c2$LLK = apply(c2[c("B1", "B2")], 1, FUN = function(x){ return(dmvnorm( c(x[1], x[2]), mean = c(0,0), sigma = S, log = T))})
        
        V1 = var(c3$B1)
        V2 = var(c3$B2)
        COV = beta*sqrt(V1*V2)
        S[1,1] = V1
        S[1,2] = COV
        S[2,1] = COV
        S[2,2] = V2
        c3$LLK = apply(c3[c("B1", "B2")], 1, FUN = function(x){ return(dmvnorm( c(x[1], x[2]), mean = c(0,0), sigma = S, log = T))})
        
	s2 = sum(c1$LLK) + sum(c2$LLK) + sum(c3$LLK)

        llk = s2
        return(-llk)
}


llk_1_causes_2 = function(theta, d){
        #
        # class 1 loci:
        #       beta ~ MVN(c(0, 0), V1 = \sigma^2_1 +\sigma^2_e, V2 = \sigma^2_2+ \sigma^2_e, C = c*sqrt(V1*V2))
        #
        # class 2 loci:
        #       beta ~ MVN( c(0, 0), V1 = \sigma^2_e, V2 = \sigma^2_2+\sigma^2_e, C = 0)
        #
        # class 3 loci:
        #       beta ~ MVN( c(0, 0), V1 = \sigma^2_1 + \sigma^2_e, V2 = \sigma^2_2+\sigma^2_e, C = c*sqrt(V1*V2))

        d = addbetas(d)
        beta = theta[1]
        S = matrix(nrow = 2, ncol = 2)
        c1 =d[d$MODEL == 1 | d$MODEL == 4,]
        c2 =d[d$MODEL == 2 | d$MODEL == "4_2",]
        c3 = d[d$MODEL ==3,]

        V1 = var(c1$B1)
        V2 = var(c1$B2)
        COV = beta*sqrt(V1*V2)
        S[1,1] = V1
        S[1,2] = COV
        S[2,1] = COV
        S[2,2] = V2
        c1$LLK = apply(c1[c("B1", "B2")], 1, FUN = function(x){ return(dmvnorm( c(x[1], x[2]), mean = c(0,0), sigma = S, log = T))})

        V1 = var(c2$B1)
        V2 = var(c2$B2)
        COV = 0
        S[1,1] = V1
        S[1,2] = COV
        S[2,1] = COV
        S[2,2] = V2
        c2$LLK = apply(c2[c("B1", "B2")], 1, FUN = function(x){ return(dmvnorm( c(x[1], x[2]), mean = c(0,0), sigma = S, log = T))})

        V1 = var(c3$B1)
        V2 = var(c3$B2)
        COV = beta*sqrt(V1*V2)
        S[1,1] = V1
        S[1,2] = COV
        S[2,1] = COV
        S[2,2] = V2
        c3$LLK = apply(c3[c("B1", "B2")], 1, FUN = function(x){ return(dmvnorm( c(x[1], x[2]), mean = c(0,0), sigma = S, log = T))})

        s2 = sum(c1$LLK) + sum(c2$LLK) + sum(c3$LLK)

        llk = s2
        return(-llk)

}


llk_indep_effects = function(theta, d){
        #
        # class 1 loci:
        #       beta ~ MVN(c(0, 0), V1 = \sigma^2_1 +\sigma^2_e, V2 = \sigma^2_e, C = 0)
        #
        # class 2 loci:
        #       beta ~ MVN( c(0, 0), V1 = \sigma^2_e, V2 = \sigma^2_2+\sigma^2_e, C = 0)
        #
        # class 3 loci:
        #       beta ~ MVN( c(0, 0), V1 = \sigma^2_1 + \sigma^2_e, V2 = \sigma^2_2+\sigma^2_e, C = c*sqrt(V1*V2))


        d = addbetas(d)
        beta = theta[1]
        S = matrix(nrow = 2, ncol = 2)
        c1 =d[d$MODEL == 1 | d$MODEL == 4,]
        c2 =d[d$MODEL == 2 | d$MODEL == "4_2",]
        c3 = d[d$MODEL ==3,]

        V1 = var(c1$B1)
        V2 = var(c1$B2)
        COV = 0
        S[1,1] = V1
        S[1,2] = COV
        S[2,1] = COV
        S[2,2] = V2
        c1$LLK = apply(c1[c("B1", "B2")], 1, FUN = function(x){ return(dmvnorm( c(x[1], x[2]), mean = c(0,0), sigma = S, log = T))})

        V1 = var(c2$B1)
        V2 = var(c2$B2)
        COV = 0
        S[1,1] = V1
        S[1,2] = COV
        S[2,1] = COV
        S[2,2] = V2
        c2$LLK = apply(c2[c("B1", "B2")], 1, FUN = function(x){ return(dmvnorm( c(x[1], x[2]), mean = c(0,0), sigma = S, log = T))})

        V1 = var(c3$B1)
        V2 = var(c3$B2)
        COV = beta*sqrt(V1*V2)
        S[1,1] = V1
        S[1,2] = COV
        S[2,1] = COV
        S[2,2] = V2
        c3$LLK = apply(c3[c("B1", "B2")], 1, FUN = function(x){ return(dmvnorm( c(x[1], x[2]), mean = c(0,0), sigma = S, log = T))})

        s2 = sum(c1$LLK) + sum(c2$LLK) + sum(c3$LLK)

        llk = s2
        return(-llk)


}



llk_indep_effects_uncor = function(d){
        #
        # class 1 loci:
        #       beta ~ MVN(c(0, 0), V1 = \sigma^2_1 +\sigma^2_e, V2 = \sigma^2_e, C = 0)
        #
        # class 2 loci:
        #       beta ~ MVN( c(0, 0), V1 = \sigma^2_e, V2 = \sigma^2_2+\sigma^2_e, C = 0)
        #
        # class 3 loci:
        #       beta ~ MVN( c(0, 0), V1 = \sigma^2_1 + \sigma^2_e, V2 = \sigma^2_2+\sigma^2_e, C = c*sqrt(V1*V2))


        d = addbetas(d)
        S = matrix(nrow = 2, ncol = 2)
        c1 =d[d$MODEL == 1 | d$MODEL == 4,]
        c2 =d[d$MODEL == 2 | d$MODEL == "4_2",]
        c3 = d[d$MODEL ==3,]

        V1 = var(c1$B1)
        V2 = var(c1$B2)
        COV = 0
        S[1,1] = V1
        S[1,2] = COV
        S[2,1] = COV
        S[2,2] = V2
        c1$LLK = apply(c1[c("B1", "B2")], 1, FUN = function(x){ return(dmvnorm( c(x[1], x[2]), mean = c(0,0), sigma = S, log = T))})

        V1 = var(c2$B1)
        V2 = var(c2$B2)
        COV = 0
        S[1,1] = V1
        S[1,2] = COV
        S[2,1] = COV
        S[2,2] = V2
        c2$LLK = apply(c2[c("B1", "B2")], 1, FUN = function(x){ return(dmvnorm( c(x[1], x[2]), mean = c(0,0), sigma = S, log = T))})

        V1 = var(c3$B1)
        V2 = var(c3$B2)
        COV = 0
        S[1,1] = V1
        S[1,2] = COV
        S[2,1] = COV
        S[2,2] = V2
        c3$LLK = apply(c3[c("B1", "B2")], 1, FUN = function(x){ return(dmvnorm( c(x[1], x[2]), mean = c(0,0), sigma = S, log = T))})

        s2 = sum(c1$LLK) + sum(c2$LLK) + sum(c3$LLK)

        llk = s2
        return(-llk)


}



llk_samepheno = function(theta, d){
        #
        # class 1 loci:
        #       beta ~ MVN( c(0, 0), V1 = \sigma^2_1, V2 = \sigma^2_2, C = C*sqrt(V1*V2))
        #
        # class 2 loci:
        #       beta ~ MVN( c(0, 0), V1 = \sigma^2_1, V2 = \sigma^2_2, C = C*sqrt(V1*V2))
        #
        # class 3 loci:
        #       beta ~ MVN( c(0, 0), V1 = \sigma^2_1, V2 = \sigma^2_2, C = c*sqrt(V1*V2))



        d = addbetas(d)
        beta = theta[1]
        S = matrix(nrow = 2, ncol = 2)
        c1 =d[d$MODEL == 1 | d$MODEL == 4,]
        c2 =d[d$MODEL == 2 | d$MODEL == "4_2",]
        c3 = d[d$MODEL ==3,]

        V1 = var(c1$B1)
        V2 = var(c1$B2)
        COV =  beta*sqrt(V1*V2)
        S[1,1] = V1
        S[1,2] = COV
        S[2,1] = COV
        S[2,2] = V2
        c1$LLK = apply(c1[c("B1", "B2")], 1, FUN = function(x){ return(dmvnorm( c(x[1], x[2]), mean = c(0,0), sigma = S, log = T))})

        V1 = var(c2$B1)
        V2 = var(c2$B2)
        COV =  beta*sqrt(V1*V2)
        S[1,1] = V1
        S[1,2] = COV
        S[2,1] = COV
        S[2,2] = V2
        c2$LLK = apply(c2[c("B1", "B2")], 1, FUN = function(x){ return(dmvnorm( c(x[1], x[2]), mean = c(0,0), sigma = S, log = T))})

        V1 = var(c3$B1)
        V2 = var(c3$B2)
        COV = beta*sqrt(V1*V2)
        S[1,1] = V1
        S[1,2] = COV
        S[2,1] = COV
        S[2,2] = V2
        c3$LLK = apply(c3[c("B1", "B2")], 1, FUN = function(x){ return(dmvnorm( c(x[1], x[2]), mean = c(0,0), sigma = S, log = T))})

        s2 = sum(c1$LLK) + sum(c2$LLK) + sum(c3$LLK)

        llk = s2
        return(-llk)


}


