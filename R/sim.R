library(MASS)

simZ = function(p1, p2, p3, V, lambda, N){
	d = data.frame(matrix(nrow = N, ncol = 5))
	for (i in 1:nrow(d)){
		r = runif(1)
		if (r < p1){
			s = sim1(V)
			d[i,1] = s[1]
			d[i,2] = 1e-5
			d[i,3] = s[2]
			d[i,4] = 1e-5
			d[i,5] = 1
		}
		else if(r < (p1+p2)){
			s = sim2(V)
                        d[i,1] = s[1]
                        d[i,2] = 1e-5
                        d[i,3] = s[2]
                        d[i,4] = 1e-5
			d[i,5] = 2
                }
		else{
			s = sim3(V, lambda)
			d[i,1] = s[1]
                        d[i,2] = 1e-5
                        d[i,3] = s[2]
                        d[i,4] = 1e-5
			d[i,5] = 3
                }

	}
	names(d) = c("Z_1", "V_1", "Z_2", "V_2", "MODEL")
	return(d)

}

sim1 = function(V){
	beta1 = rnorm(1, 0, sd = sqrt(V))
	beta2 = 0
	betahat1 = rnorm(1, beta1, sd = sqrt(1e-5))
	betahat2 = rnorm(1, beta2, sd = sqrt(1e-5))
	Z1 = betahat1/sqrt(1e-5)
	Z2 = betahat2/sqrt(1e-5)
	return(c(Z1, Z2))
}

sim2 = function(V){
	beta2 = rnorm(1, 0, sd = sqrt(V))
	beta1 = 0
	betahat1 = rnorm(1, beta1, sd = sqrt(1e-5))
	betahat2 = rnorm(1, beta2, sd = sqrt(1e-5))
	Z1 = betahat1/sqrt(1e-5)
	Z2 = betahat2/sqrt(1e-5)
	return(c(Z1, Z2))
}
sim3 = function(V, lambda){
	C = matrix(nrow = 2, ncol = 2)
	C[1,1] = V
	C[1,2] = lambda*V
	C[2,1] = lambda*V
	C[2,2] = V
	betas = mvrnorm(mu = c(0, 0), Sigma  = C)
	C[1,1] = 1e-5
	C[1,2] = 0
	C[2,1] = 0
	C[2,2] = 1e-5

	betahats = mvrnorm(mu = betas, Sigma  = C)	
	Z1 = betahats[1]/sqrt(1e-5)
	Z2 = betahats[2]/sqrt(1e-5)
	return(c(Z1, Z2))
}

