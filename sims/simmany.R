source("sim.R")
source("stephens.R")
source("approxBF.R")

a =commandArgs(T)
cc = a[1]

p1 = a[2] # 1 or 0 for whether alt SNPs influence p1
p2 = a[3] # 1 or 0 for whether alt SNPs influence p2

S = matrix(nrow = 2, ncol = 2)
S[1,1] = 1
S[1,2] = cc
S[2,1] = cc
S[2,2] = 1


fnull = 0.8
betav = 0.3

NSNP = 1000

toplot = data.frame(matrix(nrow = NSNP, ncol = 9))

for (i in 1:NSNP){
	print(i)
	f = rbeta(1, 2, 2)
	b1 = 0
	b2 = 0 
	r = runif(1)
	if ( r > fnull & p1 == 1 ){
		b1 = rnorm(1, 0, betav)
	}
	if ( r > fnull & p2 == 1){
		b2 = rnorm(1, 0, betav)
	}
	toplot[i,1] = f
	toplot[i,2] = b1
	toplot[i,3] = b2
	s = simgp(5000, c(b1, b2), f, S)
	sBF = logBF.rankone.matrix(s$g, s$p, c(0.5, 0.5))
	abf = approxBF(s$g, s$p)
	toplot[i,4:6] = abf/log(10)
	toplot[i,7:9] = sBF$lbf[c(2,4,5),1]
}


