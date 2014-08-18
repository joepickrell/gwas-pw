source("sim.R")
source("stephens.R")
source("approxBF.R")

cc = 0

p1 = 1 # 1 or 0 for whether alt SNPs influence p1
p2 = 0 # 1 or 0 for whether alt SNPs influence p2

S = matrix(nrow = 2, ncol = 2)
S[1,1] = 1
S[1,2] = cc
S[2,1] = cc
S[2,2] = 1


fnull = 0
betav = 0.3

NSNP = 1000

toplot = data.frame(matrix(nrow = NSNP, ncol = 13))

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
	toplot[i,4:6] = abf[1:3]/log(10)
	toplot[i,10:13] = abf[4:7]
	toplot[i,7:9] = sBF$lbf[c(2,4,5),1]
}

pdf("m1_c0.pdf", height = 3.5)
par(mfrow = c(1, 3))

plot(toplot[,4], toplot[,7], xlab = "log10(BF) [model 1]", ylab = "log10(Stephens BF) [model 1]")
abline(0, 1, col = "grey")
mtext("Simulations under model 1 (C = 0)", adj = 0, line = 0.5)
plot(toplot[,5], toplot[,8], xlab = "log10(BF) [model 2]", ylab = "log10(Stephens BF) [model 2]")
abline(0, 1, col = "grey")
plot(toplot[,6], toplot[,9], xlab = "log10(BF) [model 3]", ylab = "log10(Stephens BF) [model 3]")
abline(0, 1, col = "grey")

dev.off()

write.table(toplot, file = "m1_c0.toplot", quote = F, row.names = F, col.names = F)
