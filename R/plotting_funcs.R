library(MethComp)

ploteffect = function(d){
	d$B1 = d$Z_1 *sqrt(d$V_1)
	d$B2 = d$Z_2 *sqrt(d$V_2)
	rv = mean(d$V_2/ d$V_1)
	plot(d$B1, d$B2)
	lines( c(0, 0), c(-1000, 10000), col = "grey")
	lines( c(-10000, 10000), c(0, 0), col = "grey")


        for (i in 1:nrow(d)){
                B1 = d[i,]$B1
                sd1 = sqrt(d[i,]$V_1)
                B2 = d[i,]$B2
                sd2 = sqrt(d[i,]$V_2)
                lines( c( B1 - sd1, B1 + sd1), c(B2, B2))
                lines( c( B1, B1), c(B2-sd2, B2+sd2))

        }

	tmp = d[d$MODEL == 1 | d$MODEL ==4,]
	points(tmp$B1, tmp$B2, pch = 20, col = "green")
	#l1 = summary(lm(tmp$B2 ~ tmp$B1+0))
	#abline(0, l1$coef[1], col= "green")
	dm = Deming(tmp$B1, tmp$B2, rv)
	abline(dm[1], dm[2], col = "green")

	tmp = d[d$MODEL == 2 | d$MODEL =="4_2",]
	points(tmp$B1, tmp$B2, pch = 20, col = "blue")
	dm = Deming(tmp$B1, tmp$B2, rv)
        abline(dm[1], dm[2], col = "blue")
	#l1 = summary(lm(tmp$B2 ~ tmp$B1+0))
	#abline(0, l1$coef[1], col= "blue")
	
	tmp = d[d$MODEL == 3,]
	points(tmp$B1, tmp$B2, pch = 20, col = "red")
	dm = Deming(tmp$B1, tmp$B2, rv)
        abline(dm[1], dm[2], col = "red")
	#l1 = summary(lm(tmp$B2 ~ tmp$B1+0))
	#abline(0, l1$coef[1], col= "red")
	
}



