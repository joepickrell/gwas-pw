source("/Users/jkpickrell/Documents/workspace/gwas-pw/R/Deming_noint.R")

ploteffect = function(d, cols = F){
	d$B1 = d$Z_1 *sqrt(d$V_1)
	d$B2 = d$Z_2 *sqrt(d$V_2)
	rv = mean(d$V_2/ d$V_1)
	#rv = 1
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

	col = "black"
	tmp = d[d$MODEL == 1 | d$MODEL ==4,]
	if (cols){
		col = "green"
	}
	points(tmp$B1, tmp$B2, pch = 20, col = col)
	dm = Deming_noint(tmp$B1, tmp$B2, rv)
	print(dm)
	if (!is.na(dm[1])){ abline(0, dm[1], col = "green")}

	tmp = d[d$MODEL == 2 | d$MODEL =="4_2",]
	
	if (cols){
		col = "blue"
	}
	points(tmp$B1, tmp$B2, pch = 20, col = col)
	dm = Deming_noint(tmp$B1, tmp$B2, rv)
	if (!is.na(dm[1])){ abline(0, dm[1], col = "blue")}
	
	tmp = d[d$MODEL == 3,]
	if (cols){
		col = "red"
	}
	points(tmp$B1, tmp$B2, pch = 20, col = col)
	dm = Deming_noint(tmp$B1, tmp$B2, rv)
        abline(0, dm[1], col = "red")
	legend("topleft", c("M1", "M2", "M3"), pch = 20, col = c("green", "blue", "red"))	
}

ploteffect_lm = function(d){
        d$B1 = d$Z_1 *sqrt(d$V_1)
        d$B2 = d$Z_2 *sqrt(d$V_2)
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
        l1 = summary(lm(tmp$B2 ~ tmp$B1+0))
        abline(0, l1$coef[1], col= "green")

        tmp = d[d$MODEL == 2 | d$MODEL =="4_2",]
        points(tmp$B1, tmp$B2, pch = 20, col = "blue")
        l1 = summary(lm(tmp$B2 ~ tmp$B1+0))
        abline(0, l1$coef[1], col= "blue")

        tmp = d[d$MODEL == 3,]
        points(tmp$B1, tmp$B2, pch = 20, col = "red")
        l1 = summary(lm(tmp$B2 ~ tmp$B1+0))
        abline(0, l1$coef[1], col= "red")

}





