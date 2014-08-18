
pdf("allc0.pdf", height = 8,width = 8)
par(mfrow = c(3,3))
toplot = read.table("m1_c0.toplot", as.is= T)
toplot =toplot[is.finite(toplot[,8]),]
plot(toplot[,4], toplot[,7], xlab = "log10(BF) [model 1]", ylab = "log10(Stephens BF) [model 1]")
abline(0, 1, col = "grey")
mtext("Simulations under model 1 (C = 0)", adj = 0, line = 0.5)
plot(toplot[,5], toplot[,8], xlab = "log10(BF) [model 2]", ylab = "log10(Stephens BF) [model 2]")
abline(0, 1, col = "grey")
plot(toplot[,6], toplot[,9], xlab = "log10(BF) [model 3]", ylab = "log10(Stephens BF) [model 3]")
abline(0, 1, col = "grey")

toplot = read.table("m2_c0.toplot", as.is= T)
toplot =toplot[is.finite(toplot[,8]),]
plot(toplot[,4], toplot[,7], xlab = "log10(BF) [model 1]", ylab = "log10(Stephens BF) [model 1]")
abline(0, 1, col = "grey")
mtext("Simulations under model 2 (C = 0)", adj = 0, line = 0.5)
plot(toplot[,5], toplot[,8], xlab = "log10(BF) [model 2]", ylab = "log10(Stephens BF) [model 2]")
abline(0, 1, col = "grey")
plot(toplot[,6], toplot[,9], xlab = "log10(BF) [model 3]", ylab = "log10(Stephens BF) [model 3]")
abline(0, 1, col = "grey")


toplot = read.table("m3_c0.toplot", as.is= T)
toplot =toplot[is.finite(toplot[,8]),]
toplot =toplot[is.finite(toplot[,7]),]
toplot =toplot[is.finite(toplot[,9]),]
plot(toplot[,4], toplot[,7], xlab = "log10(BF) [model 1]", ylab = "log10(Stephens BF) [model 1]")
abline(0, 1, col = "grey")
mtext("Simulations under model 3 (C = 0)", adj = 0, line = 0.5)
plot(toplot[,5], toplot[,8], xlab = "log10(BF) [model 2]", ylab = "log10(Stephens BF) [model 2]")
abline(0, 1, col = "grey")
plot(toplot[,6], toplot[,9], xlab = "log10(BF) [model 3]", ylab = "log10(Stephens BF) [model 3]")
abline(0, 1, col = "grey")

dev.off()

