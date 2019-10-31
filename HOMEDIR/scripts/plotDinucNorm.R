plotDinucNorm <- function(dinuc = "CG"){
	x <- 1:146
	sc.y <- Sc.twoNorm.df[[dinuc]]
	sp.y <- Sp.twoNorm.df[[dinuc]]
	mm.y <- Mm.twoNorm.df[[dinuc]]
	ymax <- max(sc.y, sp.y, mm.y) + 0.2
	ymin <- min(sc.y, sp.y, mm.y) - 0.2
	plot(x = x, y = sc.y, type = "n", ylim = c(ymin, ymax), 
		main = paste("Normalized ", dinuc, " frequency", sep =""), xaxt = "n", 
		ylab = "Frequency", xlab = "Nucleosomal position (bp)")
	abline(v = 73.5, lwd = 2, col = "gray")
	SHLs <- seq(5.25, 80, 10.5)
	abline(v = c(73.5-SHLs, 73.5+SHLs), lty = 2, lwd = 2, col = "gray")
	lines(x = x, y = mm.y, lty = 1, lwd = 2, col = "black")
	lines(x = x, y = sc.y, lty = 1, lwd = 2, col = "green")
	lines(x = x, y = sp.y, lty = 1, lwd = 1, col = "purple")
	points(x = x, y = sp.y, pch = 1, col = "purple", cex = 0.5)
	box(lwd = 2)
	axis(side = 1, lwd.ticks = 2, at = c(1, seq(20, 140, 20)))
	axis(side = 2, lwd.ticks = 2)
}
