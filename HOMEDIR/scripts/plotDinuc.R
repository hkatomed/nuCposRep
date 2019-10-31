plotDinuc <- function(dinuc = c("CC", "CG", "GC", "GG"), 
		COL = c("blue", "purple", "magenta", "black"), species = "sc"){
	if(species == "sc")	two.df <- Sc.two.df
	if(species == "sp")	two.df <- Sp.two.df
	if(species == "mm")	two.df <- Mm.two.df

	x <- 1:146
	y <- matrix(nrow = length(dinuc), ncol = 146)
	for(i in 1:length(dinuc)) y[i,] <- two.df[[dinuc[i]]]
	ymax <- max(y) +0.01
	ymin <- min(y) -0.01

	par(mfrow = c(4, 1))
	for(i in 1:length(dinuc)){
		main <- paste("Dinucleotide frequency: ", dinuc[i], 
			" (mean: ", round(mean(y[i,]), digits = 4), ", ", 
			"species = ", species,")", sep = "") 

		plot(x = x, y = y[1,], type = "n", ylim = c(ymin, ymax), 
			main = main, xaxt = "n", yaxt = "n", 
			ylab = "Frequency", xlab = "Nucleosomal position (bp)")
		abline(v = 73.5, lwd = 2, col = "gray")
		SHLs <- seq(5.25, 80, 10.5)
		abline(v = c(73.5-SHLs, 73.5+SHLs), lty = 2, lwd = 2, col = "gray")
		abline(h = mean(y[i,]), lwd = 2, col = "orange")
		lines(x = x, y = y[i,], lty = 1, lwd = 2, col = COL[i])
		box(lwd = 2)
		axis(side = 1, lwd.ticks = 2, at = c(1, seq(20, 140, 20)))
		axis(side = 2, lwd.ticks = 2, at = c(round(max(y), digits = 2), round(min(y), digits = 2)))
	}	
}
