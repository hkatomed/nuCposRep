drawResults2.NucII <- function(TYPE, PACKAGE, ymax.HBA = 15, TSS = FALSE){

	LEN <- 1465
	SEQNAME <- "TRP1ARS1"
	par(mfrow = c(2, 1))
	RESULTS <- get(paste("results", PACKAGE, TYPE, sep = "."))
	TITLE <- paste(PACKAGE, ", Seq:", SEQNAME, ", Model:", 
		substr(TYPE, start = 1, stop = 2), sep = "")
	cat(TITLE, "\n")
	
	# NucII
	xlim <- c(800, 1500)
	plot(x = -1000:(LEN+1000-1), y = RESULTS[,3], type = "n", xlim = xlim, 
		ylim = c(-0.35, 1.2), xlab = "Position from to the translation initiation site (bp)", 
		ylab = "Prob./Occup.", xaxt = "n", yaxt = "n")
	title(TITLE)
	polygon(c(-1000, -1000:(LEN+1000-1), (LEN+1000-1)), c(0, RESULTS[,3], 0), col = 8)	# p-dyad
	# points(x = -1000:(LEN+1000-1), y = RESULTS[,4], type = "l", col = 2)
	points(x = -1000:(LEN+1000-1)+73, y = RESULTS[,2], type = "h", col = 4, lwd = 2)
	lines(x = c(0, LEN), y = c(-0.3, -0.3), lwd = 3, col = "gray")
	lines(x = c(115, 789), y = c(-0.2, -0.2), lwd = 3, col = "red")
	lines(x = c(869, 879), y = c(-0.1, -0.1), lwd = 4, col = "blue")
	lines(x = c(847, 859), y = c(-0.1, -0.1), lwd = 4, col = "green4")
	lines(x = c(810, 820), y = c(-0.1, -0.1), lwd = 4, col = "purple")
	lines(x = c(763, 788), y = c(-0.1, -0.1), lwd = 4, col = "black")
	points(x = c(104, 285, 450, 622, 940, 1132, 1308), 
		y = c(1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1), pch = 25, bg = "orange")
	points(x = c(117, 284, 449, 613, 950, 1143, 1329), 
		y = c(1.2, 1.2, 1.2, 1.2, 1.2, 1.2, 1.2), pch = 25, bg = "lightblue")
	if(TSS == TRUE) points(x = c(81), y = c(1.2), pch = 23, bg = "red")
	box(lwd = 2)
	axis(1, lwd = 2)
	axis(2, at = seq(0, 1, 0.5), lwd = 2)
	
	plot(x = -1000:(LEN+1000-1), y = RESULTS[,5], type = "n", xlim = xlim, 
		ylim = c(-23, ymax.HBA), xlab = "Position from to the translation initiation site (bp)", 
		ylab = "HBA", xaxt = "n", yaxt = "n")
	title(TITLE)
	lines(x = -1000:(LEN+1000-1), y = RESULTS[,5], col = 1, lwd = 2)
	lines(x = c(0, LEN), y = c(-22, -22), lwd = 3, col = "gray")
	lines(x = c(115, 789), y = c(-20, -20), lwd = 3, col = "red")
	lines(x = c(869, 879), y = c(-18, -18), lwd = 4, col = "blue")
	lines(x = c(847, 859), y = c(-18, -18), lwd = 4, col = "green4")
	lines(x = c(810, 820), y = c(-18, -18), lwd = 4, col = "purple")
	lines(x = c(763, 788), y = c(-18, -18), lwd = 4, col = "black")
	yvalue <- ymax.HBA-3
	points(x = c(104, 285, 450, 622, 940, 1132, 1308), 
		y = c(yvalue, yvalue, yvalue, yvalue, yvalue, yvalue, yvalue), pch = 25, bg = "orange")
	yvalue <- ymax.HBA-0
	points(x = c(117, 284, 449, 613, 950, 1143, 1329), 
		y = c(yvalue, yvalue, yvalue, yvalue, yvalue, yvalue, yvalue), pch = 25, bg = "lightblue")
	if(TSS == TRUE) points(x = c(81), y = c(yvalue), pch = 23, bg = "red")
	box(lwd = 2)
	axis(1, lwd = 2)
	axis(2, at = seq(-20, ymax.HBA, 10), lwd = 2)
}
