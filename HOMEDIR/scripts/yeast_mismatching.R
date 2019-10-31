setwd("summary_RData")

## cross-species models
COL <- adjustcolor(c("red", "orange", "gray"), alpha.f = 1)
HEIGHT <- t(as.matrix(num.Summary[,c(5,7,8)]))[,c(2,5,4,7)]
colnames(HEIGHT) <- c("Msc", "Msp", "Csc", "Csp")
barplot(height = HEIGHT, 
	main = "Cross-species models", 
	ylab = "Number of Viterbi dyads", ylim = c(0, 75000), yaxt = "n", 
	xlab = "Parameters", col = COL, lwd = 2)
box(lwd = 2)
axis(side = 2, lwd.ticks = 2, at = c(0, 25000, 50000, 75000), labels = c("0", "25k", "50k", "75k"))
