lHBA.rotation <- function(testSeq, seqName, species = "sc", ylim = c(-10, 5)){
	Data <- matrix(nrow = 9, ncol = 13)
	colnames(Data) <- LETTERS[1:13]
	rownames(Data) <- seq(-4,4,1)
	Data <- data.frame(Data, stringsAsFactors = FALSE)
	Data$seq <- ""
	seq1 <- substr(testSeq, start = 144, stop = 147)
	seq2 <- substr(testSeq, start = 1, stop = 143)
	Data$seq[1] <- paste(seq1, seq2, collapse = "", sep = "")
	seq1 <- substr(testSeq, start = 145, stop = 147)
	seq2 <- substr(testSeq, start = 1, stop = 144)
	Data$seq[2] <- paste(seq1, seq2, collapse = "", sep = "")
	seq1 <- substr(testSeq, start = 146, stop = 147)
	seq2 <- substr(testSeq, start = 1, stop = 145)
	Data$seq[3] <- paste(seq1, seq2, collapse = "", sep = "")
	seq1 <- substr(testSeq, start = 147, stop = 147)
	seq2 <- substr(testSeq, start = 1, stop = 146)
	Data$seq[4] <- paste(seq1, seq2, collapse = "", sep = "")
	Data$seq[5] <- testSeq
	seq1 <- substr(testSeq, start = 2, stop = 147)
	seq2 <- substr(testSeq, start = 1, stop = 1)
	Data$seq[6] <- paste(seq1, seq2, collapse = "", sep = "")
	seq1 <- substr(testSeq, start = 3, stop = 147)
	seq2 <- substr(testSeq, start = 1, stop = 2)
	Data$seq[7] <- paste(seq1, seq2, collapse = "", sep = "")
	seq1 <- substr(testSeq, start = 4, stop = 147)
	seq2 <- substr(testSeq, start = 1, stop = 3)
	Data$seq[8] <- paste(seq1, seq2, collapse = "", sep = "")
	seq1 <- substr(testSeq, start = 5, stop = 147)
	seq2 <- substr(testSeq, start = 1, stop = 4)
	Data$seq[9] <- paste(seq1, seq2, collapse = "", sep = "")
	
	require(nuCpos)
	for(i in 1:9){
		Data[i, 1:13] <- localHBA(Data$seq[i], species = species, silent = TRUE)
	}

	# omit terminal values
	Data$A[1:4] <- NA
	Data$M[6:9] <- NA

	x <- 1:13
	main <- paste(seqName, " -4:4, (species=\"", species, "\")", sep = "")
	plot(x = x, y = numeric(length = 13), type = "n", xaxt = "n", yaxt = "n", 
		ylab = "local HBA", ylim = ylim, xlim = c(0.5, 13.5), 
		main = main, bty = "n")

	for(i in 1:13){
		lines(x = seq(i-0.4, i+0.4, 0.1), y = as.numeric(Data[, i]))
	}

	for(i in 1:4){
		points(x = x-0.5+0.1*i, y = as.numeric(Data[i, 1:13]), 
		pch = 20, col = "black")
	}

	for(i in 6:9){
		points(x = x-0.5+0.1*i, y = as.numeric(Data[i, 1:13]), 
		pch = 20, col = "black")
	}
	points(x = x, y = as.numeric(Data[5, 1:13]), 
	pch = 20, col = "red")
	axis(1, lwd = 2, at = 1:13, labels = LETTERS[1:13], tick = FALSE)
	axis(2, lwd = 2)

	for(i in seq(0.5, 13.5, 1)){
		lines(x = c(i, i), y = c(-10, 5))
	}
}

