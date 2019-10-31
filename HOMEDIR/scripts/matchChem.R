matchChem <- function(dyad.table, pred.table, chr.num = 16, window = 1){
	require("parallel")
	chr.names <- character(0)
	file.names <- character(0)
	for(i in 1:chr.num){
		chr.name <- paste("chr", i, sep = "")
		chr.names <- c(chr.names, chr.name)
	}

	out <- integer(length = 0)
	for(i in 1:chr.num){
		cat("chr.names[i]: ", chr.names[i], "\n", sep = "")
		dyad.table.chr <- subset(dyad.table, chr == chr.names[i])
		cat("nrow(dyad.table.chr): ", nrow(dyad.table.chr), "\n")
		pred.table.chr <- subset(pred.table, chr == chr.names[i])
		cat("nrow(pred.table.chr): ", nrow(pred.table.chr), "\n")

		out.chr <- integer(length = nrow(pred.table.chr))
		out.chr[1:nrow(pred.table.chr)] <- 0
		for(j in 1:nrow(dyad.table.chr)){
			pos <- dyad.table.chr$pos[j]
			out.chr[(pos-window):(pos+window)] <- 1
		}
		out <- c(out, out.chr)
	}
	return(out)
}
