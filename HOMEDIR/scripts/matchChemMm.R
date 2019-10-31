matchChemMm <- function(dyad.table, pred.table, mm.chr = "chr19", window = 1){
	require("parallel")
	chr.name <- mm.chr

	out <- integer(length = 0)
	
		cat("chr.name: ", chr.name, "\n", sep = "")
		if(mm.chr == "chr19") chr.len <- 61342430
		dyad.table.chr <- subset(dyad.table, pos > window & pos < (chr.len - window))
		cat("nrow(dyad.table.chr): ", nrow(dyad.table.chr), "\n")
		pred.table.chr <- subset(pred.table, chr == chr.name)
		cat("nrow(pred.table.chr): ", nrow(pred.table.chr), "\n")

		out.chr <- integer(length = nrow(pred.table.chr))
		out.chr[1:nrow(pred.table.chr)] <- 0
		for(j in 1:nrow(dyad.table.chr)){
			pos <- dyad.table.chr$pos[j]
			out.chr[(pos-window):(pos+window)] <- 1
		}
		out <- c(out, out.chr)
	
	return(out)
}
