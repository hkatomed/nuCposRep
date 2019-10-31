MultiMatchChemMm <- function(prefix = "NuPoPMmMm", minwindow = 0, maxwindow = 0, mm.chr = mm.chr){

	obj.name <- paste("mm9_", mm.chr, "_Uni", sep = "")
	dyad.table.unique <- get(obj.name)
	dyad.table.unique <- data.frame(pos = dyad.table.unique)
	# obj.name <- paste("mm9_", mm.chr, "_Red", sep = "")
	obj.name <- paste("mm9_", mm.chr, "_Red30", sep = "")
	dyad.table.redundant <- get(obj.name)
	dyad.table.redundant <- data.frame(pos = dyad.table.redundant)

	pred.table <- get(paste(prefix, ".", mm.chr, ".allpred", sep = ""), envir = .GlobalEnv)

	for(i in minwindow:maxwindow){
		cat("window: ", i, "Unique\n")
		Colname <- paste("UniqueW", i, sep = "")
		pred.table[Colname] <- matchChemMm(dyad.table = dyad.table.unique, 
				pred.table = pred.table, 
				mm.chr = mm.chr, window = i)
		cat("window: ", i, "Redundant\n")
		Colname <- paste("RedundantW", i, sep = "")
		pred.table[Colname] <- matchChemMm(dyad.table = dyad.table.redundant, 
				pred.table = pred.table, 
				mm.chr = mm.chr, window = i)
	}
	obj.name <- paste(prefix, ".", mm.chr, ".allpred", sep = "")
	assign(obj.name, pred.table, envir = .GlobalEnv)
	#return(pred.table)
}
