getNearestDyadMm <- function(prefix = "NuPoPScSc", nuc.size = 147, mm.chr = "chr19"){
	
	obj.name <- paste("mm9_", mm.chr, "_Uni", sep = "")
	dyad.table.unique <- get(obj.name)
	dyad.table.unique <- data.frame(pos = dyad.table.unique)
	obj.name <- paste("mm9_", mm.chr, "_Red", sep = "")
	dyad.table.redundant <- get(obj.name)
	dyad.table.redundant <- data.frame(pos = dyad.table.redundant)

	out1 <- nearestDyadPredMm(dyad.table = dyad.table.unique, 
			prefix = prefix, nuc.size = nuc.size, mm.chr = mm.chr)
	obj.name <- paste("mm9_", mm.chr, "_Uni", prefix, ".allpred", sep = "")
	assign(obj.name, out1, envir = .GlobalEnv)

	out2 <- nearestDyadPredMm(dyad.table = dyad.table.redundant, 
			prefix = prefix, nuc.size = nuc.size, mm.chr = mm.chr)
	obj.name <- paste("mm9_", mm.chr, "_Red", prefix, ".allpred", sep = "")
	assign(obj.name, out2, envir = .GlobalEnv)
}
