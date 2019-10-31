getNearestDyadVit <- function(prefix = "NuPoPScSc", nuc.size = 147){
	if(length(grep("ScSc", prefix)) == 1 )	Target <- "Sc"
	if(length(grep("ScSp", prefix)) == 1)	Target <- "Sc"
	if(length(grep("SpSc", prefix)) == 1)	Target <- "Sp"
	if(length(grep("SpSp", prefix)) == 1)	Target <- "Sp"
	if(Target == "Sc"){
		dyad.table.unique <- nature11142_s2.sacCer3
		dyad.table.redundant <- nature11142_s3.sacCer3
		chr.num = 16
	}
	if(Target == "Sp"){
		dyad.table.unique <- sd01
		dyad.table.redundant <- sd02
		chr.num = 3
	}
	
	out1 <- nearestDyadVit(dyad.table = dyad.table.unique, 
			chr.num = chr.num, prefix = prefix, nuc.size = nuc.size)
	if(Target == "Sc")	obj.name <- paste("nature11142_s2.", prefix, sep = "")
	if(Target == "Sp")	obj.name <- paste("sd01.", prefix, sep = "")
	assign(obj.name, out1, envir = .GlobalEnv)

	out2 <- nearestDyadVit(dyad.table = dyad.table.redundant, 
			chr.num = chr.num, prefix = prefix, nuc.size = nuc.size)
	if(Target == "Sc")	obj.name <- paste("nature11142_s3.", prefix, sep = "")
	if(Target == "Sp")	obj.name <- paste("sd02.", prefix, sep = "")
	assign(obj.name, out2, envir = .GlobalEnv)
}
