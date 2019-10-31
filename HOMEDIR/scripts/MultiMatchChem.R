MultiMatchChem <- function(prefix = "NuPoPScSc", minwindow = 0, maxwindow = 0){
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
	pred.table <- get(paste(prefix, ".allpred", sep = ""), envir = .GlobalEnv)

	for(i in minwindow:maxwindow){
		cat("window: ", i, "\n")
		Colname <- paste("UniqueW", i, sep = "")
		pred.table[Colname] <- matchChem(dyad.table = dyad.table.unique, 
				pred.table = pred.table, 
				chr.num = chr.num, window = i)
		Colname <- paste("RedundantW", i, sep = "")
		pred.table[Colname] <- matchChem(dyad.table = dyad.table.redundant, 
				pred.table = pred.table, 
				chr.num = chr.num, window = i)
	}
	obj.name <- paste(prefix, ".allpred", sep = "")
	assign(obj.name, pred.table, envir = .GlobalEnv)
	#return(pred.table)
}
