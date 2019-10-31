## これは wig2ROC の一部を動かすための関数。最終的には削除すべき

wig2ROCVit <- function(prefix = "NuPoPScSc", nuc.size = 147, 
		save.dir = "../../RData"){

	if(length(grep("ScSc", prefix)) == 1 )	Target <- "Sc"
	if(length(grep("ScSp", prefix)) == 1)	Target <- "Sc"
	if(length(grep("SpSc", prefix)) == 1)	Target <- "Sp"
	if(length(grep("SpSp", prefix)) == 1)	Target <- "Sp"

	wd <- getwd()

	# generate nature11142_s2.NuPoPScSc and nature11142_s3.NuPoPScSc
	getNearestDyadVit(prefix = prefix, nuc.size = nuc.size)
	setwd(save.dir)
	if(Target == "Sc")	obj.name <- paste("nature11142_s2.", prefix, sep = "")
	if(Target == "Sc")	file.name <- paste("nature11142_s2_", prefix, ".RData", sep = "")
	if(Target == "Sp")	obj.name <- paste("sd01.", prefix, sep = "")
	if(Target == "Sp")	file.name <- paste("sd01_", prefix, ".RData", sep = "")
	save(list = obj.name, file = file.name, envir = .GlobalEnv)
	rm(list = obj.name, envir = .GlobalEnv)
	if(Target == "Sc")	obj.name <- paste("nature11142_s3.", prefix, sep = "")
	if(Target == "Sc")	file.name <- paste("nature11142_s3_", prefix, ".RData", sep = "")
	if(Target == "Sp")	obj.name <- paste("sd02.", prefix, sep = "")
	if(Target == "Sp")	file.name <- paste("sd02_", prefix, ".RData", sep = "")
	save(list = obj.name, file = file.name, envir = .GlobalEnv)
	rm(list = obj.name, envir = .GlobalEnv)
	setwd(wd)

}

