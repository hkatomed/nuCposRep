wig2ROC <- function(prefix = "NuPoPScSc", nuc.size = 147, 
		save.dir = "../../RData"){

	if(length(grep("ScSc", prefix)) == 1 )	Target <- "Sc"
	if(length(grep("ScSp", prefix)) == 1)	Target <- "Sc"
	if(length(grep("SpSc", prefix)) == 1)	Target <- "Sp"
	if(length(grep("SpSp", prefix)) == 1)	Target <- "Sp"

	wd <- getwd()

	# generate NuPoPScSc.allpred
	getPredTable(prefix = prefix, nuc.size = nuc.size)

	# revise NuPoPScSc.allpred
	MultiMatchChem(prefix = prefix, minwindow = 0, maxwindow = 9)

	# use NuPoPScSc.allpred to generate NuPoPScSc.summary, NuPoPScSc.perf.UniqueW0, NuPoPScSc.perf.RedundantW0
	pred.summary(prefix = prefix)
	setwd(save.dir)
	obj.name <- paste(prefix, ".allpred", sep = "")
	file.name <- paste(prefix, "_allpred.RData", sep = "")
	save(list = obj.name, file = file.name, envir = .GlobalEnv)
	rm(list = obj.name, envir = .GlobalEnv)
	obj.name <- paste(prefix, ".summary", sep = "")
	file.name <- paste(prefix, "_summary.RData", sep = "")
	save(list = obj.name, file = file.name, envir = .GlobalEnv)
	rm(list = obj.name, envir = .GlobalEnv)
	obj.name <- paste(prefix, ".perf.UniqueW0", sep = "")
	file.name <- paste(prefix, "_perf_UniqueW0.RData", sep = "")
	save(list = obj.name, file = file.name, envir = .GlobalEnv)
	rm(list = obj.name, envir = .GlobalEnv)
	obj.name <- paste(prefix, ".perf.RedundantW0", sep = "")
	file.name <- paste(prefix, "_perf_RedundantW0.RData", sep = "")
	save(list = obj.name, file = file.name, envir = .GlobalEnv)
	rm(list = obj.name, envir = .GlobalEnv)
	setwd(wd)

	# generate nature11142_s2.NuPoPScSc.allpred and nature11142_s3.NuPoPScSc.allpred
	getNearestDyad(prefix = prefix, nuc.size = nuc.size)
	setwd(save.dir)
	if(Target == "Sc")	obj.name <- paste("nature11142_s2.", prefix, ".allpred", sep = "")
	if(Target == "Sc")	file.name <- paste("nature11142_s2_", prefix, "_allpred.RData", sep = "")
	if(Target == "Sp")	obj.name <- paste("sd01.", prefix, ".allpred", sep = "")
	if(Target == "Sp")	file.name <- paste("sd01_", prefix, "_allpred.RData", sep = "")
	save(list = obj.name, file = file.name, envir = .GlobalEnv)
	rm(list = obj.name, envir = .GlobalEnv)
	if(Target == "Sc")	obj.name <- paste("nature11142_s3.", prefix, ".allpred", sep = "")
	if(Target == "Sc")	file.name <- paste("nature11142_s3_", prefix, "_allpred.RData", sep = "")
	if(Target == "Sp")	obj.name <- paste("sd02.", prefix, ".allpred", sep = "")
	if(Target == "Sp")	file.name <- paste("sd02_", prefix, "_allpred.RData", sep = "")
	save(list = obj.name, file = file.name, envir = .GlobalEnv)
	rm(list = obj.name, envir = .GlobalEnv)
	setwd(wd)

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

