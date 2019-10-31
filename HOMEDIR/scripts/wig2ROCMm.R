wig2ROCMm <- function(prefix = "NuPoPMmMm", nuc.size = 147, 
		save.dir = "../../RData", mm.chr = "chr19"){

	wd <- getwd()

	# generate NuPoPMmMm.allpred
	getPredTableMm(prefix = prefix, nuc.size = nuc.size, mm.chr = mm.chr)

	# revise NuPoPMmMm.allpred
	MultiMatchChemMm(prefix = prefix, minwindow = 0, maxwindow = 9, mm.chr = mm.chr)

	# use NuPoPMmMm.allpred to generate NuPoPMmMm.summary, NuPoPMmMm.perf.UniqueW0, NuPoPMmMm.perf.RedundantW0
	pred.summaryMm(prefix = prefix, mm.chr = mm.chr)

	setwd(save.dir)
	obj.name <- paste(prefix, ".", mm.chr, ".allpred", sep = "")
	file.name <- paste(prefix, "_", mm.chr, "_allpred.RData", sep = "")
	save(list = obj.name, file = file.name, envir = .GlobalEnv)
	rm(list = obj.name, envir = .GlobalEnv)
	obj.name <- paste(prefix, ".", mm.chr, ".summary", sep = "")
	file.name <- paste(prefix, "_", mm.chr, "_summary.RData", sep = "")
	save(list = obj.name, file = file.name, envir = .GlobalEnv)
	rm(list = obj.name, envir = .GlobalEnv)
	obj.name <- paste(prefix, ".", mm.chr, ".perf.UniqueW0", sep = "")
	file.name <- paste(prefix, "_", mm.chr, "_perf_UniqueW0.RData", sep = "")
	save(list = obj.name, file = file.name, envir = .GlobalEnv)
	rm(list = obj.name, envir = .GlobalEnv)
	obj.name <- paste(prefix, ".", mm.chr, ".perf.RedundantW0", sep = "")
	file.name <- paste(prefix, "_", mm.chr, "_perf_RedundantW0.RData", sep = "")
	save(list = obj.name, file = file.name, envir = .GlobalEnv)
	rm(list = obj.name, envir = .GlobalEnv)
	setwd(wd)

	## 以下の getNearestDyadMm と getNearestDyadVitMm および保存は現実的でない。
	## もしどうしてもやるなら、Unique nucleosomes のみに絞るべき。（参考：20190304_ROC_1.txt）
	# # generate mm9_chr19_Uni.NuPoPMmMm.allpred and mm9_chr19_Red.NuPoPMmMm.allpred
	# getNearestDyadMm(prefix = prefix, nuc.size = nuc.size, mm.chr = mm.chr)
	# 
	# setwd(save.dir)
	# obj.name <- paste("mm9_", mm.chr, "_Uni", prefix, ".allpred", sep = "")
	# file.name <- paste("mm9_", mm.chr, "_Uni", prefix, "_allpred.RData", sep = "")
	# save(list = obj.name, file = file.name, envir = .GlobalEnv)
	# rm(list = obj.name, envir = .GlobalEnv)
	# obj.name <- paste("mm9_", mm.chr, "_Red", prefix, ".allpred", sep = "")
	# file.name <- paste("mm9_", mm.chr, "_Red", prefix, "_allpred.RData", sep = "")
	# save(list = obj.name, file = file.name, envir = .GlobalEnv)
	# rm(list = obj.name, envir = .GlobalEnv)
	# setwd(wd)
	# 
	# # generate mm9_chr19_Uni.NuPoPMmMm and mm9_chr19_Red.NuPoPMmMm
	# getNearestDyadVitMm(prefix = prefix, nuc.size = nuc.size, mm.chr = mm.chr)
	# setwd(save.dir)
	# obj.name <- paste("mm9_", mm.chr, "_Uni", prefix, sep = "")
	# file.name <- paste("mm9_", mm.chr, "_Uni", prefix, ".RData", sep = "")
	# save(list = obj.name, file = file.name, envir = .GlobalEnv)
	# rm(list = obj.name, envir = .GlobalEnv)
	# obj.name <- paste("mm9_", mm.chr, "_Red", prefix, sep = "")
	# file.name <- paste("mm9_", mm.chr, "_Red", prefix, ".RData", sep = "")
	# save(list = obj.name, file = file.name, envir = .GlobalEnv)
	# rm(list = obj.name, envir = .GlobalEnv)
	# setwd(wd)
}

