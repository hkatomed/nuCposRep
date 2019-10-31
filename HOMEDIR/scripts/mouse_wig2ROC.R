setwd("scripts")
source("getPredTableMm.R", chdir = TRUE)
source("getNearestDyadMm.R", chdir = TRUE)
source("nearestDyadPredMm.R", chdir = TRUE)
source("getNearestDyadVitMm.R", chdir = TRUE)
source("NearestDyadVitMm.R", chdir = TRUE)
source("matchChemMm.R", chdir = TRUE)
source("MultiMatchChemMm.R", chdir = TRUE)
source("pred.summaryMm.R", chdir = TRUE)
source("mm_chemical.R", chdir = TRUE)
source("wig2ROCMm.R", chdir = TRUE)

setwd("../")
setwd("mm_genome")
setwd("predNuPoP_mm")
prefix <- "NuPoPMmMm"
wig2ROCMm(prefix = prefix, nuc.size = 147, mm.chr = "chr19")

setwd("../")
setwd("predNuCpos_mm")
prefix <- "nuCposMmMm"
wig2ROCMm(prefix = prefix, nuc.size = 147, mm.chr = "chr19")

