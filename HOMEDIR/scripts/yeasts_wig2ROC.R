setwd("scripts")
source("getPredTable.R", chdir = TRUE)
source("getNearestDyad.R", chdir = TRUE)
source("nearestDyadPred.R", chdir = TRUE)
source("getNearestDyadVit.R", chdir = TRUE)
source("NearestDyadVit.R", chdir = TRUE)
source("matchChem.R", chdir = TRUE)
source("MultiMatchChem.R", chdir = TRUE)
source("pred.summary.R", chdir = TRUE)
source("sc_chemical.R", chdir = TRUE)
source("sp_chemical.R", chdir = TRUE)
source("wig2ROC.R", chdir = TRUE)

setwd("../")
setwd("sc_genome_2011")
setwd("predNuPoP_sc")
prefix <- "NuPoPScSc"
wig2ROC(prefix = prefix, nuc.size = 147)

setwd("../")
setwd("predNuPoP_sp")
prefix <- "NuPoPScSp"
wig2ROC(prefix = prefix, nuc.size = 147)

setwd("../")
setwd("predNuCpos_sc")
prefix <- "nuCposScSc"
wig2ROC(prefix = prefix, nuc.size = 147)

setwd("../")
setwd("predNuCpos_sp")
prefix <- "nuCposScSp"
wig2ROC(prefix = prefix, nuc.size = 147)

setwd("../../")
setwd("sp_genome")
setwd("predNuPoP_sc")
prefix <- "NuPoPSpSc"
wig2ROC(prefix = prefix, nuc.size = 147)

setwd("../")
setwd("predNuPoP_sp")
prefix <- "NuPoPSpSp"
wig2ROC(prefix = prefix, nuc.size = 147)

setwd("../")
setwd("predNuCpos_sc")
prefix <- "nuCposSpSc"
wig2ROC(prefix = prefix, nuc.size = 147)

setwd("../")
setwd("predNuCpos_sp")
prefix <- "nuCposSpSp"
wig2ROC(prefix = prefix, nuc.size = 147)

