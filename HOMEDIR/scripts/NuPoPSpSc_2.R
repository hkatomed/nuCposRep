setwd("scripts")
source("nearestDyadVit.R", chdir = TRUE)
source("sp_chemical.R", chdir = TRUE)

setwd("../")
setwd("sp_genome")
setwd("predNuPoP_sc")

# non-redundant chemical dyad と Viterbi dyads を比べる
sd01.NuPoPSpSc <- nearestDyadVit(dyad.table = sd01, chr.num = 3, 
			prefix = "NuPoPSpSc", nuc.size = 147)

# redundant chemical dyad と Viterbi dyads を比べる
sd02.NuPoPSpSc <- nearestDyadVit(dyad.table = sd02, chr.num = 3, 
			prefix = "NuPoPSpSc", nuc.size = 147)

# 保存する。
setwd("../../")
setwd("RData")
save(sd01.NuPoPSpSc, file = "sd01_NuPoPSpSc.RData")
save(sd02.NuPoPSpSc, file = "sd02_NuPoPSpSc.RData")
