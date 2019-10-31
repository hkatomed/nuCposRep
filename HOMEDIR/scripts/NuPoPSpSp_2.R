setwd("scripts")
source("nearestDyadVit.R", chdir = TRUE)
source("sp_chemical.R", chdir = TRUE)

setwd("../")
setwd("sp_genome")
setwd("predNuPoP_sp")

# non-redundant chemical dyad と Viterbi dyads を比べる
sd01.NuPoPSpSp <- nearestDyadVit(dyad.table = sd01, chr.num = 3, 
			prefix = "NuPoPSpSp", nuc.size = 147)

# redundant chemical dyad と Viterbi dyads を比べる
sd02.NuPoPSpSp <- nearestDyadVit(dyad.table = sd02, chr.num = 3, 
			prefix = "NuPoPSpSp", nuc.size = 147)

# 保存する。
setwd("../../")
setwd("RData")
save(sd01.NuPoPSpSp, file = "sd01_NuPoPSpSp.RData")
save(sd02.NuPoPSpSp, file = "sd02_NuPoPSpSp.RData")
