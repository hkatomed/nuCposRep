setwd("scripts")
source("nearestDyadVit.R", chdir = TRUE)
source("sp_chemical.R", chdir = TRUE)

setwd("../")
setwd("sp_genome")
setwd("predNuCpos_sp")

# non-redundant chemical dyad と Viterbi dyads を比べる
sd01.nuCposSpSp <- nearestDyadVit(dyad.table = sd01, chr.num = 3, 
			prefix = "nuCposSpSp", nuc.size = 147)

# redundant chemical dyad と Viterbi dyads を比べる
sd02.nuCposSpSp <- nearestDyadVit(dyad.table = sd02, chr.num = 3, 
			prefix = "nuCposSpSp", nuc.size = 147)

# 保存する。
setwd("../../")
setwd("RData")
save(sd01.nuCposSpSp, file = "sd01_nuCposSpSp.RData")
save(sd02.nuCposSpSp, file = "sd02_nuCposSpSp.RData")
