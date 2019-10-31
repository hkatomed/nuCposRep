setwd("scripts")
source("nearestDyadVit.R", chdir = TRUE)
source("sc_chemical.R", chdir = TRUE)

setwd("../")
setwd("sc_genome_2011")
setwd("predNuPoP_sp")

# non-redundant chemical dyad と Viterbi dyads を比べる
nature11142_s2.NuPoPScSp <- nearestDyadVit(dyad.table = nature11142_s2.sacCer3, chr.num = 16, 
			prefix = "NuPoPScSp", nuc.size = 147)

# redundant chemical dyad と Viterbi dyads を比べる
nature11142_s3.NuPoPScSp <- nearestDyadVit(dyad.table = nature11142_s3.sacCer3, chr.num = 16, 
			prefix = "NuPoPScSp", nuc.size = 147)

# 保存する。
setwd("../../")
setwd("RData")
save(nature11142_s2.NuPoPScSp, file = "nature11142_s2_sacCer3_NuPoPScSp.RData")
save(nature11142_s3.NuPoPScSp, file = "nature11142_s3_sacCer3_NuPoPScSp.RData")
