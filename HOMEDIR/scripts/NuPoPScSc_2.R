setwd("scripts")
source("nearestDyadVit.R", chdir = TRUE)
source("sc_chemical.R", chdir = TRUE)

setwd("../")
setwd("sc_genome_2011")
setwd("predNuPoP_sc")

# non-redundant chemical dyad と Viterbi dyads を比べる
nature11142_s2.NuPoPScSc <- nearestDyadVit(dyad.table = nature11142_s2.sacCer3, chr.num = 16, 
			prefix = "NuPoPScSc", nuc.size = 147)

# redundant chemical dyad と Viterbi dyads を比べる
nature11142_s3.NuPoPScSc <- nearestDyadVit(dyad.table = nature11142_s3.sacCer3, chr.num = 16, 
			prefix = "NuPoPScSc", nuc.size = 147)

# 保存する。
setwd("../../")
setwd("RData")
save(nature11142_s2.NuPoPScSc, file = "nature11142_s2_sacCer3_NuPoPScSc.RData")
save(nature11142_s3.NuPoPScSc, file = "nature11142_s3_sacCer3_NuPoPScSc.RData")
