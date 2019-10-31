setwd("scripts")
source("nearestDyadVit.R", chdir = TRUE)
source("sc_chemical.R", chdir = TRUE)

setwd("../")
setwd("sc_genome_2011")
setwd("predNuCpos_sp")

# non-redundant chemical dyad と Viterbi dyads を比べる
nature11142_s2.nuCposScSp <- nearestDyadVit(dyad.table = nature11142_s2.sacCer3, chr.num = 16, 
			prefix = "nuCposScSp", nuc.size = 147)

# redundant chemical dyad と Viterbi dyads を比べる
nature11142_s3.nuCposScSp <- nearestDyadVit(dyad.table = nature11142_s3.sacCer3, chr.num = 16, 
			prefix = "nuCposScSp", nuc.size = 147)

# 保存する。
setwd("../../")
setwd("RData")
save(nature11142_s2.nuCposScSp, file = "nature11142_s2_sacCer3_nuCposScSp.RData")
save(nature11142_s3.nuCposScSp, file = "nature11142_s3_sacCer3_nuCposScSp.RData")
