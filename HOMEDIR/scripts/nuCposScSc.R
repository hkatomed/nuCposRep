setwd("scripts")
source("getDyadPred.R", chdir = TRUE)
source("nearestDyadPred.R", chdir = TRUE)
source("nearestDyadVit.R", chdir = TRUE)
source("sc_chemical.R", chdir = TRUE)

setwd("../")
setwd("sc_genome_2011")
setwd("predNuCpos_sc")

# すべての predicted dyad を対象とする
sc_genome_2011.nuCposScSc.allpred <- getDyadPred(chr.num = 16, prefix = "nuCposScSc", nuc.size = 147)

# non-redundant chemical dyad と比べる
nature11142_s2.nuCposScSc.allpred <- nearestDyadPred(dyad.table = nature11142_s2.sacCer3, chr.num = 16, 
			prefix = "nuCposScSc", nuc.size = 147)

# redundant chemical dyad と比べる
nature11142_s3.nuCposScSc.allpred <- nearestDyadPred(dyad.table = nature11142_s3.sacCer3, chr.num = 16, 
			prefix = "nuCposScSc", nuc.size = 147)

# 保存する。
setwd("../../")
dir.create("RData")
setwd("RData")
save(sc_genome_2011.nuCposScSc.allpred, file = "sc_genome_2011_nuCposScSc_allpred.RData")
save(nature11142_s2.nuCposScSc.allpred, file = "nature11142_s2_sacCer3_nuCposScSc_allpred.RData")
save(nature11142_s3.nuCposScSc.allpred, file = "nature11142_s3_sacCer3_nuCposScSc_allpred.RData")
