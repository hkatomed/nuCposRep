setwd("scripts")
source("getDyadPred.R", chdir = TRUE)
source("nearestDyadPred.R", chdir = TRUE)
source("nearestDyadVit.R", chdir = TRUE)
source("sp_chemical.R", chdir = TRUE)

setwd("../")
setwd("sp_genome")
setwd("predNuCpos_sc")

# すべての predicted dyad を対象とする
sp_genome.nuCposSpSc.allpred <- getDyadPred(chr.num = 3, prefix = "nuCposSpSc", nuc.size = 147)

# non-redundant chemical dyad と比べる
sd01.nuCposSpSc.allpred <- nearestDyadPred(dyad.table = sd01, chr.num = 3, 
			prefix = "nuCposSpSc", nuc.size = 147)

# redundant chemical dyad と比べる
sd02.nuCposSpSc.allpred <- nearestDyadPred(dyad.table = sd02, chr.num = 3, 
			prefix = "nuCposSpSc", nuc.size = 147)

# 保存する。
setwd("../../")
dir.create("RData")
setwd("RData")
save(sp_genome.nuCposSpSc.allpred, file = "sp_genome_nuCposSpSc_allpred.RData")
save(sd01.nuCposSpSc.allpred, file = "sd01_nuCposSpSc_allpred.RData")
save(sd02.nuCposSpSc.allpred, file = "sd02_nuCposSpSc_allpred.RData")
