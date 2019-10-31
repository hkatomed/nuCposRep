setwd("scripts")
source("getDyadPred.R", chdir = TRUE)
source("nearestDyadPred.R", chdir = TRUE)
source("nearestDyadVit.R", chdir = TRUE)
source("sp_chemical.R", chdir = TRUE)

setwd("../")
setwd("sp_genome")
setwd("predNuCpos_sp")

# すべての predicted dyad を対象とする
sp_genome.nuCposSpSp.allpred <- getDyadPred(chr.num = 3, prefix = "nuCposSpSp", nuc.size = 147)

# non-redundant chemical dyad と比べる
sd01.nuCposSpSp.allpred <- nearestDyadPred(dyad.table = sd01, chr.num = 3, 
			prefix = "nuCposSpSp", nuc.size = 147)

# redundant chemical dyad と比べる
sd02.nuCposSpSp.allpred <- nearestDyadPred(dyad.table = sd02, chr.num = 3, 
			prefix = "nuCposSpSp", nuc.size = 147)

# 保存する。
setwd("../../")
dir.create("RData")
setwd("RData")
save(sp_genome.nuCposSpSp.allpred, file = "sp_genome_nuCposSpSp_allpred.RData")
save(sd01.nuCposSpSp.allpred, file = "sd01_nuCposSpSp_allpred.RData")
save(sd02.nuCposSpSp.allpred, file = "sd02_nuCposSpSp_allpred.RData")
