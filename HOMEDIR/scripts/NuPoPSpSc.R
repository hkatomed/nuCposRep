setwd("scripts")
source("getDyadPred.R", chdir = TRUE)
source("nearestDyadPred.R", chdir = TRUE)
source("nearestDyadVit.R", chdir = TRUE)
source("sp_chemical.R", chdir = TRUE)

setwd("../")
setwd("sp_genome")
setwd("predNuPoP_sc")

# すべての predicted dyad を対象とする
sp_genome.NuPoPSpSc.allpred <- getDyadPred(chr.num = 3, prefix = "NuPoPSpSc", nuc.size = 147)

# non-redundant chemical dyad と比べる
sd01.NuPoPSpSc.allpred <- nearestDyadPred(dyad.table = sd01, chr.num = 3, 
			prefix = "NuPoPSpSc", nuc.size = 147)

# redundant chemical dyad と比べる
sd02.NuPoPSpSc.allpred <- nearestDyadPred(dyad.table = sd02, chr.num = 3, 
			prefix = "NuPoPSpSc", nuc.size = 147)

# 保存する。
setwd("../../")
dir.create("RData")
setwd("RData")
save(sp_genome.NuPoPSpSc.allpred, file = "sp_genome_NuPoPSpSc_allpred.RData")
save(sd01.NuPoPSpSc.allpred, file = "sd01_NuPoPSpSc_allpred.RData")
save(sd02.NuPoPSpSc.allpred, file = "sd02_NuPoPSpSc_allpred.RData")

