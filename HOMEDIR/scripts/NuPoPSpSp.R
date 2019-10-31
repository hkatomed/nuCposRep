setwd("scripts")
source("getDyadPred.R", chdir = TRUE)
source("nearestDyadPred.R", chdir = TRUE)
source("nearestDyadVit.R", chdir = TRUE)
source("sp_chemical.R", chdir = TRUE)

setwd("../")
setwd("sp_genome")
setwd("predNuPoP_sp")

# すべての predicted dyad を対象とする
sp_genome.NuPoPSpSp.allpred <- getDyadPred(chr.num = 3, prefix = "NuPoPSpSp", nuc.size = 147)

# non-redundant chemical dyad と比べる
sd01.NuPoPSpSp.allpred <- nearestDyadPred(dyad.table = sd01, chr.num = 3, 
			prefix = "NuPoPSpSp", nuc.size = 147)

# redundant chemical dyad と比べる
sd02.NuPoPSpSp.allpred <- nearestDyadPred(dyad.table = sd02, chr.num = 3, 
			prefix = "NuPoPSpSp", nuc.size = 147)

# 保存する。
setwd("../../")
dir.create("RData")
setwd("RData")
save(sp_genome.NuPoPSpSp.allpred, file = "sp_genome_NuPoPSpSp_allpred.RData")
save(sd01.NuPoPSpSp.allpred, file = "sd01_NuPoPSpSp_allpred.RData")
save(sd02.NuPoPSpSp.allpred, file = "sd02_NuPoPSpSp_allpred.RData")

