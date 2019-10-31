setwd("scripts")
source(file = "get.forROCR.R")

nature11142_s2.NuPoPScSc.allpred.forROCR <- get.forROCR(nature11142_s2.NuPoPScSc.allpred)
nature11142_s2.NuPoPScSp.allpred.forROCR <- get.forROCR(nature11142_s2.NuPoPScSp.allpred)
nature11142_s2.nuCposScSc.allpred.forROCR <- get.forROCR(nature11142_s2.nuCposScSc.allpred)
nature11142_s2.nuCposScSp.allpred.forROCR <- get.forROCR(nature11142_s2.nuCposScSp.allpred)

nature11142_s3.NuPoPScSc.allpred.forROCR <- get.forROCR(nature11142_s3.NuPoPScSc.allpred)
nature11142_s3.NuPoPScSp.allpred.forROCR <- get.forROCR(nature11142_s3.NuPoPScSp.allpred)
nature11142_s3.nuCposScSc.allpred.forROCR <- get.forROCR(nature11142_s3.nuCposScSc.allpred)
nature11142_s3.nuCposScSp.allpred.forROCR <- get.forROCR(nature11142_s3.nuCposScSp.allpred)

sd01.NuPoPSpSc.allpred.forROCR <- get.forROCR(sd01.NuPoPSpSc.allpred)
sd01.NuPoPSpSp.allpred.forROCR <- get.forROCR(sd01.NuPoPSpSp.allpred)
sd01.nuCposSpSc.allpred.forROCR <- get.forROCR(sd01.nuCposSpSc.allpred)
sd01.nuCposSpSp.allpred.forROCR <- get.forROCR(sd01.nuCposSpSp.allpred)

sd02.NuPoPSpSc.allpred.forROCR <- get.forROCR(sd02.NuPoPSpSc.allpred)
sd02.NuPoPSpSp.allpred.forROCR <- get.forROCR(sd02.NuPoPSpSp.allpred)
sd02.nuCposSpSc.allpred.forROCR <- get.forROCR(sd02.nuCposSpSc.allpred)
sd02.nuCposSpSp.allpred.forROCR <- get.forROCR(sd02.nuCposSpSp.allpred)

setwd("../")
setwd("ROCR_RData")
save(nature11142_s2.NuPoPScSc.allpred.forROCR, file = "nature11142_s2_NuPoPScSc_allpred_forROCR.RData")
save(nature11142_s2.NuPoPScSp.allpred.forROCR, file = "nature11142_s2_NuPoPScSp_allpred_forROCR.RData")
save(nature11142_s2.nuCposScSc.allpred.forROCR, file = "nature11142_s2_nuCposScSc_allpred_forROCR.RData")
save(nature11142_s2.nuCposScSp.allpred.forROCR, file = "nature11142_s2_nuCposScSp_allpred_forROCR.RData")

save(nature11142_s3.NuPoPScSc.allpred.forROCR, file = "nature11142_s3_NuPoPScSc_allpred_forROCR.RData")
save(nature11142_s3.NuPoPScSp.allpred.forROCR, file = "nature11142_s3_NuPoPScSp_allpred_forROCR.RData")
save(nature11142_s3.nuCposScSc.allpred.forROCR, file = "nature11142_s3_nuCposScSc_allpred_forROCR.RData")
save(nature11142_s3.nuCposScSp.allpred.forROCR, file = "nature11142_s3_nuCposScSp_allpred_forROCR.RData")

save(sd01.NuPoPSpSc.allpred.forROCR, file = "sd01_NuPoPSpSc_allpred_forROCR.RData")
save(sd01.NuPoPSpSp.allpred.forROCR, file = "sd01_NuPoPSpSp_allpred_forROCR.RData")
save(sd01.nuCposSpSc.allpred.forROCR, file = "sd01_nuCposSpSc_allpred_forROCR.RData")
save(sd01.nuCposSpSp.allpred.forROCR, file = "sd01_nuCposSpSp_allpred_forROCR.RData")

save(sd02.NuPoPSpSc.allpred.forROCR, file = "sd02_NuPoPSpSc_allpred_forROCR.RData")
save(sd02.NuPoPSpSp.allpred.forROCR, file = "sd02_NuPoPSpSp_allpred_forROCR.RData")
save(sd02.nuCposSpSc.allpred.forROCR, file = "sd02_nuCposSpSc_allpred_forROCR.RData")
save(sd02.nuCposSpSp.allpred.forROCR, file = "sd02_nuCposSpSp_allpred_forROCR.RData")
