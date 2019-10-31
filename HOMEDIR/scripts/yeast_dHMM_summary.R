setwd("RData")

load(file = "sc_genome_2011_NuPoPScSc_allpred.RData")
load(file = "nature11142_s2_sacCer3_NuPoPScSc_allpred.RData")
load(file = "nature11142_s3_sacCer3_NuPoPScSc_allpred.RData")
load(file = "nature11142_s2_sacCer3_NuPoPScSc.RData")
load(file = "nature11142_s3_sacCer3_NuPoPScSc.RData")

load(file = "sc_genome_2011_NuPoPScSp_allpred.RData")
load(file = "nature11142_s2_sacCer3_NuPoPScSp_allpred.RData")
load(file = "nature11142_s3_sacCer3_NuPoPScSp_allpred.RData")
load(file = "nature11142_s2_sacCer3_NuPoPScSp.RData")
load(file = "nature11142_s3_sacCer3_NuPoPScSp.RData")

load(file = "sc_genome_2011_nuCposScSc_allpred.RData")
load(file = "nature11142_s2_sacCer3_nuCposScSc_allpred.RData")
load(file = "nature11142_s3_sacCer3_nuCposScSc_allpred.RData")
load(file = "nature11142_s2_sacCer3_nuCposScSc.RData")
load(file = "nature11142_s3_sacCer3_nuCposScSc.RData")

load(file = "sc_genome_2011_nuCposScSp_allpred.RData")
load(file = "nature11142_s2_sacCer3_nuCposScSp_allpred.RData")
load(file = "nature11142_s3_sacCer3_nuCposScSp_allpred.RData")
load(file = "nature11142_s2_sacCer3_nuCposScSp.RData")
load(file = "nature11142_s3_sacCer3_nuCposScSp.RData")

load(file = "sp_genome_NuPoPSpSc_allpred.RData")
load(file = "sd01_NuPoPSpSc_allpred.RData")
load(file = "sd02_NuPoPSpSc_allpred.RData")
load(file = "sd01_NuPoPSpSc.RData")
load(file = "sd02_NuPoPSpSc.RData")

load(file = "sp_genome_NuPoPSpSp_allpred.RData")
load(file = "sd01_NuPoPSpSp_allpred.RData")
load(file = "sd02_NuPoPSpSp_allpred.RData")
load(file = "sd01_NuPoPSpSp.RData")
load(file = "sd02_NuPoPSpSp.RData")

load(file = "sp_genome_nuCposSpSc_allpred.RData")
load(file = "sd01_nuCposSpSc_allpred.RData")
load(file = "sd02_nuCposSpSc_allpred.RData")
load(file = "sd01_nuCposSpSc.RData")
load(file = "sd02_nuCposSpSc.RData")

load(file = "sp_genome_nuCposSpSp_allpred.RData")
load(file = "sd01_nuCposSpSp_allpred.RData")
load(file = "sd02_nuCposSpSp_allpred.RData")
load(file = "sd01_nuCposSpSp.RData")
load(file = "sd02_nuCposSpSp.RData")


num.Summary <- data.frame(species = c("sc", "sc", "sc", "sc", "sp", "sp", "sp", "sp"), 
			model = c("NuPoP_sc", "NuPoP_sp", "nuCpos_sc", "nuCpos_sp", 
				"NuPoP_sc", "NuPoP_sp", "nuCpos_sc", "nuCpos_sp"), 
			AllPred = integer(8), Viterbi = integer(8), 
			NonRedVit = integer(8), RedVit = integer(8))

num.Summary$AllPred[1] <- nrow(sc_genome_2011.NuPoPScSc.allpred)
num.Summary$AllPred[2] <- nrow(sc_genome_2011.NuPoPScSp.allpred)
num.Summary$AllPred[3] <- nrow(sc_genome_2011.nuCposScSc.allpred)
num.Summary$AllPred[4] <- nrow(sc_genome_2011.nuCposScSp.allpred)
num.Summary$AllPred[5] <- nrow(sp_genome.NuPoPSpSc.allpred)
num.Summary$AllPred[6] <- nrow(sp_genome.NuPoPSpSp.allpred)
num.Summary$AllPred[7] <- nrow(sp_genome.nuCposSpSc.allpred)
num.Summary$AllPred[8] <- nrow(sp_genome.nuCposSpSp.allpred)

num.Summary$Viterbi[1] <- nrow(subset(sc_genome_2011.NuPoPScSc.allpred, Viterbi == 1))
num.Summary$Viterbi[2] <- nrow(subset(sc_genome_2011.NuPoPScSp.allpred, Viterbi == 1))
num.Summary$Viterbi[3] <- nrow(subset(sc_genome_2011.nuCposScSc.allpred, Viterbi == 1))
num.Summary$Viterbi[4] <- nrow(subset(sc_genome_2011.nuCposScSp.allpred, Viterbi == 1))
num.Summary$Viterbi[5] <- nrow(subset(sp_genome.NuPoPSpSc.allpred, Viterbi == 1))
num.Summary$Viterbi[6] <- nrow(subset(sp_genome.NuPoPSpSp.allpred, Viterbi == 1))
num.Summary$Viterbi[7] <- nrow(subset(sp_genome.nuCposSpSc.allpred, Viterbi == 1))
num.Summary$Viterbi[8] <- nrow(subset(sp_genome.nuCposSpSp.allpred, Viterbi == 1))

num.Summary$NonRedVit[1] <- nrow(subset(nature11142_s2.NuPoPScSc[["PredVsChem"]], abs(dist) < 2))
num.Summary$NonRedVit[2] <- nrow(subset(nature11142_s2.NuPoPScSp[["PredVsChem"]], abs(dist) < 2))
num.Summary$NonRedVit[3] <- nrow(subset(nature11142_s2.nuCposScSc[["PredVsChem"]], abs(dist) < 2))
num.Summary$NonRedVit[4] <- nrow(subset(nature11142_s2.nuCposScSp[["PredVsChem"]], abs(dist) < 2))
num.Summary$NonRedVit[5] <- nrow(subset(sd01.NuPoPSpSc[["PredVsChem"]], abs(dist) < 2))
num.Summary$NonRedVit[6] <- nrow(subset(sd01.NuPoPSpSp[["PredVsChem"]], abs(dist) < 2))
num.Summary$NonRedVit[7] <- nrow(subset(sd01.nuCposSpSc[["PredVsChem"]], abs(dist) < 2))
num.Summary$NonRedVit[8] <- nrow(subset(sd01.nuCposSpSp[["PredVsChem"]], abs(dist) < 2))


num.Summary$RedVit[1] <- nrow(subset(nature11142_s3.NuPoPScSc[["PredVsChem"]], abs(dist) < 2))
num.Summary$RedVit[2] <- nrow(subset(nature11142_s3.NuPoPScSp[["PredVsChem"]], abs(dist) < 2))
num.Summary$RedVit[3] <- nrow(subset(nature11142_s3.nuCposScSc[["PredVsChem"]], abs(dist) < 2))
num.Summary$RedVit[4] <- nrow(subset(nature11142_s3.nuCposScSp[["PredVsChem"]], abs(dist) < 2))
num.Summary$RedVit[5] <- nrow(subset(sd02.NuPoPSpSc[["PredVsChem"]], abs(dist) < 2))
num.Summary$RedVit[6] <- nrow(subset(sd02.NuPoPSpSp[["PredVsChem"]], abs(dist) < 2))
num.Summary$RedVit[7] <- nrow(subset(sd02.nuCposSpSc[["PredVsChem"]], abs(dist) < 2))
num.Summary$RedVit[8] <- nrow(subset(sd02.nuCposSpSp[["PredVsChem"]], abs(dist) < 2))

num.Summary$RedOnlyVit <- num.Summary$RedVit - num.Summary$NonRedVit
num.Summary$Others <- num.Summary$Viterbi - num.Summary$NonRedVit - num.Summary$RedOnlyVit

num.Summary$NonRedVitRatio <- num.Summary$NonRedVit/num.Summary$Viterbi
num.Summary$RedVitRatio <- num.Summary$RedVit/num.Summary$Viterbi
num.Summary$RedOnlyVitRatio <- num.Summary$RedOnlyVit/num.Summary$Viterbi
num.Summary$OthersRatio <- num.Summary$Others/num.Summary$Viterbi

setwd("../summary_RData")
save(num.Summary, file = "num_Summary.RData")
