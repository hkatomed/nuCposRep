# 出芽酵母 LiftOver 済みの chemical map データをロードする。参考：20151217_NuPoP_genome_1.txt
wd <- getwd()
# setwd("/Users/hrk_kato/Documents/Bioinfomatics_since20140304/Current/20151104_NuPoP/Bioconductor_org")
# setwd("20151218_NuPoP1/sc_chemical_map_liftOver")

setwd("../RData")

load(file = "nature11142_s2_sacCer3.RData")	#  67,548 non-redundant map
load(file = "nature11142_s3_sacCer3.RData")	# 344,709 redundant map (LiftOver で 11 dyad が脱落した）
setwd(wd)

## chemical の pos が character vector になっているので、integer にする。
nature11142_s2.sacCer3$pos <- as.integer(nature11142_s2.sacCer3$pos)
nature11142_s3.sacCer3$pos <- as.integer(nature11142_s3.sacCer3$pos)
