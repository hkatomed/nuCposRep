# 分裂酵母の chemical map データをロードする。参考：20151217_NuPoP_genome_1.txt
wd <- getwd()
# setwd("/Users/hrk_kato/Documents/Bioinfomatics_since20140304/Current/20151104_NuPoP/Bioconductor_org/20151218_NuPoP1")
# setwd("sp_chemical_map_wig")

setwd("../RData")

load(file = "sd01.RData")
load(file = "sd02.RData")
setwd(wd)
