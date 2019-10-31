setwd("src")
dyn.load("HBA_3.so")
setwd("../scripts")
source(file = "HBAm.R")
source(file = "HBAc.R")
setwd("../RData")
load(file = "sysdata_NuPoP.rda")
load(file = "sysdata_nuCpos.rda")

ori601 <- "CGGGATCCTAATGACCAAGGAAAGCATGATTCTTCACACCGAGTTCATCCCTTATGTGATGGACCCTATACGCGGCCGCCCTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCTGTGCATGTATTGAACAGCGACCTTGCCGGTGCCAGTCGGATAGTGTTCCGAGCTCCC"

ori601HBA.m <- data.frame(pos = 74:209, HBA = numeric(length = 136), stringsAsFactors = FALSE)
for(i in 1:136){
	seq <- substr(x = ori601, start = i, stop = i+146)
	ori601HBA.m$HBA[i] <- HBAm(seq, silent = TRUE)
}

ori601HBA.c <- data.frame(pos = 74:209, HBA = numeric(length = 136), stringsAsFactors = FALSE)
for(i in 1:136){
	seq <- substr(x = ori601, start = i, stop = i+146)
	ori601HBA.c$HBA[i] <- HBAc(seq, silent = TRUE)
}

ori601HBA.m.mean <- ori601HBA.m
ori601HBA.c.mean <- ori601HBA.c
ori601HBA.m.mean$HBA <- ori601HBA.m.mean$HBA - mean(ori601HBA.m.mean$HBA)
ori601HBA.c.mean$HBA <- ori601HBA.c.mean$HBA - mean(ori601HBA.c.mean$HBA)

save(ori601HBA.m.mean, file = "ori601HBA_m_mean.RData")
save(ori601HBA.c.mean, file = "ori601HBA_c_mean.RData")


MMTV <- "CAAAAACTTATGGCATGAGTTATTATGAATAGCCTTTATTGGCCCAACCTTGCGGTTCCCAGGGCTTAAGTAAGTTTTTGGTTACAAACTGTTCTTAAAACGAGGATGTGAGACAAGTGGTTTCCTGACTTGGTTTGGTATCAAAGGTTCTGATCTGAGCTCTGAGTGTTCTATTTTCCTATGTTCTTTTGGAATTTATCCAAATCTTATGTAAATGCTTATGTAAACCAAGATATAAAAGAGTGCTGATTTTTTGAGTAAACTTGCAACAGTCCTAACATTCACCTCTTGTGTGTTTGTGTCTGTTCGCCATCCCGTCTCCGCTCGTCACTTATCCTTCACTTTCCAGAGGGTCCCCCCGCAGACCCCGGCGACCCTCAGGTCGGCCGACTGCGGCACAGTTTTTTGCTCCTTTTTCTAGATGTAATTTTTAAAGCTTATTTTTTAACTTTCACATGTGCTACACTCACATGTGCAATGAGTGA"

DYAD <- c(139, 335)

MMTVHBA.m <- data.frame(pos = 74:(nchar(MMTV)-73), HBA = numeric(length = length(seq(74, nchar(MMTV)-73))), stringsAsFactors = FALSE)
for(i in 1:length(seq(74, nchar(MMTV)-73))){
	seq <- substr(x = MMTV, start = i, stop = i+146)
	MMTVHBA.m$HBA[i] <- HBAm(seq, silent = TRUE)
}

MMTVHBA.c <- data.frame(pos = 74:(nchar(MMTV)-73), HBA = numeric(length = length(seq(74, nchar(MMTV)-73))), stringsAsFactors = FALSE)
for(i in 1:length(seq(74, nchar(MMTV)-73))){
	seq <- substr(x = MMTV, start = i, stop = i+146)
	MMTVHBA.c$HBA[i] <- HBAc(seq, silent = TRUE)
}

MMTVHBA.m.mean <- MMTVHBA.m
MMTVHBA.c.mean <- MMTVHBA.c
MMTVHBA.m.mean$HBA <- MMTVHBA.m.mean$HBA - mean(MMTVHBA.m.mean$HBA)
MMTVHBA.c.mean$HBA <- MMTVHBA.c.mean$HBA - mean(MMTVHBA.c.mean$HBA)

save(MMTVHBA.m.mean, file = "MMTVHBA_m_mean.RData")
save(MMTVHBA.c.mean, file = "MMTVHBA_c_mean.RData")
save(MMTV, file = "MMTV.RData")
save(DYAD, file = "MMTV_DYAD.RData")


rDNA <- "GAGGAATTCCAACGAATAACTTCCAGGGATTTATAAGCCGATGACGTCATAACATCCCTGACCCTTTAAATAGCTTAACTTTCATCAAGCAAGAGCCTACGACCATACCATGCTGAATATACCGGTTCTCGTCCGATCACCGAAGTCAAGCAGCATAGGGCTCGGTTAGTACTTGGATGGGAGACCGCCTGGGAATACCGAATTCCCC"

DYAD <- 85

rDNAHBA.m <- data.frame(pos = 74:(nchar(rDNA)-73), HBA = numeric(length = length(seq(74, nchar(rDNA)-73))), stringsAsFactors = FALSE)
for(i in 1:length(seq(74, nchar(rDNA)-73))){
	seq <- substr(x = rDNA, start = i, stop = i+146)
	rDNAHBA.m$HBA[i] <- HBAm(seq, silent = TRUE)
}

rDNAHBA.c <- data.frame(pos = 74:(nchar(rDNA)-73), HBA = numeric(length = length(seq(74, nchar(rDNA)-73))), stringsAsFactors = FALSE)
for(i in 1:length(seq(74, nchar(rDNA)-73))){
	seq <- substr(x = rDNA, start = i, stop = i+146)
	rDNAHBA.c$HBA[i] <- HBAc(seq, silent = TRUE)
}

rDNAHBA.m.mean <- rDNAHBA.m
rDNAHBA.c.mean <- rDNAHBA.c
rDNAHBA.m.mean$HBA <- rDNAHBA.m.mean$HBA - mean(rDNAHBA.m.mean$HBA)
rDNAHBA.c.mean$HBA <- rDNAHBA.c.mean$HBA - mean(rDNAHBA.c.mean$HBA)

save(rDNAHBA.m.mean, file = "rDNAHBA_m_mean.RData")
save(rDNAHBA.c.mean, file = "rDNAHBA_c_mean.RData")
save(rDNA, file = "rDNA.RData")
save(DYAD, file = "rDNA_DYAD.RData")


