
library(nuCpos)
library(rasterVis)

ori601 <- "CGGGATCCTAATGACCAAGGAAAGCATGATTCTTCACACCGAGTTCATCCCTTATGTGATGGACCCTATACGCGGCCGCCCTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCTGTGCATGTATTGAACAGCGACCTTGCCGGTGCCAGTCGGATAGTGTTCCGAGCTCCC"

ori601df <- data.frame(pos = 74:209, HBA = numeric(length = 136), stringsAsFactors = FALSE)
for(i in 1:136){
	seq <- substr(x = ori601, start = i, stop = i+146)
	ori601df$HBA[i] <- HBA(seq, species = "sc", silent = TRUE)
}

ori601df$lHBA_A <- as.numeric(NA)
ori601df$lHBA_B <- as.numeric(NA)
ori601df$lHBA_C <- as.numeric(NA)
ori601df$lHBA_D <- as.numeric(NA)
ori601df$lHBA_E <- as.numeric(NA)
ori601df$lHBA_F <- as.numeric(NA)
ori601df$lHBA_G <- as.numeric(NA)
ori601df$lHBA_H <- as.numeric(NA)
ori601df$lHBA_I <- as.numeric(NA)
ori601df$lHBA_J <- as.numeric(NA)
ori601df$lHBA_K <- as.numeric(NA)
ori601df$lHBA_L <- as.numeric(NA)
ori601df$lHBA_M <- as.numeric(NA)

for(i in 1:136){
	seq <- substr(x = ori601, start = i, stop = i+146)
	ori601df[i, 3:15] <- localHBA(seq, species = "sc", silent = TRUE)
}

## levelplot で描画する。
COL = BuRdTheme()$regions$col

grid <- as.matrix(ori601df[, 3:15])
grid <- grid[, ncol(grid):1]

gridNum <- 6
gridNew <- matrix(nrow = nrow(grid), ncol = ncol(grid)*gridNum)
for(i in 1:ncol(grid)) {
	row.range <- ((i-1)*gridNum + 1):((i-1)*gridNum + gridNum)
	gridNew[, row.range] <- grid[,i]
}

setwd("RData")
LIST <- c("COL", "gridNum", "gridNew", "ori601df")
save(list = LIST, file = "Widom_lHBA.RData")

testSeq <- "CTGGAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCCTGT"
seqName <- "Widom 601"
save(list = c("testSeq", "seqName"), file = "lHBA_Widom601.RData")

testSeq <- "ATCAATATCCACCTGCAGATACTACCAAAAGTGTATTTGGAAACTGCTCCATCAAAAGGCATGTTCAGCTGGAATCCAGCTGAACATGCCTTTTGATGGAGCAGTTTCCAAATACACTTTTGGTAGTATCTGCAGGTGGATATTGAT"
seqName <- "NCP147"
save(list = c("testSeq", "seqName"), file = "lHBA_NCP147.RData")

testSeq <- "AATCAGAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGCGCTGTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATACATCGATA"
seqName <- "Widom 601LR"
save(list = c("testSeq", "seqName"), file = "lHBA_Widom601LR.RData")

testSeq <- "AATCACAATCCCGGTGCCGAGGCCGCTCAATTGGTCGTAGACAGCTCTAGCACCGCTTAAACGCACGTACGGAATCCGTACGTGCGTTTAAGCGGTGCTAGAGCTGTCTACGACCAATTGAGCGGCCTCGGCACCGGGATTGTGATA"
seqName <- "Widom 601L"
save(list = c("testSeq", "seqName"), file = "lHBA_Widom601L.RData")

testSeq <- "AATCGATGAATATATCTGACACGTGCCTGGAGACTAGGGAGTAATCCCCTTGGCGGTTAAAACGCGGGGGAGAATCTCCCCCGCGTTTTAACCGCCAAGGGGATTACTCCCTAGTCTCCAGGCACGTGTCAGATATATTCATCGATA"
seqName <- "Widom 601R"
save(list = c("testSeq", "seqName"), file = "lHBA_Widom601R.RData")


