
library(Biostrings)
library(nuCpos)

setwd("RData")
load(file = "sp_genome.RData")

dyads <- c(17, 164, 332, 470, 629, 806)
genomic.dyads <- 116575 - (dyads-1)
nuc.name <- c("+1", "+2", "+3", "+4", "+5", "+6")
ura4.nucleosomes <- data.frame(nuc.name, chr = "chr3", 
	pos = genomic.dyads, strand = "-", stringsAsFactors = FALSE)

ura4.nucleosomes$seq <- ""
for(i in 1:6){
	pos <- ura4.nucleosomes$pos[i]
	seq <- reverseComplement(sp.genome[["chr3"]][(pos-73):(pos+73)])
	ura4.nucleosomes$seq[i] <- as.character(seq)
}

seqUra4_0.5k <- reverseComplement(sp.genome[["chr3"]][(115781-500):(116575+500)])
seqUra4_0.5k <- as.character(seqUra4_0.5k)

seqUra4_0.5kdf <- data.frame(pos = 74:1722, HBA = numeric(length = 1649), stringsAsFactors = FALSE)
for(i in 1:1649){
	seq <- substr(x = seqUra4_0.5k, start = i, stop = i+146)
	seqUra4_0.5kdf$HBA[i] <- HBA(seq, species = "sp", silent = TRUE)
}

seqUra4_0.5kdf$lHBA_A <- as.numeric(NA)
seqUra4_0.5kdf$lHBA_B <- as.numeric(NA)
seqUra4_0.5kdf$lHBA_C <- as.numeric(NA)
seqUra4_0.5kdf$lHBA_D <- as.numeric(NA)
seqUra4_0.5kdf$lHBA_E <- as.numeric(NA)
seqUra4_0.5kdf$lHBA_F <- as.numeric(NA)
seqUra4_0.5kdf$lHBA_G <- as.numeric(NA)
seqUra4_0.5kdf$lHBA_H <- as.numeric(NA)
seqUra4_0.5kdf$lHBA_I <- as.numeric(NA)
seqUra4_0.5kdf$lHBA_J <- as.numeric(NA)
seqUra4_0.5kdf$lHBA_K <- as.numeric(NA)
seqUra4_0.5kdf$lHBA_L <- as.numeric(NA)
seqUra4_0.5kdf$lHBA_M <- as.numeric(NA)

for(i in 1:1649){
	seq <- substr(x = seqUra4_0.5k, start = i, stop = i+146)
	seqUra4_0.5kdf[i, 3:15] <- localHBA(seq, species = "sp", silent = TRUE)
}

grid <- as.matrix(seqUra4_0.5kdf[, 3:15])
grid <- grid[, ncol(grid):1]

gridNum <- 48
gridNew <- matrix(nrow = nrow(grid), ncol = ncol(grid)*gridNum)
for(i in 1:ncol(grid)) {
	row.range <- ((i-1)*gridNum + 1):((i-1)*gridNum + gridNum)
	gridNew[, row.range] <- grid[,i]
}

ura4.nucleosomes$HBA <- as.numeric(NA)
for(i in 1:6){
	seq <- ura4.nucleosomes$seq[i]
	ura4.nucleosomes$HBA[i] <- HBA(seq, species = "sp", silent = TRUE)
}

ura4.nucleosomes$lHBA_A <- as.numeric(NA)
ura4.nucleosomes$lHBA_B <- as.numeric(NA)
ura4.nucleosomes$lHBA_C <- as.numeric(NA)
ura4.nucleosomes$lHBA_D <- as.numeric(NA)
ura4.nucleosomes$lHBA_E <- as.numeric(NA)
ura4.nucleosomes$lHBA_F <- as.numeric(NA)
ura4.nucleosomes$lHBA_G <- as.numeric(NA)
ura4.nucleosomes$lHBA_H <- as.numeric(NA)
ura4.nucleosomes$lHBA_I <- as.numeric(NA)
ura4.nucleosomes$lHBA_J <- as.numeric(NA)
ura4.nucleosomes$lHBA_K <- as.numeric(NA)
ura4.nucleosomes$lHBA_L <- as.numeric(NA)
ura4.nucleosomes$lHBA_M <- as.numeric(NA)

for(i in 1:6){
	seq <- ura4.nucleosomes$seq[i]
	ura4.nucleosomes[i, 7:19] <- localHBA(seq, species = "sp", silent = TRUE)
}

seqUra4_0.5k <- reverseComplement(sp.genome[["chr3"]][(115781-500):(116575+500)])
ura4WT <- substring(seqUra4_0.5k, first = 501, last = 500+795)

ura4Dyad <- "ATGGACGCGCGCGTGTTCCAGTCATACAGCGCGCGCGCGGAAGGGATGAAAAATCCCATTGCCAAGGAATTGTTGGCTTTGATGGAAGAAAAGCAAAGCAACTTGTCAGTCGCGGTCGATTTGACGAAGAAATCCGAAATACTCGAGCTAGTCGACAAGATCGGGCCATACGTATGCGTCATTAAAACGCATATTGACGTTGTCGAGGATTTCGACCAGGATATGGTAGAAAAACTGGTGGCCTTAGGTAAAAAGCATCGTTTTCTTATCTTTGAGGATCGCAAATTCGCAGACATTGGAAATACCGTTAAATTGCAGTACGCGAGCGGCGTATATAAGATCGCGAGCTGGGCGCACATCACAAATTGCCATACAGTGCCAGGCGAGGGTATTATACAAGGCCTCAAAGAAGTTGGTTTACCTTTGGGACGTGGTCTCTTGCTTCTAGCAGAGATGAGCAGTAAGGGAAGCCTAGCGACGGGCAGTTATACGGAAAAAACCTTAGAATGGTTTGAGAAGCATACCGATTTTTGCTTTGGCTTTATAGCTGGTCGTCGATTTCCTAACCTTCAAAGCGACTACATAACTATGTCCCCTGGTATCGGACTAGACGTCAAGGGCGATGGATTAGGCCAACAGTACAGGACGCCGGAGGAAGTGATTGTAAACTGCGGTAGCGATATCATCATTGTTGGTCGTGGAGTCTATGGAGCTGGTCGTAATCCTGTTGTCGAAGCCAAGAGATATAGAGAAGCTGGTTGGAAGGCATATCAGCAAAGATTGAGCCAACACTAA"

ura4Linker <- "ATGGATGCTAGAGTATTTCAAAGCTATTCAGCTAGAGCTGAGGGGATGAAAAATCCCATTGCCAAGGAATTGTTGGCGCTAATGGAGGAGAAACAGTCAAATCTAAGCGTAGCTGTAGACCTAACGAAGAAATCCGAAATCTTAGAATTGGTAGATAAAATTGGACCCTATGTCTGTGTTATCAAGACACATATTGACGTTGTCGAGGATTTCGACCAGGATATGGTAGAGAAGTTAGTAGCACTCGGCAAGAAACACAGGTTCTTGATATTCGAAGATCGCAAATTCGCAGACATTGGAAATACCGTCAAGCTACAATATGCATCTGGTGTGTACAAAATTGCTTCTTGGGCTCATATCACAAATTGCCATACAGTACCGGGTGAAGGCATCATCCAGGGATTAAAGGAGGTCGGCCTGCCGTTGGGACGTGGTCTCTTGCTTTTGGCTGAAATGTCTTCCAAAGGCTCTTTGGCTACTGGTTCCTACACAGAGAAAACCTTAGAATGGTTTGAGAAGCATACCGATTTTTGCTTTGGCTTTATCGCGGGCAGGAGGTTCCCGAATTTGCAGTCAGATTATATCACGATGAGCCCTGGTATCGGCTTGGATGTTAAAGGAGACGGGCTGGGACAGCAATATCGTACTCCTGAAGAAGTGATTGTAAACTGCGGTAGCGATATCATCATTGTTGGTCGTGGAGTCTATGGAGCCGGCAGGAACCCGGTCGTAGAGGCTAAACGCTACCGCGAGGCCGGCTGGAAGGCATATCAGCAAAGACTTTCTCAGCATTAA"

ura4Int <- "ATGGATGCTAGAGTATTTCAAAGCTATTCAGCTAGAGCTGAGGGGATGAAGAACCCGATCGCTAAAGAGCTATTGGCTTTGATGGAAGAAAAGCAAAGCAACTTGTCAGTCGCGGTCGATCTAACTAAAAAGAGTGAGATACTCGAATTGGTAGATAAAATTGGACCCTATGTCTGTGTTATCAAGACACATATTGACGTCGTAGAAGACTTTGATCAAGACATGGTAGAAAAACTGGTGGCCTTAGGTAAAAAGCATCGTTTTCTTATCTTTGAGGATAGAAAGTTTGCGGATATCGGCAACACCGTCAAGCTACAATATGCATCTGGTGTGTACAAAATTGCTTCTTGGGCTCACATAACGAACTGTCACACCGTACCAGGCGAGGGTATTATACAAGGCCTCAAAGAAGTTGGTTTACCTCTAGGCAGGGGCTTACTATTGCTAGCTGAAATGTCTTCCAAAGGCTCTTTGGCTACTGGTTCCTACACAGAGAAAACCTTAGAATGGTTCGAAAAACACACTGACTTCTGCTTTGGCTTTATAGCTGGTCGTCGATTTCCTAACCTTCAAAGCGACTACATAACTATGAGTCCGGGCATAGGACTAGACGTTAAAGGAGACGGGCTGGGACAGCAATATCGTACTCCTGAAGAAGTGATTGTAAACTGTGGCTCAGACATTATTATCGTCGGTCGTGGAGTCTATGGAGCTGGTCGTAATCCTGTTGTCGAAGCCAAGAGATATAGAGAAGCTGGCTGGAAAGCGTACCAACAGCGCCTTTCTCAGCATTAA"

seqUra4_0.5k.WT <- seqUra4_0.5k
seqUra4_0.5k.Dyad <- sub(pattern = ura4WT, 
		replacement = ura4Dyad, seqUra4_0.5k)
seqUra4_0.5k.Linker <- sub(pattern = ura4WT, 
		replacement = ura4Linker, seqUra4_0.5k)
seqUra4_0.5k.Int <- sub(pattern = ura4WT, 
		replacement = ura4Int, seqUra4_0.5k)

save(list = c("seqUra4_0.5k.WT", "seqUra4_0.5k.Dyad", "seqUra4_0.5k.Linker",
	"seqUra4_0.5k.Int"), file = "ura4_seqs.RData")

for(i in 1:6){
	pos <- dyads[i]
	seq <- substring(seqUra4_0.5k.Dyad, first = 500+pos-73, last = 500+pos+73)
	ura4.nucleosomes$seqDyad[i] <- seq
}

for(i in 1:6){
	pos <- dyads[i]
	seq <- substring(seqUra4_0.5k.Linker, first = 500+pos-73, last = 500+pos+73)
	ura4.nucleosomes$seqLinker[i] <- seq
}

for(i in 1:6){
	pos <- dyads[i]
	seq <- substring(seqUra4_0.5k.Int, first = 500+pos-73, last = 500+pos+73)
	ura4.nucleosomes$seqInt[i] <- seq
}

ura4.nucleosomes$HBA.Dyad <- as.numeric(NA)
for(i in 1:6){
	seq <- ura4.nucleosomes$seqDyad[i]
	ura4.nucleosomes$HBA.Dyad[i] <- HBA(seq, species = "sp", silent = TRUE)
}

ura4.nucleosomes$lHBA_A.Dyad <- as.numeric(NA)
ura4.nucleosomes$lHBA_B.Dyad <- as.numeric(NA)
ura4.nucleosomes$lHBA_C.Dyad <- as.numeric(NA)
ura4.nucleosomes$lHBA_D.Dyad <- as.numeric(NA)
ura4.nucleosomes$lHBA_E.Dyad <- as.numeric(NA)
ura4.nucleosomes$lHBA_F.Dyad <- as.numeric(NA)
ura4.nucleosomes$lHBA_G.Dyad <- as.numeric(NA)
ura4.nucleosomes$lHBA_H.Dyad <- as.numeric(NA)
ura4.nucleosomes$lHBA_I.Dyad <- as.numeric(NA)
ura4.nucleosomes$lHBA_J.Dyad <- as.numeric(NA)
ura4.nucleosomes$lHBA_K.Dyad <- as.numeric(NA)
ura4.nucleosomes$lHBA_L.Dyad <- as.numeric(NA)
ura4.nucleosomes$lHBA_M.Dyad <- as.numeric(NA)

for(i in 1:6){
	seq <- ura4.nucleosomes$seqDyad[i]
	ura4.nucleosomes[i, 24:36] <- localHBA(seq, species = "sp", silent = TRUE)
}

ura4.nucleosomes$HBA.Linker <- as.numeric(NA)
for(i in 1:6){
	seq <- ura4.nucleosomes$seqLinker[i]
	ura4.nucleosomes$HBA.Linker[i] <- HBA(seq, species = "sp", silent = TRUE)
}

ura4.nucleosomes$lHBA_A.Linker <- as.numeric(NA)
ura4.nucleosomes$lHBA_B.Linker <- as.numeric(NA)
ura4.nucleosomes$lHBA_C.Linker <- as.numeric(NA)
ura4.nucleosomes$lHBA_D.Linker <- as.numeric(NA)
ura4.nucleosomes$lHBA_E.Linker <- as.numeric(NA)
ura4.nucleosomes$lHBA_F.Linker <- as.numeric(NA)
ura4.nucleosomes$lHBA_G.Linker <- as.numeric(NA)
ura4.nucleosomes$lHBA_H.Linker <- as.numeric(NA)
ura4.nucleosomes$lHBA_I.Linker <- as.numeric(NA)
ura4.nucleosomes$lHBA_J.Linker <- as.numeric(NA)
ura4.nucleosomes$lHBA_K.Linker <- as.numeric(NA)
ura4.nucleosomes$lHBA_L.Linker <- as.numeric(NA)
ura4.nucleosomes$lHBA_M.Linker <- as.numeric(NA)

for(i in 1:6){
	seq <- ura4.nucleosomes$seqLinker[i]
	ura4.nucleosomes[i, 38:50] <- localHBA(seq, species = "sp", silent = TRUE)
}

ura4.nucleosomes$HBA.Int <- as.numeric(NA)
for(i in 1:6){
	seq <- ura4.nucleosomes$seqInt[i]
	ura4.nucleosomes$HBA.Int[i] <- HBA(seq, species = "sp", silent = TRUE)
}

ura4.nucleosomes$lHBA_A.Int <- as.numeric(NA)
ura4.nucleosomes$lHBA_B.Int <- as.numeric(NA)
ura4.nucleosomes$lHBA_C.Int <- as.numeric(NA)
ura4.nucleosomes$lHBA_D.Int <- as.numeric(NA)
ura4.nucleosomes$lHBA_E.Int <- as.numeric(NA)
ura4.nucleosomes$lHBA_F.Int <- as.numeric(NA)
ura4.nucleosomes$lHBA_G.Int <- as.numeric(NA)
ura4.nucleosomes$lHBA_H.Int <- as.numeric(NA)
ura4.nucleosomes$lHBA_I.Int <- as.numeric(NA)
ura4.nucleosomes$lHBA_J.Int <- as.numeric(NA)
ura4.nucleosomes$lHBA_K.Int <- as.numeric(NA)
ura4.nucleosomes$lHBA_L.Int <- as.numeric(NA)
ura4.nucleosomes$lHBA_M.Int <- as.numeric(NA)

for(i in 1:6){
	seq <- ura4.nucleosomes$seqInt[i]
	ura4.nucleosomes[i, 52:64] <- localHBA(seq, species = "sp", silent = TRUE)
}

save(list = c("dyads", "grid", "gridNum", "gridNew", 
	"ura4.nucleosomes", "seqUra4_0.5kdf"), file = "lHBA_ura4_WT.RData")


seqUra4_0.5k.Dyaddf <- data.frame(pos = 74:1722, HBA = numeric(length = 1649), stringsAsFactors = FALSE)
for(i in 1:1649){
	seq <- substr(x = seqUra4_0.5k.Dyad, start = i, stop = i+146)
	seqUra4_0.5k.Dyaddf$HBA[i] <- HBA(seq, species = "sp", silent = TRUE)
}

seqUra4_0.5k.Dyaddf$lHBA_A <- as.numeric(NA)
seqUra4_0.5k.Dyaddf$lHBA_B <- as.numeric(NA)
seqUra4_0.5k.Dyaddf$lHBA_C <- as.numeric(NA)
seqUra4_0.5k.Dyaddf$lHBA_D <- as.numeric(NA)
seqUra4_0.5k.Dyaddf$lHBA_E <- as.numeric(NA)
seqUra4_0.5k.Dyaddf$lHBA_F <- as.numeric(NA)
seqUra4_0.5k.Dyaddf$lHBA_G <- as.numeric(NA)
seqUra4_0.5k.Dyaddf$lHBA_H <- as.numeric(NA)
seqUra4_0.5k.Dyaddf$lHBA_I <- as.numeric(NA)
seqUra4_0.5k.Dyaddf$lHBA_J <- as.numeric(NA)
seqUra4_0.5k.Dyaddf$lHBA_K <- as.numeric(NA)
seqUra4_0.5k.Dyaddf$lHBA_L <- as.numeric(NA)
seqUra4_0.5k.Dyaddf$lHBA_M <- as.numeric(NA)

for(i in 1:1649){
	seq <- substr(x = seqUra4_0.5k.Dyad, start = i, stop = i+146)
	seqUra4_0.5k.Dyaddf[i, 3:15] <- localHBA(seq, species = "sp", silent = TRUE)
}

grid <- as.matrix(seqUra4_0.5k.Dyaddf[, 3:15])
grid <- grid[, ncol(grid):1]

gridNum <- 48
gridNew <- matrix(nrow = nrow(grid), ncol = ncol(grid)*gridNum)
for(i in 1:ncol(grid)) {
	row.range <- ((i-1)*gridNum + 1):((i-1)*gridNum + gridNum)
	gridNew[, row.range] <- grid[,i]
}

save(list = c("grid", "gridNum", "gridNew", 
	"ura4.nucleosomes", "seqUra4_0.5k.Dyaddf"), file = "lHBA_ura4_Dyad.RData")


seqUra4_0.5k.Linkerdf <- data.frame(pos = 74:1722, HBA = numeric(length = 1649), stringsAsFactors = FALSE)
for(i in 1:1649){
	seq <- substr(x = seqUra4_0.5k.Linker, start = i, stop = i+146)
	seqUra4_0.5k.Linkerdf$HBA[i] <- HBA(seq, species = "sp", silent = TRUE)
}

seqUra4_0.5k.Linkerdf$lHBA_A <- as.numeric(NA)
seqUra4_0.5k.Linkerdf$lHBA_B <- as.numeric(NA)
seqUra4_0.5k.Linkerdf$lHBA_C <- as.numeric(NA)
seqUra4_0.5k.Linkerdf$lHBA_D <- as.numeric(NA)
seqUra4_0.5k.Linkerdf$lHBA_E <- as.numeric(NA)
seqUra4_0.5k.Linkerdf$lHBA_F <- as.numeric(NA)
seqUra4_0.5k.Linkerdf$lHBA_G <- as.numeric(NA)
seqUra4_0.5k.Linkerdf$lHBA_H <- as.numeric(NA)
seqUra4_0.5k.Linkerdf$lHBA_I <- as.numeric(NA)
seqUra4_0.5k.Linkerdf$lHBA_J <- as.numeric(NA)
seqUra4_0.5k.Linkerdf$lHBA_K <- as.numeric(NA)
seqUra4_0.5k.Linkerdf$lHBA_L <- as.numeric(NA)
seqUra4_0.5k.Linkerdf$lHBA_M <- as.numeric(NA)

for(i in 1:1649){
	seq <- substr(x = seqUra4_0.5k.Linker, start = i, stop = i+146)
	seqUra4_0.5k.Linkerdf[i, 3:15] <- localHBA(seq, species = "sp", silent = TRUE)
}

grid <- as.matrix(seqUra4_0.5k.Linkerdf[, 3:15])
grid <- grid[, ncol(grid):1]
gridNum <- 48
gridNew <- matrix(nrow = nrow(grid), ncol = ncol(grid)*gridNum)
for(i in 1:ncol(grid)) {
	row.range <- ((i-1)*gridNum + 1):((i-1)*gridNum + gridNum)
	gridNew[, row.range] <- grid[,i]
}

save(list = c("grid", "gridNum", "gridNew", 
	"ura4.nucleosomes", "seqUra4_0.5k.Linkerdf"), file = "lHBA_ura4_Linker.RData")


seqUra4_0.5k.Intdf <- data.frame(pos = 74:1722, HBA = numeric(length = 1649), stringsAsFactors = FALSE)
for(i in 1:1649){
	seq <- substr(x = seqUra4_0.5k.Int, start = i, stop = i+146)
	seqUra4_0.5k.Intdf$HBA[i] <- HBA(seq, species = "sp", silent = TRUE)
}

seqUra4_0.5k.Intdf$lHBA_A <- as.numeric(NA)
seqUra4_0.5k.Intdf$lHBA_B <- as.numeric(NA)
seqUra4_0.5k.Intdf$lHBA_C <- as.numeric(NA)
seqUra4_0.5k.Intdf$lHBA_D <- as.numeric(NA)
seqUra4_0.5k.Intdf$lHBA_E <- as.numeric(NA)
seqUra4_0.5k.Intdf$lHBA_F <- as.numeric(NA)
seqUra4_0.5k.Intdf$lHBA_G <- as.numeric(NA)
seqUra4_0.5k.Intdf$lHBA_H <- as.numeric(NA)
seqUra4_0.5k.Intdf$lHBA_I <- as.numeric(NA)
seqUra4_0.5k.Intdf$lHBA_J <- as.numeric(NA)
seqUra4_0.5k.Intdf$lHBA_K <- as.numeric(NA)
seqUra4_0.5k.Intdf$lHBA_L <- as.numeric(NA)
seqUra4_0.5k.Intdf$lHBA_M <- as.numeric(NA)

for(i in 1:1649){
	seq <- substr(x = seqUra4_0.5k.Int, start = i, stop = i+146)
	seqUra4_0.5k.Intdf[i, 3:15] <- localHBA(seq, species = "sp", silent = TRUE)
}

grid <- as.matrix(seqUra4_0.5k.Intdf[, 3:15])
grid <- grid[, ncol(grid):1]
gridNum <- 48
gridNew <- matrix(nrow = nrow(grid), ncol = ncol(grid)*gridNum)
for(i in 1:ncol(grid)) {
	row.range <- ((i-1)*gridNum + 1):((i-1)*gridNum + gridNum)
	gridNew[, row.range] <- grid[,i]
}

save(list = c("grid", "gridNum", "gridNew", 
	"ura4.nucleosomes", "seqUra4_0.5k.Intdf"), file = "lHBA_ura4_Int.RData")




