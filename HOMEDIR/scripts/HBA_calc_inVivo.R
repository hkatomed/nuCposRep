setwd("src")
dyn.load("HBA_3.so")
setwd("../scripts")
source(file = "HBAm.R")
source(file = "HBAc.R")
setwd("../RData")
load(file = "sysdata_NuPoP.rda")
load(file = "sysdata_nuCpos.rda")

setwd("../Chereji")
library(openxlsx)
Chereji <- read.xlsx(xlsxFile = "13059_2018_1398_MOESM2_ESM.xlsx", colNames = TRUE)

names(Chereji)[8] <- "plus1"
names(Chereji)[9] <- "minus1"

chr.lengths <- c(230218, 813184, 316620, 1531933, 439888, 576874, 
			270161, 1090940, 562643, 745751, 666816, 1078177, 
			924431, 784333, 1091291, 948066)
chr.names <- c("chr1", "chr2", "chr3", "chr4", "chr9", "chr5", 
		"chr6", "chr7", "chr8", "chr10", "chr11", "chr12", 
		"chr13", "chr14", "chr15", "chr16") 
chr.names2 <- c("chrI", "chrII", "chrIII", "chrIV", "chrIX", "chrV", 
		"chrVI", "chrVII", "chrVIII", "chrX", "chrXI", "chrXII", 
		"chrXIII", "chrXIV", "chrXV", "chrXVI") 

for(i in 1:16)	Chereji$Chr[Chereji$Chr == chr.names2[i]] <- chr.names[i]

library(Biostrings)
setwd("../RData")
load(file = "sc_genome_2011.RData")


Chereji$plus1.seq <- ""
Chereji$minus1.seq <- ""

for(i in 1:nrow(Chereji)){
	if(Chereji$Strand[i] == 1){
	Chereji$plus1.seq[i] <- as.character(
		sc.genome.2011[[Chereji$Chr[i]]][(Chereji$plus1[i]-173):(Chereji$plus1[i]+173)])
	Chereji$minus1.seq[i] <- as.character(
		sc.genome.2011[[Chereji$Chr[i]]][(Chereji$minus1[i]-173):(Chereji$minus1[i]+173)])
	}
	if(Chereji$Strand[i] == -1){
	Chereji$plus1.seq[i] <- as.character(reverseComplement(
		sc.genome.2011[[Chereji$Chr[i]]][(Chereji$plus1[i]-173):(Chereji$plus1[i]+173)]))
	Chereji$minus1.seq[i] <- as.character(reverseComplement(
		sc.genome.2011[[Chereji$Chr[i]]][(Chereji$minus1[i]-173):(Chereji$minus1[i]+173)]))
	}
}

Chereji$plus1.347 <- FALSE
Chereji$minus1.347 <- FALSE

for(i in 1:nrow(Chereji)){
	if(nchar(Chereji$plus1.seq[i]) == 347) Chereji$plus1.347[i] <- TRUE
	if(nchar(Chereji$minus1.seq[i]) == 347) Chereji$minus1.347[i] <- TRUE
}

Chereji.plus1.HBAm <- matrix(nrow = nrow(Chereji), ncol = 201)
rownames(Chereji.plus1.HBAm) <- Chereji$ORF
colnames(Chereji.plus1.HBAm) <- -100:100

Chereji.minus1.HBAm <- matrix(nrow = nrow(Chereji), ncol = 201)
rownames(Chereji.minus1.HBAm) <- Chereji$ORF
colnames(Chereji.minus1.HBAm) <- -100:100

Chereji.plus1.HBAc <- matrix(nrow = nrow(Chereji), ncol = 201)
rownames(Chereji.plus1.HBAc) <- Chereji$ORF
colnames(Chereji.plus1.HBAc) <- -100:100

Chereji.minus1.HBAc <- matrix(nrow = nrow(Chereji), ncol = 201)
rownames(Chereji.minus1.HBAc) <- Chereji$ORF
colnames(Chereji.minus1.HBAc) <- -100:100

cat("nrow(Chereji):", nrow(Chereji), "\n", sep = "")
for(i in 1:nrow(Chereji)){
	if(i%%100 == 0) cat(i, ",", sep = "")
	plus1.seq <- Chereji$plus1.seq[i]
	minus1.seq <- Chereji$minus1.seq[i]
	for(j in 1:201){
		p1seq <- substring(plus1.seq, first = j, last = j + 146)
		Chereji.plus1.HBAm[i, j] <- HBAm(p1seq, silent = TRUE)
		Chereji.plus1.HBAc[i, j] <- HBAc(p1seq, silent = TRUE)
		m1seq <- substring(minus1.seq, first = j, last = j + 146)
		Chereji.minus1.HBAm[i, j] <- HBAm(m1seq, silent = TRUE)
		Chereji.minus1.HBAc[i, j] <- HBAc(m1seq, silent = TRUE)
	}
}

save(Chereji.plus1.HBAm, file = "Chereji_plus1_HBAm.RData")
save(Chereji.plus1.HBAc, file = "Chereji_plus1_HBAc.RData")
save(Chereji.minus1.HBAm, file = "Chereji_minus1_HBAm.RData")
save(Chereji.minus1.HBAc, file = "Chereji_minus1_HBAc.RData")

Chereji.plus1.AT <- matrix(nrow = nrow(Chereji), ncol = 201)
rownames(Chereji.plus1.AT) <- Chereji$ORF
colnames(Chereji.plus1.AT) <- -100:100

Chereji.minus1.AT <- matrix(nrow = nrow(Chereji), ncol = 201)
rownames(Chereji.minus1.AT) <- Chereji$ORF
colnames(Chereji.minus1.AT) <- -100:100

cat("nrow(Chereji):", nrow(Chereji), "\n", sep = "")
for(i in 1:nrow(Chereji)){
	if(i%%100 == 0) cat(i, ",", sep = "")
	plus1.seq <- Chereji$plus1.seq[i]
	minus1.seq <- Chereji$minus1.seq[i]
	for(j in 1:201){
		p1seq <- substring(plus1.seq, first = j, last = j + 146)
		Chereji.plus1.AT[i, j] <- letterFrequency(DNAString(p1seq), 
			letters = c("AT"), OR = "|", as.prob = TRUE)
		m1seq <- substring(minus1.seq, first = j, last = j + 146)
		Chereji.minus1.AT[i, j] <- letterFrequency(DNAString(m1seq), 
			letters = c("AT"), OR = "|", as.prob = TRUE)
	}
}

save(Chereji.plus1.AT, file = "Chereji_plus1_AT.RData")
save(Chereji.minus1.AT, file = "Chereji_minus1_AT.RData")

