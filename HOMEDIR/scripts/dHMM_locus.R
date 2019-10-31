
library(Biostrings)
library(nuCpos)
library(NuPoP)

wd <- getwd()
setwd("src")
dyn.load("HBA_3.so")
dyn.load("HBA_3sp.so")
setwd("../scripts")
source(file = "HBAm.R")
source(file = "HBAc.R")
setwd("../RData")
load(file = "sysdata_NuPoP.rda")
load(file = "sysdata_nuCpos.rda")
setwd("../seq")
TRP1ARS1 <- readDNAStringSet(filepath = "TRP1ARS1x1.fasta")

TRP1ARS1 <- paste(as.character(TRP1ARS1), as.character(TRP1ARS1), 
	as.character(TRP1ARS1), as.character(TRP1ARS1), as.character(TRP1ARS1), 
	collapse = "", sep = "")

TRP1ARS1 <- DNAStringSet(TRP1ARS1)
names(TRP1ARS1) <- "TRP1ARS1"

setwd("../")
dir.create("TRP1ARS1")
setwd("TRP1ARS1")
dir.create("nuCpos")
dir.create("NuPoP")
writeXStringSet(TRP1ARS1, filepath = "nuCpos/TRP1ARS1.fasta")
writeXStringSet(TRP1ARS1, filepath = "NuPoP/TRP1ARS1.fasta")

setwd("nuCpos")
predNuCpos(file = "TRP1ARS1.fasta", species = "sc", smoothHBA = FALSE, ActLikePredNuPoP = TRUE)
results.nuCpos.Sc <- readNuPoP(file = "TRP1ARS1.fasta_Prediction4.txt", startPos = (2931-1000), 
	endPos = (4395+1000))
setwd("../NuPoP")
predNuPoP(file = "TRP1ARS1.fasta", species = 7, model = 4)
results.NuPoP.Sc <- readNuPoP(file = "TRP1ARS1.fasta_Prediction4.txt", startPos = (2931-1000), 
	endPos = (4395+1000))

TRP1ARS1 <- as.character(TRP1ARS1)
TRP1ARS1HBA.m <- data.frame(pos = 74:(nchar(TRP1ARS1)-73), 
	HBA = numeric(length = length(seq(74, nchar(TRP1ARS1)-73))), stringsAsFactors = FALSE)
for(i in 1:length(seq(74, nchar(TRP1ARS1)-73))){
	seq <- substr(x = TRP1ARS1, start = i, stop = i+146)
	TRP1ARS1HBA.m$HBA[i] <- HBAm(seq, silent = TRUE)
}

results.NuPoP.Sc$Affinity <- TRP1ARS1HBA.m[1858:5322,2]

TRP1ARS1.AT <- numeric(length = nchar(TRP1ARS1))
TRP1ARS1.AT[1:nchar(TRP1ARS1)] <- NA 

for(i in 74:(nchar(TRP1ARS1)-73)){
	testseq <- substring(TRP1ARS1, first = i-73, last = i+73)
	TRP1ARS1.AT[i] <- letterFrequency(DNAString(testseq), 
			letters = c("AT"), OR = "|", as.prob = TRUE)
}

setwd(wd)
setwd("RData")
save(results.nuCpos.Sc, file = "results_nuCpos_Sc_TRP1ARS1.RData")
save(results.NuPoP.Sc, file = "results_NuPoP_Sc_TRP1ARS1.RData")
save(TRP1ARS1.AT, file = "TRP1ARS1_AT.RData")



setwd(wd)
setwd("RData")
load(file = "sp_genome.RData")

seqUra4 <- reverseComplement(sp.genome[["chr3"]][(115781-3000):(116575+3000)])
seqUra4 <- DNAStringSet(seqUra4)
names(seqUra4) <- "ura4_3kb"

setwd("../")
dir.create("ura4")
setwd("ura4")
dir.create("nuCpos")
dir.create("NuPoP")
writeXStringSet(seqUra4, filepath = "nuCpos/ura4_3kb.fasta")
writeXStringSet(seqUra4, filepath = "NuPoP/ura4_3kb.fasta")


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

setwd("nuCpos")
predNuCpos(file = "ura4_3kb.fasta", species = "sp", smoothHBA = FALSE, ActLikePredNuPoP = TRUE)
results.nuCpos <- readNuPoP(file = "ura4_3kb.fasta_Prediction4.txt", startPos = 2000, endPos = 4795)
setwd("../NuPoP")
predNuPoP(file = "ura4_3kb.fasta", species = 9, model = 4)
results.NuPoP <- readNuPoP(file = "ura4_3kb.fasta_Prediction4.txt", startPos = 2000, endPos = 4795)

seqUra4 <- as.character(seqUra4)
seqUra4HBA.m <- data.frame(pos = 74:(nchar(seqUra4)-73), 
	HBA = numeric(length = length(seq(74, nchar(seqUra4)-73))), stringsAsFactors = FALSE)
for(i in 1:length(seq(74, nchar(seqUra4)-73))){
	seq <- substr(x = seqUra4, start = i, stop = i+146)
	seqUra4HBA.m$HBA[i] <- HBAm(seq, silent = TRUE, species = "sp")
}

results.NuPoP$Affinity <- seqUra4HBA.m[1927:4722,2]

seqUra4.AT <- numeric(length = nchar(seqUra4))
seqUra4.AT[1:nchar(seqUra4)] <- NA 

for(i in 74:(nchar(seqUra4)-73)){
	testseq <- substring(seqUra4, first = i-73, last = i+73)
	seqUra4.AT[i] <- letterFrequency(DNAString(testseq), 
			letters = c("AT"), OR = "|", as.prob = TRUE)
}

setwd(wd)
setwd("RData")
save(results.NuPoP, file = "results_NuPoP_ura4.RData")
save(results.nuCpos, file = "results_nuCpos_ura4.RData")
save(dyads, file = "ura4_dyads.RData")
save(seqUra4.AT, file = "seqUra4_AT.RData")

