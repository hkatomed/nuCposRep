setwd("genome")
library(Biostrings)
sc.genome.2011 <- readDNAStringSet(filepath = "S288C_reference_sequence_R64-1-1_20110203.fsa") # bug fixed
sc.genome.2011 <- sc.genome.2011[1:16]
names(sc.genome.2011) <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
	"chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16")
sc.genome <- sc.genome.2011
sp.genome <- readDNAStringSet(filepath = "Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa.gz")
sp.genome <- sp.genome[1:3]
names(sp.genome) <- c("chr1", "chr2", "chr3")

setwd("../RData")
save(sc.genome.2011, file = "sc_genome_2011.RData")
save(sc.genome, file = "sc_genome.RData")
save(sp.genome, file = "sp_genome.RData")


