setwd("genome")
library("Biostrings")
sc.genome <- readDNAStringSet(filepath = "Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz")
sc.genome <- sc.genome[1:16]
names(sc.genome) <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", 
	"chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16")

setwd("../")
dir.create("sc_genome_2011")
setwd("sc_genome_2011")
dir.create("predNuPoP_sc")
dir.create("predNuPoP_sp")
dir.create("predNuCpos_sc")
dir.create("predNuCpos_sp")

for(i in 1:length(sc.genome)){
chrDNA <- sc.genome[i]
filename <- paste("sc_genome_2011_chr", i, ".fasta", sep = "")
setwd("predNuPoP_sc")
writeXStringSet(chrDNA, filepath = filename)
setwd("../predNuPoP_sp")
writeXStringSet(chrDNA, filepath = filename)
setwd("../predNuCpos_sc")
writeXStringSet(chrDNA, filepath = filename)
setwd("../predNuCpos_sp")
writeXStringSet(chrDNA, filepath = filename)
setwd("../")
}
