setwd("genome")
library("Biostrings")
sp.genome <- readDNAStringSet(filepath = "Schizosaccharomyces_pombe.ASM294v2.dna.toplevel.fa.gz")
sp.genome <- sp.genome[1:3]
names(sp.genome) <- c("chr1", "chr2", "chr3")

setwd("../")
dir.create("sp_genome")
setwd("sp_genome")
dir.create("predNuPoP_sc")
dir.create("predNuPoP_sp")
dir.create("predNuCpos_sc")
dir.create("predNuCpos_sp")

for(i in 1:length(sp.genome)){
chrDNA <- sp.genome[i]
filename <- paste("sp_genome_chr", i, ".fasta", sep = "")
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
