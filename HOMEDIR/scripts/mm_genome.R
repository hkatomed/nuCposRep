setwd("genome")
library("Biostrings")
mm9 <- readDNAStringSet(filepath = "Mus_musculus.NCBIM37.67.dna.toplevel.fa.gz")
mm9 <- mm9[191:211]
names(mm9) <- c("chrY", "chr19", "chr18", "chr17", "chr16", "chr15", 
	"chr13", "chr12", "chr11", "chr9", "chr14", "chr10", "chr8", "chr6", "chr7", 
	"chr5", "chr4", "chr3", "chrX", "chr2", "chr1")

setwd("../")
dir.create("mm_genome")
setwd("mm_genome")
dir.create("predNuPoP_mm")
dir.create("predNuCpos_mm")

setwd("predNuPoP_mm")
writeXStringSet(mm9["chr19"], filepath = "mm_genome_chr19.fasta")
setwd("../predNuCpos_mm")
writeXStringSet(mm9["chr19"], filepath = "mm_genome_chr19.fasta")

mm9.unmasked <- mm9
setwd("../../RData")
save(mm9.unmasked, file = "mm9_unmasked.RData")
