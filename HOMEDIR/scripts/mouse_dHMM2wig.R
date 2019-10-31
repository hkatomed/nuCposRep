setwd("scripts")
source(file = "dHMM2wig.R")

setwd("../")
setwd("mm_genome")
setwd("predNuPoP_mm")
dHMM2wig(genome.name = "mm_genome", species = "mm", func.name = "predNuPoP", 
		nuc.size = 147, prefix = "NuPoPMmMm", mm.chr = "chr19")
setwd("../predNuCpos_mm")
dHMM2wig(genome.name = "mm_genome", species = "mm", func.name = "predNuCpos", 
		nuc.size = 147, prefix = "nuCposMmMm", mm.chr = "chr19")

