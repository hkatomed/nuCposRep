setwd("scripts")
source(file = "dHMM2wig.R")

setwd("../")
setwd("sc_genome_2011")
setwd("predNuPoP_sc")
dHMM2wig(genome.name = "sc_genome_2011", species = "sc", func.name = "predNuPoP", 
		nuc.size = 147, prefix = "NuPoPScSc")
setwd("../predNuPoP_sp")
dHMM2wig(genome.name = "sc_genome_2011", species = "sp", func.name = "predNuPoP", 
		nuc.size = 147, prefix = "NuPoPScSp")
setwd("../predNuCpos_sc")
dHMM2wig(genome.name = "sc_genome_2011", species = "sc", func.name = "predNuCpos", 
		nuc.size = 147, prefix = "nuCposScSc")
setwd("../predNuCpos_sp")
dHMM2wig(genome.name = "sc_genome_2011", species = "sp", func.name = "predNuCpos", 
		nuc.size = 147, prefix = "nuCposScSp")

setwd("../../sp_genome")
setwd("predNuPoP_sc")
dHMM2wig(genome.name = "sp_genome", species = "sc", func.name = "predNuPoP", 
		nuc.size = 147, prefix = "NuPoPSpSc")
setwd("../predNuPoP_sp")
dHMM2wig(genome.name = "sp_genome", species = "sp", func.name = "predNuPoP", 
		nuc.size = 147, prefix = "NuPoPSpSp")
setwd("../predNuCpos_sc")
dHMM2wig(genome.name = "sp_genome", species = "sc", func.name = "predNuCpos", 
		nuc.size = 147, prefix = "nuCposSpSc")
setwd("../predNuCpos_sp")
dHMM2wig(genome.name = "sp_genome", species = "sp", func.name = "predNuCpos", 
		nuc.size = 147, prefix = "nuCposSpSp")
