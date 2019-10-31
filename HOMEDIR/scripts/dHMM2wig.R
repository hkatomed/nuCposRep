## wig を保存する関数（20190109 変更なし）
save.wig <- function(pos, value, file.name, chr.name = "chr1", func.name = "predNuPoP", 
				description = "Predicted nucleosome occupancy"){
		# wig のヘッダーをつくる関数
		get.wig.head <- function(chr.name = "chr1", func.name = "predNuPoP", 
						description = "Predicted nucleosome occupancy"){
			wig.1st.line <- paste("track type= wiggle_0 name= \"", 
				chr.name, "_", func.name, 
				"\" description= \"", description, "\"", 
				sep = "")
			wig.2nd.line <- paste("variableStep chrom=", chr.name, " span=1", sep="")
			wig.head <- c(wig.1st.line, wig.2nd.line)
			return(wig.head)
		} # get.wig.head
	wig.head <- get.wig.head(chr.name = chr.name, func.name = func.name, description = description)
	wig.body <- paste(pos, value, sep = "\t")
	wig.total <- c(wig.head, wig.body)
	write(wig.total, file = file.name)
}

# 各カラムから wig をつくる関数。save.wig() を使う。（20190109 変更なし）
save.wigs <- function(NuPoP.out, chr.name = "chr1", func.name = "predNuPoP", nuc.size = 147, prefix){

	descriptions <- c("P_start", "Occupancy", "Viterbi", "Affinity")
	pos <- NuPoP.out$Position
	for(i in 1:4){
		# if(i == 4) NuPoP.out <- subset(NuPoP.out, is.na(Affinity) == FALSE)
		NuPoP.out$Affinity[is.na(NuPoP.out$Affinity) == TRUE] <- 0
		NuPoP.out$Affinity[is.nan(NuPoP.out$Affinity) == TRUE] <- 0
		file.name <- paste(prefix, "_", nuc.size, "bp_", chr.name, "_", descriptions[i], ".wig", sep = "")
		save.wig(pos = pos, value = NuPoP.out[,i+1], file.name = file.name, 
				chr.name = chr.name, func.name = func.name, 
				description = descriptions[i])
		if(i == 1){
			file.name <- paste(prefix, "_", nuc.size, "bp_", chr.name, "_", "P_dyad", ".wig", sep = "")
			value <- NuPoP.out[,2]
			cat("length(value):", length(value), ";")
			value2 <- value * 0
			cat("length(value):", length(value), "\n")
			value2[(nuc.size %/% 2 + 1):(length(value))] <- 
					value[1:(length(value) - nuc.size %/% 2)]
			save.wig(pos = pos, value = value2, file.name = file.name, 
					chr.name = chr.name, func.name = func.name, 
					description = "P_dyad")
		}
	}
	# P_start から P_dyad をつくる
}

# まとめて全ゲノムを解析する関数（20190109: use.original, apply.147.linker, do.NuPoP を廃止した）
dHMM2wig <- function(genome.name = "sc_genome_2011", species = "sc", func.name = "predNuPoP", 
		nuc.size = 147, prefix = "ori", mm.chr = "chr19"){
	if(species == "sc") species.num <- 7
	if(species == "sp") species.num <- 9
	if(species == "mm") species.num <- 2
	if(genome.name == "sc_genome_2011") chr.num = 16
	if(genome.name == "sp_genome") chr.num = 3
	if(genome.name == "mm_genome") chr.num = 1
	chr.names <- character(0)
	file.names <- character(0)
	file.names2 <- character(0)
	file.names3 <- character(0)
	for(i in 1:chr.num){
		chr.name <- paste("chr", i, sep = "")
		chr.names <- c(chr.names, chr.name)
		if(genome.name == "mm_genome") chr.names <- mm.chr
		file.name <- paste(genome.name, "_", chr.names[i], ".fasta", sep = "")
		file.names <- c(file.names, file.name)
		file.name2 <- paste(file.name, "_Prediction4.txt", sep = "")
		file.names2 <- c(file.names2, file.name2)
		file.name3 <- paste(prefix, "_", nuc.size, "bp_", chr.names[i], ".txt", sep = "")
		file.names3 <- c(file.names3, file.name3)
	}
	for(i in 1:chr.num){
		if(func.name == "predNuPoP"){
			require("NuPoP")
			predNuPoP(file = file.names[i], species = species.num)
		}
		if(func.name == "predNuCpos"){
			require("nuCpos")
			predNuCpos(file = file.names[i], inseq = "", species = species, 
					ActLikePredNuPoP = TRUE)
		}

		file.rename(from = file.names2[i], to = file.names3[i])
		
		NuPoP.out <- read.table(file = file.names3[i], stringsAsFactors = FALSE, header = TRUE)
		save.wigs(NuPoP.out, chr.name = chr.names[i], func.name = func.name, 
			nuc.size = nuc.size, prefix = prefix)
		gc(reset = TRUE)	# ガベージコレクションを実行
	}
}
