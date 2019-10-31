#########
# Viterbi path の 1 (nucleosome state) が始まる位置を得る必要がある。
# 20151224_NuPoP_genome_1.txt の関数を改善する。func.name のかわりに nuc.size
nearestDyadVit <- function(dyad.table, chr.num = 16, prefix = "sc2", nuc.size = 147){
	require("parallel")
	chr.names <- character(0)
	file.names <- character(0)
	for(i in 1:chr.num){
		chr.name <- paste("chr", i, sep = "")
		chr.names <- c(chr.names, chr.name)
		file.name <- paste(prefix, "_", nuc.size, "bp_", chr.names[i], ".txt", sep = "")
		file.names <- c(file.names, file.name)

	}
	vit.nuc.num <- 0
	out <- data.frame(chr = character(0), pred_pos = integer(0), P_dyad = numeric(0), dist = integer(0), 
				chem_pos = integer(0), NCP_score = numeric(0), NCP_SN = numeric(0))
	out2 <- data.frame(chr = character(0), pred_pos = integer(0), P_dyad = numeric(0), dist = integer(0), 
				chem_pos = integer(0), NCP_score = numeric(0), NCP_SN = numeric(0))

	for(i in 1:chr.num){
		cat("chr.names[i]: ", chr.names[i], "\n", sep = "")
		dyad.table.chr <- subset(dyad.table, chr == chr.names[i])
		cat("nrow(dyad.table.chr): ", nrow(dyad.table.chr), "\n")
		pred.table.chr <- read.table(file.names[i], skip = 1, stringsAsFactors = FALSE)
		colnames(pred.table.chr) <- c("Position", "P_start", "Occup", "Viterbi", "Affinity")

		cat("nrow(pred.table.chr) before adding 2:", nrow(pred.table.chr), "\n")
		# l <- nrow(pred.table.chr)
		# for(k in 1:(nrow(pred.table.chr)-146)){
		# 	if(pred.table.chr$Viterbi[k] == 1) {
		# 		pred.table.chr$Viterbi[(k+1):(k+146)] <- 2
		# 	}	
		# 	k <- k + 146
		# }

		Viterbi.new <- paste(pred.table.chr$Viterbi, sep = "", collapse = "")
		if(nuc.size == 147){
			Vit.1 <- paste("11111111111111111111111111111111111111111111111111", 
					"11111111111111111111111111111111111111111111111111", 
					"11111111111111111111111111111111111111111111111", sep = "")
			Vit.2 <- paste("12222222222222222222222222222222222222222222222222",
					"22222222222222222222222222222222222222222222222222",
					"22222222222222222222222222222222222222222222222", sep = "")
		}else if(nuc.size == 137){
			Vit.1 <- paste("11111111111111111111111111111111111111111111111111", 
					"11111111111111111111111111111111111111111111111111", 
					"1111111111111111111111111111111111111", sep = "")
			Vit.2 <- paste("12222222222222222222222222222222222222222222222222",
					"22222222222222222222222222222222222222222222222222",
					"2222222222222222222222222222222222222", sep = "")
		}else if(nuc.size == 127){
			Vit.1 <- paste("11111111111111111111111111111111111111111111111111", 
					"11111111111111111111111111111111111111111111111111", 
					"111111111111111111111111111", sep = "")
			Vit.2 <- paste("12222222222222222222222222222222222222222222222222",
					"22222222222222222222222222222222222222222222222222",
					"222222222222222222222222222", sep = "")
		}else if(nuc.size == 117){
			Vit.1 <- paste("11111111111111111111111111111111111111111111111111", 
					"11111111111111111111111111111111111111111111111111", 
					"11111111111111111", sep = "")
			Vit.2 <- paste("12222222222222222222222222222222222222222222222222",
					"22222222222222222222222222222222222222222222222222",
					"22222222222222222", sep = "")
		}else if(nuc.size == 107){
			Vit.1 <- paste("11111111111111111111111111111111111111111111111111", 
					"11111111111111111111111111111111111111111111111111", 
					"1111111", sep = "")
			Vit.2 <- paste("12222222222222222222222222222222222222222222222222",
					"22222222222222222222222222222222222222222222222222",
					"2222222", sep = "")
		}
		
		pred.table.chr$Viterbi <- as.integer(strsplit(gsub(Vit.1, Vit.2, Viterbi.new), split = "")[[1]])

		pred.table.chr <- subset(pred.table.chr, Viterbi == 1)
		cat("nrow(pred.table.chr):", nrow(pred.table.chr), "\n")

		# ChemVsPred
		cat("ChemVsPred...\n")
		lines <- nrow(dyad.table.chr)
		out.chr <- data.frame(chr = character(lines), pred_pos = integer(lines), 
				P_dyad = numeric(lines), dist = integer(lines), 
				chem_pos = integer(lines), NCP_score = numeric(lines), NCP_SN = numeric(lines))
		out.chr$chr <- chr.names[i]
		pred.pos <- as.integer(pred.table.chr$Position + 73)	# 5'-end から dyad の位置を求める
		pred.P_dyad <- pred.table.chr$P_start
		cat("head(pred.pos, n = 20): ", head(pred.pos, n = 20), "\n")
		cat("head(pred.P_dyad, n = 20): ", head(pred.P_dyad, n = 20), "\n")
		cat("length(pred.pos): ", length(pred.pos), "\n")
		cat("pred.pos[1]: ", pred.pos[1], "\n")
		cat("pred.pos[length(pred.pos)]: ", pred.pos[length(pred.pos)], "\n")

			get.nearest.indexes <- function(pos){
				return(which(abs(pred.pos - pos) == min(abs(pred.pos - pos)))[1])
			}
		nearest.indexs <- integer(lines)
		nearest.indexs <- unlist(mcMap(f = get.nearest.indexes, pos = as.integer(dyad.table.chr$pos),
						mc.cores = 20))

		out.chr$pred_pos <- pred.pos[nearest.indexs]
		out.chr$P_dyad <- pred.P_dyad[nearest.indexs]
		out.chr$chem_pos <- dyad.table.chr$pos
		out.chr$NCP_score <- dyad.table.chr$NCP_score
		out.chr$NCP_SN <- dyad.table.chr$NCP_SN
		out.chr$dist <- out.chr$pred_pos - as.integer(dyad.table.chr$pos)

		out <- rbind(out, out.chr)
		vit.nuc.num <- vit.nuc.num + nrow(pred.table.chr)

		# PredVsChem
		cat("PredVsChem...\n")
		lines <- nrow(pred.table.chr)
		out2.chr <- data.frame(chr = character(lines), pred_pos = numeric(lines), 
				P_dyad = numeric(lines), dist = numeric(lines), 
				chem_pos = numeric(lines), NCP_score = numeric(lines), NCP_SN = numeric(lines))
		out2.chr$chr <- chr.names[i]
		chem.pos <- as.numeric(dyad.table.chr$pos)
		NCP_score <- as.numeric(dyad.table.chr$NCP_score)
		NCP_SN <- as.numeric(dyad.table.chr$NCP_SN)


			get.nearest.indexes <- function(pos){
				return(which(abs(chem.pos - pos) == min(abs(chem.pos - pos)))[1])
			}
		nearest.indexs <- integer(lines)
		# nearest.indexs <- unlist(mcMap(f = get.nearest.indexes, 
		# 				pos = as.integer(pred.table.chr$Position),
		# 				mc.cores = 20))
		nearest.indexs <- unlist(mcMap(f = get.nearest.indexes, 
						pos = as.integer(pred.table.chr$Position + 73),
						mc.cores = 20))

		out2.chr$pred_pos <- as.integer(pred.table.chr$Position + 73)
		out2.chr$P_dyad <- pred.table.chr$P_start
		out2.chr$chem_pos <- as.integer(chem.pos[nearest.indexs])
		out2.chr$NCP_score <- NCP_score[nearest.indexs]
		out2.chr$NCP_SN <-  NCP_SN[nearest.indexs]
		# out2.chr$dist <- as.integer(out2.chr$chem_pos - pred.table.chr$Position)
		# out2.chr$dist <- as.integer(out2.chr$chem_pos - pred.table.chr$Position + 73)
		out2.chr$dist <- as.integer(out2.chr$chem_pos - (pred.table.chr$Position + 73))

		out2 <- rbind(out2, out2.chr)
	}

	cat("vit.nuc.num: ", vit.nuc.num, "\n", sep = "")

	return(list(ChemVsPred = out, PredVsChem = out2))
}
