# 予測された dyad のうち、chemical dyad に最も近いものを拾う関数。
# 20160213_NuPoP_1.txt や 20160213_NuPoP_2.txt のものを改変。func.name のかわりに nuc.size
nearestDyadPredMm <- function(dyad.table, prefix = "NuPoPMmMm", nuc.size = 147, mm.chr = "chr19"){
	require("parallel")
	chr.name <- mm.chr
	file.name <- paste(prefix, "_", nuc.size, "bp_", chr.name, ".txt", sep = "")
	
	out <- data.frame(chr = character(0), pred_pos = integer(0), P_dyad = numeric(0), dist = integer(0), 
				chem_pos = integer(0), NCP_score = numeric(0), NCP_SN = numeric(0), 
				Occup = numeric(0), Viterbi = integer(0), Affinity = numeric(0))

		cat("chr.name: ", chr.name, "\n", sep = "")
		dyad.table.chr <- dyad.table
		cat("nrow(dyad.table.chr): ", nrow(dyad.table.chr), "\n")
		pred.table.chr <- read.table(file.name, skip = 1, stringsAsFactors = FALSE)
		colnames(pred.table.chr) <- c("Position", "P_start", "Occup", "Viterbi", "Affinity")

		# Viterbi path の Nuc の最初の塩基を 1 とし、Nuc の残りを 2 とする。
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

		# cat("nrow(pred.table.chr) after adding 2:", nrow(pred.table.chr), "\n")
		# pred.table.chr <- pred.table.chr[1:l,]
		# cat("nrow(pred.table.chr) after correction:", nrow(pred.table.chr), "\n")
	
		# pred.table.chr に P_dyad と Viterbi2 カラムを入れる。
		# 5'-end から dyad の位置を求める。Viterbi2 == 1 の dyad は Most probable dyad。
		pred.table.chr$P_dyad <- as.integer(0)
		pred.table.chr$P_dyad[74:nrow(pred.table.chr)] <- pred.table.chr$P_start[1:(nrow(pred.table.chr)-73)]
		pred.table.chr$Viterbi2 <- as.integer(0)
		pred.table.chr$Viterbi2[74:nrow(pred.table.chr)] <- pred.table.chr$Viterbi[1:(nrow(pred.table.chr)-73)]

		# 
		cat("nrow(pred.table.chr) before selection:", nrow(pred.table.chr), "\n")

		# ChemVsPredAll
		cat("ChemVsPred...\n")
		lines <- nrow(dyad.table.chr)
		out.chr <- data.frame(chr = character(lines), pred_pos = integer(lines), 
				P_dyad = numeric(lines), dist = integer(lines), 
				chem_pos = integer(lines), NCP_score = numeric(lines), NCP_SN = numeric(lines), 
				Occup = numeric(lines), Viterbi = integer(lines), Affinity = numeric(lines))
		out.chr$chr <- chr.name
		pred.table.chr <- subset(pred.table.chr, P_dyad > 0)	# P_start が 0 ではないものを選抜
		cat("nrow(pred.table.chr) after selection:", nrow(pred.table.chr), "\n")
		pred.pos <- as.integer(pred.table.chr$Position)
		pred.P_dyad <- pred.table.chr$P_dyad
		pred.Occup <- pred.table.chr$Occup
		pred.Viterbi <- pred.table.chr$Viterbi2
		pred.Affinity <- pred.table.chr$Affinity
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
						mc.cores = 4))

		out.chr$pred_pos <- pred.pos[nearest.indexs]
		out.chr$P_dyad <- pred.P_dyad[nearest.indexs]
		out.chr$chem_pos <- dyad.table.chr$pos
		out.chr$NCP_score <- dyad.table.chr$NCP_score
		out.chr$NCP_SN <- dyad.table.chr$NCP_SN
		out.chr$dist <- out.chr$pred_pos - as.integer(dyad.table.chr$pos)
		out.chr$Occup <- pred.Occup[nearest.indexs]
		out.chr$Viterbi <- pred.Viterbi[nearest.indexs]
		out.chr$Affinity <- pred.Affinity[nearest.indexs]

		out <- rbind(out, out.chr)
	
	return(list(ChemVsPredAll = out))
}
