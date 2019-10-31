# 予測されたすべての dyad を拾う関数。
# 20160213_NuPoP_1.txt や 20160213_NuPoP_2.txt のものを改変。func.name のかわりに nuc.size
getDyadPred <- function(chr.num = 16, prefix = "sc2", nuc.size = 147){
	chr.names <- character(0)
	file.names <- character(0)
	for(i in 1:chr.num){
		chr.name <- paste("chr", i, sep = "")
		chr.names <- c(chr.names, chr.name)
		file.name <- paste(prefix, "_", nuc.size, "bp_", chr.names[i], ".txt", sep = "")
		file.names <- c(file.names, file.name)

	}
	# out <- data.frame(chr = character(0), pred_pos = numeric(0), P_dyad = numeric(0), 
	# 			Occup = numeric(0), Viterbi = numeric(0), Affinity = numeric(0))
	out <- data.frame(chr = character(0), pred_pos = integer(0), P_dyad = numeric(0), 
				Occup = numeric(0), Viterbi = integer(0), Affinity = numeric(0))

	for(i in 1:chr.num){
		# 染色体単位で読み込む。
		cat("chr.names[i]: ", chr.names[i], "\n", sep = "")
		pred.table.chr <- read.table(file.names[i], skip = 1, stringsAsFactors = FALSE)
		colnames(pred.table.chr) <- c("Position", "P_start", "Occup", "Viterbi", "Affinity")

		# Viterbi path の Nuc の最初の塩基を 1 とし、Nuc の残りを 2 とする。
		# for(k in 1:nrow(pred.table.chr)){
		# 	vit <- pred.table.chr$Viterbi[k]
		# 	if(vit == 1) {
		# 		pred.table.chr$Viterbi[(k+1):(k+146)] <- 2
		# 	}	
		# }
		cat("nrow(pred.table.chr) before adding 2:", nrow(pred.table.chr), "\n")
		l <- nrow(pred.table.chr)

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

		cat("nrow(pred.table.chr) after adding 2:", nrow(pred.table.chr), "\n")
		pred.table.chr <- pred.table.chr[1:l,]
		cat("nrow(pred.table.chr) after correction:", nrow(pred.table.chr), "\n")
	
		# pred.table.chr に P_dyad と Viterbi2 カラムを入れる。
		# 5'-end から dyad の位置を求める。Viterbi2 == 1 の dyad は Most probable dyad。
		pred.table.chr$P_dyad <- as.integer(0)
		pred.table.chr$P_dyad[74:nrow(pred.table.chr)] <- pred.table.chr$P_start[1:(nrow(pred.table.chr)-73)]
		pred.table.chr$Viterbi2 <- as.integer(0)
		pred.table.chr$Viterbi2[74:nrow(pred.table.chr)] <- pred.table.chr$Viterbi[1:(nrow(pred.table.chr)-73)]

		# P_dyad が 0 で無いものを選抜する。
		cat("nrow(pred.table.chr) before selection:", nrow(pred.table.chr), "\n")
		pred.table.chr <- subset(pred.table.chr, P_dyad > 0)	# P_dyad が 0 ではないものを選抜
		cat("nrow(pred.table.chr) after selection:", nrow(pred.table.chr), "\n")
		
		lines <- nrow(pred.table.chr)
		out.chr <- data.frame(chr = character(lines), pred_pos = numeric(lines), 
				P_dyad = numeric(lines), 
				Occup = numeric(lines), Viterbi = numeric(lines), Affinity = numeric(lines))
		out.chr$chr <- chr.names[i]
		out.chr$pred_pos <- pred.table.chr$Position
		out.chr$P_dyad <- pred.table.chr$P_dyad
		out.chr$Occup <- pred.table.chr$Occup
		out.chr$Viterbi <- pred.table.chr$Viterbi2
		out.chr$Affinity <- pred.table.chr$Affinity
		out <- rbind(out, out.chr)
	}
	return(out)
}
