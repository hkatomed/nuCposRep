getPredTableMm <- function(prefix = "NuPoPMmMm", nuc.size = 147, mm.chr = "chr19"){

	chr.name <- mm.chr
	file.name <- paste(prefix, "_", nuc.size, "bp_", chr.name, ".txt", sep = "")

	out <- data.frame(chr = character(0), pred_pos = integer(0), P_dyad = numeric(0), 
				Occup = numeric(0), Viterbi = integer(0), Affinity = numeric(0))

	nrowall <- 0

		# 染色体単位で読み込む。
		cat("chr.name: ", chr.name, "\n", sep = "")
		pred.table.chr <- read.table(file.name, skip = 1, stringsAsFactors = FALSE)
		colnames(pred.table.chr) <- c("Position", "P_start", "Occup", "Viterbi", "Affinity")

		cat("nrow(pred.table.chr):", nrow(pred.table.chr), "\n")
		nrowall <- nrowall + nrow(pred.table.chr)

		Viterbi.new <- paste(pred.table.chr$Viterbi, sep = "", collapse = "")
		if(nuc.size == 147){
			Vit.1 <- paste("11111111111111111111111111111111111111111111111111", 
					"11111111111111111111111111111111111111111111111111", 
					"11111111111111111111111111111111111111111111111", sep = "")
			Vit.2 <- paste("10000000000000000000000000000000000000000000000000",
					"00000000000000000000000000000000000000000000000000",
					"00000000000000000000000000000000000000000000000", sep = "")
		}else if(nuc.size == 137){
			Vit.1 <- paste("11111111111111111111111111111111111111111111111111", 
					"11111111111111111111111111111111111111111111111111", 
					"1111111111111111111111111111111111111", sep = "")
			Vit.2 <- paste("10000000000000000000000000000000000000000000000000",
					"00000000000000000000000000000000000000000000000000",
					"0000000000000000000000000000000000000", sep = "")
		}else if(nuc.size == 127){
			Vit.1 <- paste("11111111111111111111111111111111111111111111111111", 
					"11111111111111111111111111111111111111111111111111", 
					"111111111111111111111111111", sep = "")
			Vit.2 <- paste("10000000000000000000000000000000000000000000000000",
					"00000000000000000000000000000000000000000000000000",
					"000000000000000000000000000", sep = "")
		}else if(nuc.size == 117){
			Vit.1 <- paste("11111111111111111111111111111111111111111111111111", 
					"11111111111111111111111111111111111111111111111111", 
					"11111111111111111", sep = "")
			Vit.2 <- paste("10000000000000000000000000000000000000000000000000",
					"00000000000000000000000000000000000000000000000000",
					"00000000000000000", sep = "")
		}else if(nuc.size == 107){
			Vit.1 <- paste("11111111111111111111111111111111111111111111111111", 
					"11111111111111111111111111111111111111111111111111", 
					"1111111", sep = "")
			Vit.2 <- paste("10000000000000000000000000000000000000000000000000",
					"00000000000000000000000000000000000000000000000000",
					"0000000", sep = "")
		}
					
		Viterbi.new <- as.integer(strsplit(gsub(Vit.1, Vit.2, Viterbi.new), split = "")[[1]])

		# pred.table.chr に P_dyad と Viterbi2 カラムを入れる。
		# 5'-end から dyad の位置を求める。Viterbi2 == 1 の dyad は Most probable dyad。
		pred.table.chr$P_dyad <- as.integer(0)
		pred.table.chr$P_dyad[74:nrow(pred.table.chr)] <- pred.table.chr$P_start[1:(nrow(pred.table.chr)-73)]
		pred.table.chr$Viterbi2 <- as.integer(0)
		pred.table.chr$Viterbi2[74:nrow(pred.table.chr)] <- Viterbi.new[1:(nrow(pred.table.chr)-73)]
		
		lines <- nrow(pred.table.chr)
		out.chr <- data.frame(chr = character(lines), pred_pos = numeric(lines), 
				P_dyad = numeric(lines), 
				Occup = numeric(lines), Viterbi = numeric(lines), Affinity = numeric(lines))
		out.chr$chr <- chr.name
		out.chr$pred_pos <- pred.table.chr$Position
		out.chr$P_dyad <- pred.table.chr$P_dyad
		out.chr$Occup <- pred.table.chr$Occup
		out.chr$Viterbi <- pred.table.chr$Viterbi2
		out.chr$Affinity <- pred.table.chr$Affinity
		out <- rbind(out, out.chr)

	cat("nrow(out):", nrow(out), "\n")
	cat("nrowall:", nrowall, "\n")
	obj.name <- paste(prefix, ".", mm.chr, ".allpred", sep = "")
	assign(obj.name, out, envir = .GlobalEnv)
	# return(out)
}
