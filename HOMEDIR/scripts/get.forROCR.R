get.forROCR <- function(allpred){
	out.allpred <- allpred
	allpred <- allpred[[1]]
	allpred$dup <- as.integer(0)
	allpred$all <- seq(1,nrow(allpred))
	
	cat("nrow(allpred): ", nrow(allpred), "\n")
	# pred_pos で同じ値が続くときに、距離が近い方を採用する。
	dup.count <- as.integer(0)
	for(i in 2:nrow(allpred)){
		if(i %% 1000 == 0) cat(i, ",", sep = "")
	
		if(i == 2){
			if(allpred$pred_pos[i-1] == allpred$pred_pos[i]){
			dup.count <- dup.count + 1
			allpred$dup[i-1:i] <- dup.count
			}
		}
		if(i > 2){
			if(allpred$pred_pos[i-1] == allpred$pred_pos[i]){
				if(allpred$pred_pos[i-2] == allpred$pred_pos[i-1]){
					allpred$dup[(i-1):i] <- dup.count
				}
				if(allpred$pred_pos[i-2] != allpred$pred_pos[i-1]){
					dup.count <- dup.count + 1
					allpred$dup[(i-1):i] <- dup.count
				}
			}
		}
	}
	cat("\n\n")
	
	single <- subset(allpred, dup == 0)
	dup <- subset(allpred, dup != 0)
	
	dup.new <- dup[1:max(allpred$dup),]
	for(i in 1:max(allpred$dup)){
		selected <- subset(dup, dup == i)
		dup.new[i,] <- selected[order(selected$dist)[1],]
	}
	
	allpred.forROCR <- out.allpred
	allpred.forROCR[[1]] <- rbind(single, dup.new)
	allpred.forROCR[[1]] <- allpred.forROCR[[1]][order(allpred.forROCR[[1]]$all),]
	return(allpred.forROCR)
}
