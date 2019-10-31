pred.summaryMm <- function(prefix = "NuPoPMmMm", mm.chr = "chr19"){
	Target <- "Mm"
	Model <- "Mm"
	
	obj.name <- paste("mm9_", mm.chr, "_Uni", sep = "")
	dyad.table.unique <- get(obj.name)
	dyad.table.unique <- data.frame(pos = dyad.table.unique)
	# obj.name <- paste("mm9_", mm.chr, "_Red", sep = "")
	obj.name <- paste("mm9_", mm.chr, "_Red30", sep = "")
	dyad.table.redundant <- get(obj.name)
	dyad.table.redundant <- data.frame(pos = dyad.table.redundant)
	
	if(length(grep("NuPoP", prefix)) == 1) Package <- "NuPoP"
	if(length(grep("nuCpos", prefix)) == 1) Package <- "nuCpos"
	allpred <- get(paste(prefix, ".", mm.chr, ".allpred", sep = ""), envir = .GlobalEnv)
	Genome.size <- nrow(allpred)
	Num.Unique <- nrow(dyad.table.unique)
	Num.Redundant <- nrow(dyad.table.redundant)
	Predicted <- nrow(subset(allpred, P_dyad > 0))
	Viterbi <- nrow(subset(allpred, Viterbi == 1))
	
	UniqueW0 <- nrow(subset(allpred, Viterbi == 1 & UniqueW0 == 1))
	RedundantW0 <- nrow(subset(allpred, Viterbi == 1 & RedundantW0 == 1))
	Rate.UniqueW0 <- UniqueW0/Viterbi
	Rate.RedundantW0 <- RedundantW0/Viterbi
	
	UniqueW1 <- nrow(subset(allpred, Viterbi == 1 & UniqueW1 == 1))
	RedundantW1 <- nrow(subset(allpred, Viterbi == 1 & RedundantW1 == 1))
	Rate.UniqueW1 <- UniqueW1/Viterbi
	Rate.RedundantW1 <- RedundantW1/Viterbi
	
	UniqueW2 <- nrow(subset(allpred, Viterbi == 1 & UniqueW2 == 1))
	RedundantW2 <- nrow(subset(allpred, Viterbi == 1 & RedundantW2 == 1))
	Rate.UniqueW2 <- UniqueW2/Viterbi
	Rate.RedundantW2 <- RedundantW2/Viterbi
	
	UniqueW3 <- nrow(subset(allpred, Viterbi == 1 & UniqueW3 == 1))
	RedundantW3 <- nrow(subset(allpred, Viterbi == 1 & RedundantW3 == 1))
	Rate.UniqueW3 <- UniqueW3/Viterbi
	Rate.RedundantW3 <- RedundantW3/Viterbi

	UniqueW4 <- nrow(subset(allpred, Viterbi == 1 & UniqueW4 == 1))
	RedundantW4 <- nrow(subset(allpred, Viterbi == 1 & RedundantW4 == 1))
	Rate.UniqueW4 <- UniqueW4/Viterbi
	Rate.RedundantW4 <- RedundantW4/Viterbi

	UniqueW5 <- nrow(subset(allpred, Viterbi == 1 & UniqueW5 == 1))
	RedundantW5 <- nrow(subset(allpred, Viterbi == 1 & RedundantW5 == 1))
	Rate.UniqueW5 <- UniqueW5/Viterbi
	Rate.RedundantW5 <- RedundantW5/Viterbi
	
	UniqueW6 <- nrow(subset(allpred, Viterbi == 1 & UniqueW6 == 1))
	RedundantW6 <- nrow(subset(allpred, Viterbi == 1 & RedundantW6 == 1))
	Rate.UniqueW6 <- UniqueW6/Viterbi
	Rate.RedundantW6 <- RedundantW6/Viterbi
	
	UniqueW7 <- nrow(subset(allpred, Viterbi == 1 & UniqueW7 == 1))
	RedundantW7 <- nrow(subset(allpred, Viterbi == 1 & RedundantW7 == 1))
	Rate.UniqueW7 <- UniqueW7/Viterbi
	Rate.RedundantW7 <- RedundantW7/Viterbi
	
	UniqueW8 <- nrow(subset(allpred, Viterbi == 1 & UniqueW8 == 1))
	RedundantW8 <- nrow(subset(allpred, Viterbi == 1 & RedundantW8 == 1))
	Rate.UniqueW8 <- UniqueW8/Viterbi
	Rate.RedundantW8 <- RedundantW8/Viterbi

	UniqueW9 <- nrow(subset(allpred, Viterbi == 1 & UniqueW9 == 1))
	RedundantW9 <- nrow(subset(allpred, Viterbi == 1 & RedundantW9 == 1))
	Rate.UniqueW9 <- UniqueW9/Viterbi
	Rate.RedundantW9 <- RedundantW9/Viterbi

	require(ROCR)
	allpred2 <- subset(allpred, P_dyad > 0)

	pred.UniqueW0 <- prediction(predictions = allpred2$P_dyad, 
				labels = allpred2$UniqueW0)
	AUC.UniqueW0 <- performance(pred.UniqueW0, measure = "auc")@y.values[[1]]
	pred.RedundantW0 <- prediction(predictions = allpred2$P_dyad, 
				labels = allpred2$RedundantW0)
	AUC.RedundantW0 <- performance(pred.RedundantW0, measure = "auc")@y.values[[1]]

	perf.UniqueW0 <- performance(pred.UniqueW0, measure = "tpr", x.measure = "fpr")
	perf.RedundantW0 <- performance(pred.RedundantW0, measure = "tpr", x.measure = "fpr")

	obj.name <- paste(prefix, ".", mm.chr, ".perf.UniqueW0", sep = "")
	assign(obj.name, perf.UniqueW0, envir = .GlobalEnv)
	obj.name <- paste(prefix, ".", mm.chr, ".perf.RedundantW0", sep = "")
	assign(obj.name, perf.RedundantW0, envir = .GlobalEnv)

	pred.UniqueW1 <- prediction(predictions = allpred2$P_dyad, 
				labels = allpred2$UniqueW1)
	AUC.UniqueW1 <- performance(pred.UniqueW1, measure = "auc")@y.values[[1]]
	pred.RedundantW1 <- prediction(predictions = allpred2$P_dyad, 
				labels = allpred2$RedundantW1)
	AUC.RedundantW1 <- performance(pred.RedundantW1, measure = "auc")@y.values[[1]]

	pred.UniqueW2 <- prediction(predictions = allpred2$P_dyad, 
				labels = allpred2$UniqueW2)
	AUC.UniqueW2 <- performance(pred.UniqueW2, measure = "auc")@y.values[[1]]
	pred.RedundantW2 <- prediction(predictions = allpred2$P_dyad, 
				labels = allpred2$RedundantW2)
	AUC.RedundantW2 <- performance(pred.RedundantW2, measure = "auc")@y.values[[1]]

	pred.UniqueW3 <- prediction(predictions = allpred2$P_dyad, 
				labels = allpred2$UniqueW3)
	AUC.UniqueW3 <- performance(pred.UniqueW3, measure = "auc")@y.values[[1]]
	pred.RedundantW3 <- prediction(predictions = allpred2$P_dyad, 
				labels = allpred2$RedundantW3)
	AUC.RedundantW3 <- performance(pred.RedundantW3, measure = "auc")@y.values[[1]]

	pred.UniqueW4 <- prediction(predictions = allpred2$P_dyad, 
				labels = allpred2$UniqueW4)
	AUC.UniqueW4 <- performance(pred.UniqueW4, measure = "auc")@y.values[[1]]
	pred.RedundantW4 <- prediction(predictions = allpred2$P_dyad, 
				labels = allpred2$RedundantW4)
	AUC.RedundantW4 <- performance(pred.RedundantW4, measure = "auc")@y.values[[1]]

	pred.UniqueW5 <- prediction(predictions = allpred2$P_dyad, 
				labels = allpred2$UniqueW5)
	AUC.UniqueW5 <- performance(pred.UniqueW5, measure = "auc")@y.values[[1]]
	pred.RedundantW5 <- prediction(predictions = allpred2$P_dyad, 
				labels = allpred2$RedundantW5)
	AUC.RedundantW5 <- performance(pred.RedundantW5, measure = "auc")@y.values[[1]]

	pred.UniqueW6 <- prediction(predictions = allpred2$P_dyad, 
				labels = allpred2$UniqueW6)
	AUC.UniqueW6 <- performance(pred.UniqueW6, measure = "auc")@y.values[[1]]
	pred.RedundantW6 <- prediction(predictions = allpred2$P_dyad, 
				labels = allpred2$RedundantW6)
	AUC.RedundantW6 <- performance(pred.RedundantW6, measure = "auc")@y.values[[1]]

	pred.UniqueW7 <- prediction(predictions = allpred2$P_dyad, 
				labels = allpred2$UniqueW7)
	AUC.UniqueW7 <- performance(pred.UniqueW7, measure = "auc")@y.values[[1]]
	pred.RedundantW7 <- prediction(predictions = allpred2$P_dyad, 
				labels = allpred2$RedundantW7)
	AUC.RedundantW7 <- performance(pred.RedundantW7, measure = "auc")@y.values[[1]]

	pred.UniqueW8 <- prediction(predictions = allpred2$P_dyad, 
				labels = allpred2$UniqueW8)
	AUC.UniqueW8 <- performance(pred.UniqueW8, measure = "auc")@y.values[[1]]
	pred.RedundantW8 <- prediction(predictions = allpred2$P_dyad, 
				labels = allpred2$RedundantW8)
	AUC.RedundantW8 <- performance(pred.RedundantW8, measure = "auc")@y.values[[1]]

	pred.UniqueW9 <- prediction(predictions = allpred2$P_dyad, 
				labels = allpred2$UniqueW9)
	AUC.UniqueW9 <- performance(pred.UniqueW9, measure = "auc")@y.values[[1]]
	pred.RedundantW9 <- prediction(predictions = allpred2$P_dyad, 
				labels = allpred2$RedundantW9)
	AUC.RedundantW9 <- performance(pred.RedundantW9, measure = "auc")@y.values[[1]]


	out <- data.frame(prefix, Package, Target, Model, 
		Genome.size, Num.Unique, Num.Redundant, Predicted, Viterbi, 
		UniqueW0, RedundantW0, Rate.UniqueW0, Rate.RedundantW0, AUC.UniqueW0, AUC.RedundantW0, 
		UniqueW1, RedundantW1, Rate.UniqueW1, Rate.RedundantW1, AUC.UniqueW1, AUC.RedundantW1, 
		UniqueW2, RedundantW2, Rate.UniqueW2, Rate.RedundantW2, AUC.UniqueW2, AUC.RedundantW2, 
		UniqueW3, RedundantW3, Rate.UniqueW3, Rate.RedundantW3, AUC.UniqueW3, AUC.RedundantW3, 
		UniqueW4, RedundantW4, Rate.UniqueW4, Rate.RedundantW4, AUC.UniqueW4, AUC.RedundantW4, 
		UniqueW5, RedundantW5, Rate.UniqueW5, Rate.RedundantW5, AUC.UniqueW5, AUC.RedundantW5, 
		UniqueW6, RedundantW6, Rate.UniqueW6, Rate.RedundantW6, AUC.UniqueW6, AUC.RedundantW6, 
		UniqueW7, RedundantW7, Rate.UniqueW7, Rate.RedundantW7, AUC.UniqueW7, AUC.RedundantW7, 
		UniqueW8, RedundantW8, Rate.UniqueW8, Rate.RedundantW8, AUC.UniqueW8, AUC.RedundantW8, 
		UniqueW9, RedundantW9, Rate.UniqueW9, Rate.RedundantW9, AUC.UniqueW9, AUC.RedundantW9)
	colnames(out) <- c("Prefix", "Package", "Target", "Model", 
		"Genome.size", "Num.Unique", "Num.Redundant", "Predicted", "Viterbi", 
		"UniqueW0", "RedundantW0", "Rate.UniqueW0", "Rate.RedundantW0", "AUC.UniqueW0", "AUC.RedundantW0", 
		"UniqueW1", "RedundantW1", "Rate.UniqueW1", "Rate.RedundantW1", "AUC.UniqueW1", "AUC.RedundantW1", 
		"UniqueW2", "RedundantW2", "Rate.UniqueW2", "Rate.RedundantW2", "AUC.UniqueW2", "AUC.RedundantW2", 
		"UniqueW3", "RedundantW3", "Rate.UniqueW3", "Rate.RedundantW3", "AUC.UniqueW3", "AUC.RedundantW3", 
		"UniqueW4", "RedundantW4", "Rate.UniqueW4", "Rate.RedundantW4", "AUC.UniqueW4", "AUC.RedundantW4", 
		"UniqueW5", "RedundantW5", "Rate.UniqueW5", "Rate.RedundantW5", "AUC.UniqueW5", "AUC.RedundantW5", 
		"UniqueW6", "RedundantW6", "Rate.UniqueW6", "Rate.RedundantW6", "AUC.UniqueW6", "AUC.RedundantW6", 
		"UniqueW7", "RedundantW7", "Rate.UniqueW7", "Rate.RedundantW7", "AUC.UniqueW7", "AUC.RedundantW7", 
		"UniqueW8", "RedundantW8", "Rate.UniqueW8", "Rate.RedundantW8", "AUC.UniqueW8", "AUC.RedundantW8", 
		"UniqueW9", "RedundantW9", "Rate.UniqueW9", "Rate.RedundantW9", "AUC.UniqueW9", "AUC.RedundantW9")
	
	obj.name <- paste(prefix, ".", mm.chr, ".summary", sep = "")
	assign(obj.name, out, envir = .GlobalEnv)
	#return(out)
}
