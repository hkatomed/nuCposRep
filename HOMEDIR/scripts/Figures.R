setwd("RData")

#####################################################
# Figure S1
# source: 20190213_linker_length_distribution_1.txt
#####################################################

load(file = "../RData/sysdata_NuPoP.rda")
load(file = "../RData/sysdata_nuCpos.rda")

## Plot the linker length distribution of NuPoP
plot(x = 1:500, y = numeric(length = 500), ylim = c(0, 0.03), type = "n", 
	xlab = "Length (bp)", ylab = "Probability")
lines(x = 1:500, y = Pd[, 2], lwd = 2, col = "black")
lines(x = 1:500, y = Pd[, 7], lwd = 2, col = "green")
lines(x = 1:500, y = Pd[, 9], lwd = 2, col = "purple")
box(lwd = 2)
axis(side = 1, lwd.ticks = 2)
axis(side = 2, lwd.ticks = 2)

## Plot the linker length distribution of nuCpos
plot(x = 1:500, y = numeric(length = 500), ylim = c(0, 0.065), type = "n", 
	xlab = "Length (bp)", ylab = "Probability")
lines(x = 1:500, y = chem.mm9.LinkerDNA.prob_SMA, lwd = 2, col = "black")
lines(x = 1:500, y = nature11142_s2.linker.147.prob_SMA, lwd = 2, col = "green")
lines(x = 1:500, y = sd01.linker.147.prob_SMA, lwd = 2, col = "purple")
box(lwd = 2)
axis(side = 1, lwd.ticks = 2)
axis(side = 2, lwd.ticks = 2)

## Plot the linker length distribution of NuPoP (magnified veiw)
plot(x = 1:100, y = numeric(length = 100), ylim = c(0, 0.03), type = "n", 
	xlab = "Length (bp)", ylab = "Probability", xaxt = "n", yaxt = "n")
lines(x = 1:100, y = Pd[1:100, 2], lwd = 2, col = "black")
lines(x = 1:100, y = Pd[1:100, 7], lwd = 2, col = "green")
lines(x = 1:100, y = Pd[1:100, 9], lwd = 2, col = "purple")
box(lwd = 2)
axis(side = 1, lwd.ticks = 2, at = c(0, 25, 50, 75, 100))
axis(side = 2, lwd.ticks = 2, at = c(0, 0.015, 0.03))

## Plot the linker length distribution of nuCpos (magnified veiw)
plot(x = 1:100, y = numeric(length = 100), ylim = c(0, 0.065), type = "n", 
	xlab = "Length (bp)", ylab = "Probability", xaxt = "n", yaxt = "n")
lines(x = 1:100, y = chem.mm9.LinkerDNA.prob_SMA[1:100], lwd = 2, col = "black")
lines(x = 1:100, y = nature11142_s2.linker.147.prob_SMA[1:100], lwd = 2, col = "green")
lines(x = 1:100, y = sd01.linker.147.prob_SMA[1:100], lwd = 2, col = "purple")
box(lwd = 2)
axis(side = 1, lwd.ticks = 2, at = c(0, 25, 50, 75, 100))
axis(side = 2, lwd.ticks = 2, at = c(0, 0.03, 0.06))


#####################################################
## Figure S2
## source: 20190225_ROC_1.txt
#####################################################

load(file = "../RData/Rate_Redundant_yeasts.RData")
load(file = "../RData/Rate_Unique_yeasts.RData")

## Figure S2, matched-species model (Redundant)
plot(x = 0:9, y = numeric(length = 10), type = "n", ylim = c(0, 0.7), 
	xaxt = "n", xlab = "Matching window (bp)", ylab = "Matching Rate", 
	main = "matched-species model (Redundant)")
lines(x = 0:9, y = Rate.Redundant[Rate.Redundant$Prefix == "NuPoPScSc", 2:11], 
	lwd = 2, lty = 1, col = "blue")
lines(x = 0:9, y = Rate.Redundant[Rate.Redundant$Prefix == "NuPoPSpSp", 2:11], 
	lwd = 2, lty = 2, col = "blue")
lines(x = 0:9, y = Rate.Redundant[Rate.Redundant$Prefix == "nuCposScSc", 2:11], 
	lwd = 2, lty = 1, col = "red")
lines(x = 0:9, y = Rate.Redundant[Rate.Redundant$Prefix == "nuCposSpSp", 2:11], 
	lwd = 2, lty = 2, col = "red")
box(lwd = 2)
axis(side = 1, lwd.ticks = 2, at = 0:9, labels = as.character(0:9))
axis(side = 2, lwd.ticks = 2)

## Figure S2, Cross-species model (Redundant)
plot(x = 0:9, y = numeric(length = 10), type = "n", ylim = c(0, 0.7), 
	xaxt = "n", xlab = "Matching window (bp)", ylab = "Matching Rate", 
	main = "Cross-species model (Redundant)")
lines(x = 0:9, y = Rate.Redundant[Rate.Redundant$Prefix == "NuPoPScSp", 2:11], 
	lwd = 2, lty = 1, col = "blue")
lines(x = 0:9, y = Rate.Redundant[Rate.Redundant$Prefix == "NuPoPSpSc", 2:11], 
	lwd = 2, lty = 2, col = "blue")
lines(x = 0:9, y = Rate.Redundant[Rate.Redundant$Prefix == "nuCposScSp", 2:11], 
	lwd = 2, lty = 1, col = "red")
lines(x = 0:9, y = Rate.Redundant[Rate.Redundant$Prefix == "nuCposSpSc", 2:11], 
	lwd = 2, lty = 2, col = "red")
box(lwd = 2)
axis(side = 1, lwd.ticks = 2, at = 0:9, labels = as.character(0:9))
axis(side = 2, lwd.ticks = 2)

## Figure S2, matched-species model (Unique)
plot(x = 0:9, y = numeric(length = 10), type = "n", ylim = c(0, 0.4), 
	xaxt = "n", xlab = "Matching window (bp)", ylab = "Matching Rate", 
	main = "matched-species model (Unique)")
lines(x = 0:9, y = Rate.Unique[Rate.Unique$Prefix == "NuPoPScSc", 2:11], 
	lwd = 2, lty = 1, col = "blue")
lines(x = 0:9, y = Rate.Unique[Rate.Unique$Prefix == "NuPoPSpSp", 2:11], 
	lwd = 2, lty = 2, col = "blue")
lines(x = 0:9, y = Rate.Unique[Rate.Unique$Prefix == "nuCposScSc", 2:11], 
	lwd = 2, lty = 1, col = "red")
lines(x = 0:9, y = Rate.Unique[Rate.Unique$Prefix == "nuCposSpSp", 2:11], 
	lwd = 2, lty = 2, col = "red")
box(lwd = 2)
axis(side = 1, lwd.ticks = 2, at = 0:9, labels = as.character(0:9))
axis(side = 2, lwd.ticks = 2)

## Figure S2, Cross-species model (Unique)
plot(x = 0:9, y = numeric(length = 10), type = "n", ylim = c(0, 0.4), 
	xaxt = "n", xlab = "Matching window (bp)", ylab = "Matching Rate", 
	main = "Cross-species model (Unique)")
lines(x = 0:9, y = Rate.Unique[Rate.Unique$Prefix == "NuPoPScSp", 2:11], 
	lwd = 2, lty = 1, col = "blue")
lines(x = 0:9, y = Rate.Unique[Rate.Unique$Prefix == "NuPoPSpSc", 2:11], 
	lwd = 2, lty = 2, col = "blue")
lines(x = 0:9, y = Rate.Unique[Rate.Unique$Prefix == "nuCposScSp", 2:11], 
	lwd = 2, lty = 1, col = "red")
lines(x = 0:9, y = Rate.Unique[Rate.Unique$Prefix == "nuCposSpSc", 2:11], 
	lwd = 2, lty = 2, col = "red")
box(lwd = 2)
axis(side = 1, lwd.ticks = 2, at = 0:9, labels = as.character(0:9))
axis(side = 2, lwd.ticks = 2)


###################################################
## Figure 1C and Table S1
## source: 20190225_ROC_1.txt
###################################################

load(file = "../RData/AUC_summary_yeasts.RData")
load(file = "../RData/num_summary_yeasts.RData")

## Figure 1C, matched-species models
COL <- adjustcolor(c("red", "orange", "gray"), alpha.f = 1)
HEIGHT <- t(as.matrix(num.Summary[,c(10,14,15)]))[,c(1,6,3,8)]
colnames(HEIGHT) <- c("Msc", "Msp", "Csc", "Csp")
barplot(height = HEIGHT, 
	main = "matched-species models", 
	ylab = "Number of Viterbi dyads", ylim = c(0, 75000), yaxt = "n", 
	xlab = "Parameters", col = COL, lwd = 2)
box(lwd = 2)
axis(side = 2, lwd.ticks = 2, at = c(0, 25000, 50000, 75000), labels = c("0", "25k", "50k", "75k"))

## Figure 1C, cross-species models
COL <- adjustcolor(c("red", "orange", "gray"), alpha.f = 1)
HEIGHT <- t(as.matrix(num.Summary[,c(10,14,15)]))[,c(2,5,4,7)]
colnames(HEIGHT) <- c("Msc", "Msp", "Csc", "Csp")
barplot(height = HEIGHT, 
	main = "Cross-species models", 
	ylab = "Number of Viterbi dyads", ylim = c(0, 75000), yaxt = "n", 
	xlab = "Parameters", col = COL, lwd = 2)
box(lwd = 2)
axis(side = 2, lwd.ticks = 2, at = c(0, 25000, 50000, 75000), labels = c("0", "25k", "50k", "75k"))


#####################################################
## Figure 1A
## source: 20190225_ROC_1.txt
#####################################################

prefixes <- c("NuPoPScSc", "NuPoPScSp", "nuCposScSc", "nuCposScSp", 
	"NuPoPSpSc", "NuPoPSpSp", "nuCposSpSc", "nuCposSpSp")

for(i in 1:length(prefixes)){
	prefix <- prefixes[i]
	filename <- paste(prefix, "_perf_RedundantW0.RData", sep = "")
	load(file = filename)
	filename <- paste(prefix, "_perf_UniqueW0.RData", sep = "")
	load(file = filename)
}

# Figure 1A NEW, Matched-species models (Target: Sc, Model: Sc)
plot(x = NuPoPScSc.perf.UniqueW0@x.values[[1]], y = NuPoPScSc.perf.UniqueW0@y.values[[1]], 
	type = "n", lwd = 2, xaxt = "n", yaxt = "n", 
	main = "Sc: Matched-species models", 
	xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(a=0, b= 1, lwd = 4, col = "orange")
lines(x = NuPoPScSc.perf.UniqueW0@x.values[[1]], y = NuPoPScSc.perf.UniqueW0@y.values[[1]], 
	lty = 1, lwd = 4, col = "blue")
# lines(x = NuPoPScSc.perf.RedundantW0@x.values[[1]], y = NuPoPScSc.perf.RedundantW0@y.values[[1]], 
# 	lty = 2, lwd = 4, col = "blue")
lines(x = nuCposScSc.perf.UniqueW0@x.values[[1]], y = nuCposScSc.perf.UniqueW0@y.values[[1]], 
	lty = 1, lwd = 4, col = "red")
# lines(x = nuCposScSc.perf.RedundantW0@x.values[[1]], y = nuCposScSc.perf.RedundantW0@y.values[[1]], 
# 	lty = 2, lwd = 4, col = "red")
box(lwd = 4)
axis(side = 1, lwd.ticks = 4, at = seq(0, 1, 0.2), labels = as.character(seq(0, 1, 0.2)))
axis(side = 2, lwd.ticks = 4, at = seq(0, 1, 0.2), labels = as.character(seq(0, 1, 0.2)))

# Figure 1A NEW, Matched-species models (Target: Sp, Model: Sp)
plot(x = NuPoPSpSp.perf.UniqueW0@x.values[[1]], y = NuPoPSpSp.perf.UniqueW0@y.values[[1]], 
	type = "n", lwd = 2, xaxt = "n", yaxt = "n", 
	main = "Sp: Matched-species models", 
	xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(a=0, b= 1, lwd = 4, col = "orange")
lines(x = NuPoPSpSp.perf.UniqueW0@x.values[[1]], y = NuPoPSpSp.perf.UniqueW0@y.values[[1]], 
	lty = 1, lwd = 4, col = "blue")
# lines(x = NuPoPSpSp.perf.RedundantW0@x.values[[1]], y = NuPoPSpSp.perf.RedundantW0@y.values[[1]], 
# 	lty = 2, lwd = 4, col = "blue")
lines(x = nuCposSpSp.perf.UniqueW0@x.values[[1]], y = nuCposSpSp.perf.UniqueW0@y.values[[1]], 
	lty = 1, lwd = 4, col = "red")
# lines(x = nuCposSpSp.perf.RedundantW0@x.values[[1]], y = nuCposSpSp.perf.RedundantW0@y.values[[1]], 
# 	lty = 2, lwd = 4, col = "red")
box(lwd = 4)
axis(side = 1, lwd.ticks = 4, at = seq(0, 1, 0.2), labels = as.character(seq(0, 1, 0.2)))
axis(side = 2, lwd.ticks = 4, at = seq(0, 1, 0.2), labels = as.character(seq(0, 1, 0.2)))

# Figure 1A OLD, Matched-species models (Target: Sc, Model: Sc)
plot(x = NuPoPScSc.perf.UniqueW0@x.values[[1]], y = NuPoPScSc.perf.UniqueW0@y.values[[1]], 
	type = "n", lwd = 2, xaxt = "n", yaxt = "n", 
	main = "Sc: Matched-species models", 
	xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(a=0, b= 1, lwd = 2, col = "orange")
lines(x = NuPoPScSc.perf.UniqueW0@x.values[[1]], y = NuPoPScSc.perf.UniqueW0@y.values[[1]], 
	lty = 1, lwd = 2, col = "blue")
lines(x = NuPoPScSc.perf.RedundantW0@x.values[[1]], y = NuPoPScSc.perf.RedundantW0@y.values[[1]], 
	lty = 2, lwd = 2, col = "blue")
lines(x = nuCposScSc.perf.UniqueW0@x.values[[1]], y = nuCposScSc.perf.UniqueW0@y.values[[1]], 
	lty = 1, lwd = 2, col = "red")
lines(x = nuCposScSc.perf.RedundantW0@x.values[[1]], y = nuCposScSc.perf.RedundantW0@y.values[[1]], 
	lty = 2, lwd = 2, col = "red")
box(lwd = 2)
axis(side = 1, lwd.ticks = 2, at = seq(0, 1, 0.2), labels = as.character(seq(0, 1, 0.2)))
axis(side = 2, lwd.ticks = 2, at = seq(0, 1, 0.2), labels = as.character(seq(0, 1, 0.2)))

# Figure 1A OLD, Cross-species models (Target: Sc, Model: Sp)
plot(x = NuPoPScSp.perf.UniqueW0@x.values[[1]], y = NuPoPScSp.perf.UniqueW0@y.values[[1]], 
	type = "n", lwd = 2, xaxt = "n", yaxt = "n", 
	main = "Sc: Cross-species models", 
	xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(a=0, b= 1, lwd = 2, col = "orange")
lines(x = NuPoPScSp.perf.UniqueW0@x.values[[1]], y = NuPoPScSp.perf.UniqueW0@y.values[[1]], 
	lty = 1, lwd = 2, col = "blue")
lines(x = NuPoPScSp.perf.RedundantW0@x.values[[1]], y = NuPoPScSp.perf.RedundantW0@y.values[[1]], 
	lty = 2, lwd = 2, col = "blue")
lines(x = nuCposScSp.perf.UniqueW0@x.values[[1]], y = nuCposScSp.perf.UniqueW0@y.values[[1]], 
	lty = 1, lwd = 2, col = "red")
lines(x = nuCposScSp.perf.RedundantW0@x.values[[1]], y = nuCposScSp.perf.RedundantW0@y.values[[1]], 
	lty = 2, lwd = 2, col = "red")
box(lwd = 2)
axis(side = 1, lwd.ticks = 2, at = seq(0, 1, 0.2), labels = as.character(seq(0, 1, 0.2)))
axis(side = 2, lwd.ticks = 2, at = seq(0, 1, 0.2), labels = as.character(seq(0, 1, 0.2)))

# Figure 1A OLD, Cross-species models (Target: Sp, Model: Sc)
plot(x = NuPoPSpSc.perf.UniqueW0@x.values[[1]], y = NuPoPSpSc.perf.UniqueW0@y.values[[1]], 
	type = "n", lwd = 2, xaxt = "n", yaxt = "n", 
	main = "Sp: Cross-species models", 
	xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(a=0, b= 1, lwd = 2, col = "orange")
lines(x = NuPoPSpSc.perf.UniqueW0@x.values[[1]], y = NuPoPSpSc.perf.UniqueW0@y.values[[1]], 
	lty = 1, lwd = 2, col = "blue")
lines(x = NuPoPSpSc.perf.RedundantW0@x.values[[1]], y = NuPoPSpSc.perf.RedundantW0@y.values[[1]], 
	lty = 2, lwd = 2, col = "blue")
lines(x = nuCposSpSc.perf.UniqueW0@x.values[[1]], y = nuCposSpSc.perf.UniqueW0@y.values[[1]], 
	lty = 1, lwd = 2, col = "red")
lines(x = nuCposSpSc.perf.RedundantW0@x.values[[1]], y = nuCposSpSc.perf.RedundantW0@y.values[[1]], 
	lty = 2, lwd = 2, col = "red")
box(lwd = 2)
axis(side = 1, lwd.ticks = 2, at = seq(0, 1, 0.2), labels = as.character(seq(0, 1, 0.2)))
axis(side = 2, lwd.ticks = 2, at = seq(0, 1, 0.2), labels = as.character(seq(0, 1, 0.2)))

# Figure 1A OLD, Matched-species models (Target: Sp, Model: Sp)
plot(x = NuPoPSpSp.perf.UniqueW0@x.values[[1]], y = NuPoPSpSp.perf.UniqueW0@y.values[[1]], 
	type = "n", lwd = 2, xaxt = "n", yaxt = "n", 
	main = "Sp: Matched-species models", 
	xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(a=0, b= 1, lwd = 2, col = "orange")
lines(x = NuPoPSpSp.perf.UniqueW0@x.values[[1]], y = NuPoPSpSp.perf.UniqueW0@y.values[[1]], 
	lty = 1, lwd = 2, col = "blue")
lines(x = NuPoPSpSp.perf.RedundantW0@x.values[[1]], y = NuPoPSpSp.perf.RedundantW0@y.values[[1]], 
	lty = 2, lwd = 2, col = "blue")
lines(x = nuCposSpSp.perf.UniqueW0@x.values[[1]], y = nuCposSpSp.perf.UniqueW0@y.values[[1]], 
	lty = 1, lwd = 2, col = "red")
lines(x = nuCposSpSp.perf.RedundantW0@x.values[[1]], y = nuCposSpSp.perf.RedundantW0@y.values[[1]], 
	lty = 2, lwd = 2, col = "red")
box(lwd = 2)
axis(side = 1, lwd.ticks = 2, at = seq(0, 1, 0.2), labels = as.character(seq(0, 1, 0.2)))
axis(side = 2, lwd.ticks = 2, at = seq(0, 1, 0.2), labels = as.character(seq(0, 1, 0.2)))



#####################################################
## Figure 1B
## source: 20190225_ROC_1.txt
#####################################################

load(file = "../RData/nature11142_s3.NuPoPScSc.PredVsChem.plot400.RData")
load(file = "../RData/nature11142_s3.nuCposScSc.PredVsChem.plot400.RData")

ylim <- c(0,0.3)
xlim <- c(0,400)

## Figure 1B, NuPoP
plot(x = seq(0,400,1), y = nature11142_s3.NuPoPScSc.PredVsChem.plot400$prob[1:401], type = "l", 
	main = "nature11142_s3.NuPoPScSc.PredVsChem.plot400", sub = "n = 65,295", 
	ylab = "Probability", 
	xlab = "Distance from Viterbi dyad (bp)", ylim = ylim, xlim = xlim, 
	lwd = 2)
box(lwd = 2)
axis(side = 1, lwd.ticks = 2)
axis(side = 2, lwd.ticks = 2)

## Figure 1B, nuCpos
plot(x = seq(0,400,1), y = nature11142_s3.nuCposScSc.PredVsChem.plot400$prob[1:401], type = "l", 
	main = "nature11142_s3.nuCposScSc.PredVsChem.plot400", sub = "n = 65,295", 
	ylab = "Probability", 
	xlab = "Distance from Viterbi dyad (bp)", ylim = ylim, xlim = xlim, 
	lwd = 2)
box(lwd = 2)
axis(side = 1, lwd.ticks = 2)
axis(side = 2, lwd.ticks = 2)

ylim <- c(0,0.3)
xlim <- c(0,60)

## Figure 1B, NuPoP (magnified view)
plot(x = seq(0,400,1), y = nature11142_s3.NuPoPScSc.PredVsChem.plot400$prob[1:401], type = "l", 
	main = "nature11142_s3.NuPoPScSc.PredVsChem.plot400", sub = "n = 65,295", 
	ylab = "Probability", 
	xlab = "Distance from Viterbi dyad (bp)", ylim = ylim, xlim = xlim, 
	lwd = 2)
box(lwd = 2)
axis(side = 1, lwd.ticks = 2)
axis(side = 2, lwd.ticks = 2)

## Figure 1B, nuCpos (magnified view)
plot(x = seq(0,400,1), y = nature11142_s3.nuCposScSc.PredVsChem.plot400$prob[1:401], type = "l", 
	main = "nature11142_s3.nuCposScSc.PredVsChem.plot400", sub = "n = 65,295", 
	ylab = "Probability", 
	xlab = "Distance from Viterbi dyad (bp)", ylim = ylim, xlim = xlim, 
	lwd = 2)
box(lwd = 2)
axis(side = 1, lwd.ticks = 2)
axis(side = 2, lwd.ticks = 2)


#####################################################
## Figure S9
## source: 20190306_ROC_1.txt
#####################################################

load(file = "../RData/Rate_Redundant_mouse.RData")
load(file = "../RData/Rate_Unique_mouse.RData")

## Figure S9, Mouse model (Redundant)
plot(x = 0:9, y = numeric(length = 10), type = "n", ylim = c(0, 0.7), 
	xaxt = "n", xlab = "Matching window (bp)", ylab = "Matching Rate", 
	main = "Mouse model (Redundant)")
lines(x = 0:9, y = Rate.Redundant[Rate.Redundant$Prefix == "NuPoPMmMm", 2:11], 
	lwd = 2, lty = 1, col = "blue")
lines(x = 0:9, y = Rate.Redundant[Rate.Redundant$Prefix == "nuCposMmMm", 2:11], 
	lwd = 2, lty = 1, col = "red")
box(lwd = 2)
axis(side = 1, lwd.ticks = 2, at = 0:9, labels = as.character(0:9))
axis(side = 2, lwd.ticks = 2)

## Figure S9, Mouse model (Unique)
plot(x = 0:9, y = numeric(length = 10), type = "n", ylim = c(0, 0.2), 
	xaxt = "n", xlab = "Matching window (bp)", ylab = "Matching Rate", 
	main = "Mouse model (Unique)")
lines(x = 0:9, y = Rate.Unique[Rate.Unique$Prefix == "NuPoPMmMm", 2:11], 
	lwd = 2, lty = 1, col = "blue")
lines(x = 0:9, y = Rate.Unique[Rate.Unique$Prefix == "nuCposMmMm", 2:11], 
	lwd = 2, lty = 1, col = "red")
box(lwd = 2)
axis(side = 1, lwd.ticks = 2, at = 0:9, labels = as.character(0:9))
axis(side = 2, lwd.ticks = 2)


#####################################################
## Figure 7D and 7C
## source: 20190306_ROC_1.txt
#####################################################

load(file = "../RData/num_summary_mouse.RData")

## Figure 7D, mouse models
COL <- adjustcolor(c("red", "orange", "gray"), alpha.f = 1)
HEIGHT <- t(as.matrix(num.Summary[,c(10,14,15)]))[,c(1,2)]
colnames(HEIGHT) <- c("Mmm", "Cmm")
barplot(height = HEIGHT, 
	main = "Mouse models", 
	ylab = "Number of Viterbi dyads", ylim = c(0, 300000), yaxt = "n", xlim = c(0.5, 4.2), 
	xlab = "Parameters", col = COL, lwd = 2, space = 0.9)
box(lwd = 2)
axis(side = 2, lwd.ticks = 2, at = c(0, 100000, 200000, 300000), labels = c("0", "100k", "200k", "300k"))


prefixes <- c("NuPoPMmMm", "nuCposMmMm")
mm.chr <- "chr19"
for(i in 1:length(prefixes)){
	prefix <- prefixes[i]
	filename <- paste(prefix, "_", mm.chr, "_perf_RedundantW0.RData", sep = "")
	load(file = filename)
	filename <- paste(prefix, "_", mm.chr, "_perf_UniqueW0.RData", sep = "")
	load(file = filename)
}

## Figure 7C NEW, mouse model
plot(x = NuPoPMmMm.chr19.perf.UniqueW0@x.values[[1]], y = NuPoPMmMm.chr19.perf.UniqueW0@y.values[[1]], 
	type = "n", lwd = 2, xaxt = "n", yaxt = "n", 
	main = "Mouse models", 
	xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(a=0, b= 1, lwd = 4, col = "orange")
lines(x = NuPoPMmMm.chr19.perf.UniqueW0@x.values[[1]], y = NuPoPMmMm.chr19.perf.UniqueW0@y.values[[1]], 
	lty = 1, lwd = 4, col = "blue")
# lines(x = NuPoPMmMm.chr19.perf.RedundantW0@x.values[[1]], y =NuPoPMmMm.chr19.perf.RedundantW0@y.values[[1]], 
# 	lty = 2, lwd = 4, col = "blue")
lines(x = nuCposMmMm.chr19.perf.UniqueW0@x.values[[1]], y = nuCposMmMm.chr19.perf.UniqueW0@y.values[[1]], 
	lty = 1, lwd = 4, col = "red")
# lines(x = nuCposMmMm.chr19.perf.RedundantW0@x.values[[1]], y = nuCposMmMm.chr19.perf.RedundantW0@y.values[[1]], 
# 	lty = 2, lwd = 4, col = "red")
box(lwd = 4)
axis(side = 1, lwd.ticks = 4, at = seq(0, 1, 0.2), labels = as.character(seq(0, 1, 0.2)))
axis(side = 2, lwd.ticks = 4, at = seq(0, 1, 0.2), labels = as.character(seq(0, 1, 0.2)))


## Figure 7C OLD, mouse model
plot(x = NuPoPMmMm.chr19.perf.UniqueW0@x.values[[1]], y = NuPoPMmMm.chr19.perf.UniqueW0@y.values[[1]], 
	type = "n", lwd = 2, xaxt = "n", yaxt = "n", 
	main = "Mouse models", 
	xlab = "False Positive Rate", ylab = "True Positive Rate")
abline(a=0, b= 1, lwd = 2, col = "orange")
lines(x = NuPoPMmMm.chr19.perf.UniqueW0@x.values[[1]], y = NuPoPMmMm.chr19.perf.UniqueW0@y.values[[1]], 
	lty = 1, lwd = 2, col = "blue")
lines(x = NuPoPMmMm.chr19.perf.RedundantW0@x.values[[1]], y =NuPoPMmMm.chr19.perf.RedundantW0@y.values[[1]], 
	lty = 2, lwd = 2, col = "blue")
lines(x = nuCposMmMm.chr19.perf.UniqueW0@x.values[[1]], y = nuCposMmMm.chr19.perf.UniqueW0@y.values[[1]], 
	lty = 1, lwd = 2, col = "red")
lines(x = nuCposMmMm.chr19.perf.RedundantW0@x.values[[1]], y = nuCposMmMm.chr19.perf.RedundantW0@y.values[[1]], 
	lty = 2, lwd = 2, col = "red")
box(lwd = 2)
axis(side = 1, lwd.ticks = 2, at = seq(0, 1, 0.2), labels = as.character(seq(0, 1, 0.2)))
axis(side = 2, lwd.ticks = 2, at = seq(0, 1, 0.2), labels = as.character(seq(0, 1, 0.2)))



####################################################
# Figure 2A
# source: 20190213_Widom_601_HBA_1.txt
####################################################

load(file = "../RData/ori601HBA_m_mean.RData")
load(file = "../RData/ori601HBA_c_mean.RData")

## Figure 2A, MNase-seq-based model
plot(x = 74:209, y = ori601HBA.m.mean$HBA, type = "n", lwd = 2, 
	xlab = "Dyad position (bp)", ylab = "Normalized MNase-seq-based HBA", ylim = c(-10, 10),
	main = "MNase-seq-based HBA scores on the original 601 sequence")
abline(v = 154, lwd = 2, col = "orange")
lines(x = 74:209, y = ori601HBA.m.mean$HBA, type = "l", lwd = 2, col = "black")
box(lwd = 2)
axis(side = 1, lwd.ticks = 2)
axis(side = 2, lwd.ticks = 2)

## Figure 2A, chemical-map-based model
plot(x = 74:209, y = ori601HBA.c.mean$HBA, type = "n", lwd = 2, 
	xlab = "Dyad position (bp)", ylab = "Normalized Chemical-map-based HBA", ylim = c(-10, 10),
	main = "Chemical-map-based HBA scores on the original 601 sequence")
abline(v = 154, lwd = 2, col = "orange")
lines(x = 74:209, y = ori601HBA.c.mean$HBA, type = "l", lwd = 2, col = "black")
box(lwd = 2)
axis(side = 1, lwd.ticks = 2)
axis(side = 2, lwd.ticks = 2)


####################################################
# Figure 2B
# source: 20190213_Widom_601_HBA_1.txt
####################################################

load(file = "../RData/MMTVHBA_m_mean.RData")
load(file = "../RData/MMTVHBA_c_mean.RData")
load(file = "../RData/MMTV.RData")
load(file = "../RData/MMTV_DYAD.RData")

## Figure 2B, MNase-seq-based model
plot(x = 74:(nchar(MMTV)-73), y = MMTVHBA.m.mean$HBA, type = "n", lwd = 2, 
	xlab = "Dyad position (bp)", ylab = "Normalized MNase-seq-based HBA", ylim = c(-15, 20),
	main = "MNase-seq-based HBA scores on the original 601 sequence", yaxt = "n")
abline(v = DYAD, lwd = 2, col = "orange")
lines(x = 74:(nchar(MMTV)-73), y = MMTVHBA.m.mean$HBA, type = "l", lwd = 2, col = "black")
box(lwd = 2)
axis(side = 1, lwd.ticks = 2)
axis(side = 2, lwd.ticks = 2, at = c(-10, 0, 10, 20))

## Figure 2B, chemical-map-based model
plot(x = 74:(nchar(MMTV)-73), y = MMTVHBA.c.mean$HBA, type = "n", lwd = 2, 
	xlab = "Dyad position (bp)", ylab = "Normalized Chemical-map-based HBA", ylim = c(-15, 20),
	main = "Chemical-map-based HBA scores on the original 601 sequence", yaxt = "n")
abline(v = DYAD, lwd = 2, col = "orange")
lines(x = 74:(nchar(MMTV)-73), y = MMTVHBA.c.mean$HBA, type = "l", lwd = 2, col = "black")
box(lwd = 2)
axis(side = 1, lwd.ticks = 2)
axis(side = 2, lwd.ticks = 2, at = c(-10, 0, 10, 20))


####################################################
# Figure 2C
# source: 20190214_5SrDNA_HBA_1.txt
####################################################

load(file = "../RData/rDNAHBA_m_mean.RData")
load(file = "../RData/rDNAHBA_c_mean.RData")
load(file = "../RData/rDNA.RData")
load(file = "../RData/rDNA_DYAD.RData")

## Figure 2C, MNase-seq-based model
plot(x = 74:(nchar(rDNA)-73), y = rDNAHBA.m.mean$HBA, type = "n", lwd = 2, 
	xlab = "Dyad position (bp)", ylab = "Normalized MNase-seq-based HBA", ylim = c(-10, 10),
	main = "MNase-seq-based HBA scores on the original 601 sequence")
abline(v = DYAD, lwd = 2, col = "orange")
lines(x = 74:(nchar(rDNA)-73), y = rDNAHBA.m.mean$HBA, type = "l", lwd = 2, col = "black")
box(lwd = 2)
axis(side = 1, lwd.ticks = 2)
axis(side = 2, lwd.ticks = 2)

## Figure 2C, chemical-map-based model
plot(x = 74:(nchar(rDNA)-73), y = rDNAHBA.c.mean$HBA, type = "n", lwd = 2, 
	xlab = "Dyad position (bp)", ylab = "Normalized Chemical-map-based HBA", ylim = c(-10, 10),
	main = "Chemical-map-based HBA scores on the original 601 sequence")
abline(v = DYAD, lwd = 2, col = "orange")
lines(x = 74:(nchar(rDNA)-73), y = rDNAHBA.c.mean$HBA, type = "l", lwd = 2, col = "black")
box(lwd = 2)
axis(side = 1, lwd.ticks = 2)
axis(side = 2, lwd.ticks = 2)


####################################################
# Figure 3A
# source: 20190205_compare_HBA_1.txt
####################################################

load(file = "../RData/Chereji_plus1_HBAm.RData")
load(file = "../RData/Chereji_plus1_HBAc.RData")
load(file = "../RData/Chereji_minus1_HBAm.RData")
load(file = "../RData/Chereji_minus1_HBAc.RData")

## Figure 3A (+1)
ylim = c(-8, 2)
ylab = "Mean HBA"
main = "+1 nucleosome"
xlab = "Distance from +1 nucleosome dyad"
plot(x = -100:100, y = numeric(length = 201), ylim = ylim, type = "n", xaxt = "n", yaxt = "n", 
	main = main, xlab = xlab, ylab = ylab)
abline(v = 0, lwd = 2, col = "gray")
box(lwd = 2)
axis(side = 1, lwd.ticks = 2)
axis(side = 2, lwd.ticks = 2)
lines(x = -100:100, y = colMeans(Chereji.plus1.HBAm[, 1:201], na.rm = TRUE), lty = 1, lwd = 2, col = "blue")
lines(x = -100:100, y = colMeans(Chereji.plus1.HBAc[, 1:201], na.rm = TRUE), lty = 1, lwd = 2, col = "red")

## Figure 3A (-1)
main = "-1 nucleosome"
xlab = "Distance from -1 nucleosome dyad"
plot(x = -100:100, y = numeric(length = 201), ylim = ylim, type = "n", xaxt = "n", yaxt = "n",
	main = main, xlab = xlab, ylab = ylab)
abline(v = 0, lwd = 2, col = "gray")
box(lwd = 2)
axis(side = 1, lwd.ticks = 2)
axis(side = 2, lwd.ticks = 2)
lines(x = -100:100, y = colMeans(Chereji.minus1.HBAm[, 1:201], na.rm = TRUE), lty = 1, lwd = 2, col = "blue")
lines(x = -100:100, y = colMeans(Chereji.minus1.HBAc[, 1:201], na.rm = TRUE), lty = 1, lwd = 2, col = "red")

####################################################
# Figure 3B
# source: 20190205_compare_HBA_1.txt
####################################################

load(file = "../RData/Chereji_plus1_AT.RData")
load(file = "../RData/Chereji_minus1_AT.RData")

## Figure 3B (+1)
ylim = c(0.6, 0.65)
ylab = "Mean A/T content"
main = "+1 nucleosome"
xlab = "Distance from +1 nucleosome dyad"
plot(x = -100:100, y = numeric(length = 201), ylim = ylim, type = "n", xaxt = "n", yaxt = "n", 
	main = main, xlab = xlab, ylab = ylab)
abline(v = 0, lwd = 2, col = "gray")
box(lwd = 2)
axis(side = 1, lwd.ticks = 2)
axis(side = 2, lwd.ticks = 2)
lines(x = -100:100, y = colMeans(Chereji.plus1.AT[, 1:201], na.rm = TRUE), lty = 1, lwd = 2, col = "black")

## Figure 3B (-1)
main = "-1 nucleosome"
xlab = "Distance from -1 nucleosome dyad"
plot(x = -100:100, y = numeric(length = 201), ylim = ylim, type = "n", xaxt = "n", yaxt = "n",
	main = main, xlab = xlab, ylab = ylab)
abline(v = 0, lwd = 2, col = "gray")
box(lwd = 2)
axis(side = 1, lwd.ticks = 2)
axis(side = 2, lwd.ticks = 2)
lines(x = -100:100, y = colMeans(Chereji.minus1.AT[, 1:201], na.rm = TRUE), lty = 1, lwd = 2, col = "black")


############################################################
# Figure 3C and S3
# 20190716_TRP1ARS1.txt, 20190828_TRP1ARS1_NucII.txt
############################################################

load(file = "../RData/results_nuCpos_Sc_TRP1ARS1.RData")
load(file = "../RData/results_NuPoP_Sc_TRP1ARS1.RData")
load(file = "../RData/TRP1ARS1_AT.RData")
source(file = "../scripts/drawResults2.R")
source(file = "../scripts/drawResults2_NucII.R")

## Figure 3C, nuCpos
drawResults2(TYPE = "Sc", PACKAGE = "nuCpos")

## Figure 3C, NuPoP
drawResults2(TYPE = "Sc", PACKAGE = "NuPoP", ymax.HBA = 25)

## Figure 3C, AT
LEN <- 1465
SEQNAME <- "TRP1ARS1"
RESULTS <- TRP1ARS1.AT[((LEN*2+1)-1000):((LEN*3)+1000)]
	
par(mfrow = c(2, 1))
TITLE <- paste("Seq:", SEQNAME, sep = "")
cat(TITLE, "\n")
	
plot(x = -1000:(LEN+1000-1), y = RESULTS, type = "n", xlim = c(-200, (LEN+199)), 
		ylim = c(0.45, 0.75), xlab = "Position from to the translation initiation site (bp)", 
		ylab = "AT-content", xaxt = "n", yaxt = "n")
title(TITLE)
lines(x = -1000:(LEN+1000-1), y = RESULTS, col = "black", lwd = 2)
lines(x = c(0, LEN), y = c(0.46, 0.46), lwd = 3, col = "gray")
lines(x = c(115, 789), y = c(0.47, 0.47), lwd = 3, col = "red")
lines(x = c(869, 879), y = c(0.48, 0.48), lwd = 4, col = "blue")
lines(x = c(847, 859), y = c(0.48, 0.48), lwd = 4, col = "green4")
lines(x = c(810, 820), y = c(0.48, 0.48), lwd = 4, col = "purple")
lines(x = c(763, 788), y = c(0.48, 0.48), lwd = 4, col = "black")
points(x = c(104, 285, 450, 622, 940, 1132, 1308), 
		y = c(0.72, 0.72, 0.72, 0.72, 0.72, 0.72, 0.72), pch = 25, bg = "orange")
points(x = c(117, 284, 449, 613, 950, 1143, 1329), 
		y = c(0.74, 0.74, 0.74, 0.74, 0.74, 0.74, 0.74), pch = 25, bg = "lightblue")
box(lwd = 2)
axis(1, lwd = 2)
axis(2, at = seq(0, 1, 0.1), lwd = 2)


## Figure S3, nuCpos
drawResults2.NucII(TYPE = "Sc", PACKAGE = "nuCpos")
drawResults2.NucII(TYPE = "Sc", PACKAGE = "NuPoP")


############################################################
# Figure S4
# 20190719_ura4.txt, 20190830_ura4.txt
############################################################

load(file = "../RData/results_NuPoP_ura4.RData")
load(file = "../RData/results_nuCpos_ura4.RData")
load(file = "../RData/ura4_dyads.RData")
load(file = "../RData/seqUra4_AT.RData")

## Figure S4, nuCpos
SOFTWARE <- "nuCpos"
LEN <- 795
SEQNAME <- "seqUra4"
par(mfrow = c(2, 1))
RESULTS <- get(paste("results.", SOFTWARE, sep = ""))
TITLE <- paste("Seq: ura4, software: ", SOFTWARE, sep = "")
plot(x = -1000:(LEN+1000), y = RESULTS[,3], type = "n", xlim = c(-500, (LEN+499)), 
		ylim = c(-0.15, 1.15), xlab = "Position from to the translation initiation site (bp)", 
		ylab = "Prob./Occup.", xaxt = "n", yaxt = "n")
title(TITLE)
polygon(x = c(-1000, -1000:(LEN+1000), (LEN+1000)), 
		y = c(0, RESULTS[,3], 0), col = 8)	# p-dyad
points(x = -1000:(LEN+1000)+73, y = RESULTS[,2], type = "h", col = 4, lwd = 2)
lines(x = c(-151, -1), y = c(-0.1, -0.1), lwd = 3, col = "gray")
lines(x = c((LEN-1+1), (LEN-1+192)), y = c(-0.1, -0.1), lwd = 3, col = "gray")
arrows(x0 = 0, y0 = -0.1, x1 = LEN-1, y1 = -0.1, length = 0.1, lwd = 3, col = "red")
points(x = dyads, 
		y = c(1.1, 1.1, 1.1, 1.1, 1.1, 1.1), pch = 25, bg = "purple")
box(lwd = 2)
axis(1, lwd = 2)
axis(2, at = seq(0, 1, 0.5), lwd = 2)
plot(x = -1000:(LEN+1000), y = RESULTS[,5], type = "n", xlim = c(-500, (LEN+499)), 
		ylim = c(-22, 22), xlab = "Position from to the translation initiation site (bp)", 
		ylab = "HBA", xaxt = "n", yaxt = "n")
title(TITLE)
lines(x = -1000:(LEN+1000), y = RESULTS[,5], col = 1, lwd = 2)
yvalue <- 20
points(x = dyads, 
		y = c(yvalue, yvalue, yvalue, yvalue, yvalue, yvalue), pch = 25, bg = "purple")
yvalue2 <- -20
lines(x = c(-151, -1), y = c(yvalue2, yvalue2), lwd = 3, col = "gray")
lines(x = c((LEN-1+1), (LEN-1+192)), y = c(yvalue2, yvalue2), lwd = 3, col = "gray")
arrows(x0 = 0, y0 = yvalue2, x1 = LEN-1, y1 = yvalue2, length = 0.1, lwd = 3, col = "red")
box(lwd = 2)
axis(1, lwd = 2)
axis(2, at = seq(-20, 20, 10), lwd = 2)

## Figure S4, NuPoP
SOFTWARE <- "NuPoP"
LEN <- 795
SEQNAME <- "seqUra4"
par(mfrow = c(2, 1))
RESULTS <- get(paste("results.", SOFTWARE, sep = ""))
TITLE <- paste("Seq: ura4, software: ", SOFTWARE, sep = "")
plot(x = -1000:(LEN+1000), y = RESULTS[,3], type = "n", xlim = c(-500, (LEN+499)), 
		ylim = c(-0.15, 1.15), xlab = "Position from to the translation initiation site (bp)", 
		ylab = "Prob./Occup.", xaxt = "n", yaxt = "n")
title(TITLE)
polygon(x = c(-1000, -1000:(LEN+1000), (LEN+1000)), 
		y = c(0, RESULTS[,3], 0), col = 8)	# p-dyad
points(x = -1000:(LEN+1000)+73, y = RESULTS[,2], type = "h", col = 4, lwd = 2)
lines(x = c(-151, -1), y = c(-0.1, -0.1), lwd = 3, col = "gray")
lines(x = c((LEN-1+1), (LEN-1+192)), y = c(-0.1, -0.1), lwd = 3, col = "gray")
arrows(x0 = 0, y0 = -0.1, x1 = LEN-1, y1 = -0.1, length = 0.1, lwd = 3, col = "red")
points(x = dyads, 
		y = c(1.1, 1.1, 1.1, 1.1, 1.1, 1.1), pch = 25, bg = "purple")
box(lwd = 2)
axis(1, lwd = 2)
axis(2, at = seq(0, 1, 0.5), lwd = 2)
plot(x = -1000:(LEN+1000), y = RESULTS[,5], type = "n", xlim = c(-500, (LEN+499)), 
		ylim = c(-22, 22), xlab = "Position from to the translation initiation site (bp)", 
		ylab = "HBA", xaxt = "n", yaxt = "n")
title(TITLE)
lines(x = -1000:(LEN+1000), y = RESULTS[,5], col = 1, lwd = 2)
yvalue <- 20
points(x = dyads, 
		y = c(yvalue, yvalue, yvalue, yvalue, yvalue, yvalue), pch = 25, bg = "purple")
yvalue2 <- -20
lines(x = c(-151, -1), y = c(yvalue2, yvalue2), lwd = 3, col = "gray")
lines(x = c((LEN-1+1), (LEN-1+192)), y = c(yvalue2, yvalue2), lwd = 3, col = "gray")
arrows(x0 = 0, y0 = yvalue2, x1 = LEN-1, y1 = yvalue2, length = 0.1, lwd = 3, col = "red")
box(lwd = 2)
axis(1, lwd = 2)
axis(2, at = seq(-20, 20, 10), lwd = 2)

## Figure S4, AT
LEN <- 795
SEQNAME <- "seqUra4"
RESULTS <- seqUra4.AT[2000:4795]
par(mfrow = c(2, 1))
TITLE <- paste("Seq:", SEQNAME, sep = "")
cat(TITLE, "\n")
plot(x = -1000:(LEN+1000), y = RESULTS, type = "n", xlim = c(-500, (LEN+499)), 
		ylim = c(0.45, 0.80), xlab = "Position from to the translation initiation site (bp)", 
		ylab = "AT-content", xaxt = "n", yaxt = "n")
title(TITLE)
lines(x = -1000:(LEN+1000), y = RESULTS, col = "black", lwd = 2)
yvalue <- 0.78
points(x = dyads, 
		y = c(yvalue, yvalue, yvalue, yvalue, yvalue, yvalue), pch = 25, bg = "purple")
yvalue2 <- 0.46
lines(x = c(-151, -1), y = c(yvalue2, yvalue2), lwd = 3, col = "gray")
lines(x = c((LEN-1+1), (LEN-1+192)), y = c(yvalue2, yvalue2), lwd = 3, col = "gray")
arrows(x0 = 0, y0 = yvalue2, x1 = LEN-1, y1 = yvalue2, length = 0.1, lwd = 3, col = "red")
box(lwd = 2)
axis(1, lwd = 2)
axis(2, at = seq(0, 1, 0.1), lwd = 2)
plot(x = -1000:(LEN+1000), y = RESULTS, type = "n")

########################################################
## Figure 4
## source: 20190226_BAR1_3.txt
########################################################

load(file = "../RData/mutNuCpos_BAR1.RData")
source(file = "../scripts/draw.HBA.R")
source(file = "../scripts/draw.occup.R")

ylabels = c(-20, 0, 20)
ylim = c(-30, 20)
draw.HBA("wt", annotation = annotation, ylabels = ylabels, ylim = ylim)
draw.occup("wt", annotation = annotation)

draw.HBA("A20", annotation = annotation, ylabels = ylabels, ylim = ylim)
draw.occup("A20", annotation = annotation)

draw.HBA("A30", annotation = annotation, ylabels = ylabels, ylim = ylim)
draw.occup("A30", annotation = annotation)

draw.HBA("CG4", annotation = annotation, ylabels = ylabels, ylim = ylim)
draw.occup("CG4", annotation = annotation)

draw.HBA("CG5", annotation = annotation, ylabels = ylabels, ylim = ylim)
draw.occup("CG5", annotation = annotation)

draw.HBA("CTG12", annotation = annotation, ylabels = ylabels, ylim = ylim)
draw.occup("CTG12", annotation = annotation)

draw.HBA("Sac5", annotation = annotation, ylabels = ylabels, ylim = ylim)
draw.occup("Sac5", annotation = annotation)

draw.HBA("Sac6", annotation = annotation, ylabels = ylabels, ylim = ylim)
draw.occup("Sac6", annotation = annotation)


########################################################
## Figure S5
## source: 20190802_telo_1.txt
########################################################

load(file = "../RData/mutNuCpos_TALS.RData")
source(file = "../scripts/draw.HBA.R")
source(file = "../scripts/draw.occup.R")

ylabels = c(-20, -10, 0, 10)
ylim = c(-20, 10)
draw.HBA("wt", wt.length = 1811, center = SITE, window = 300, 
	annotation = annotation, ylabels = ylabels, ylim = ylim)
draw.occup("wt", wt.length = 1811, center = SITE, window = 300, annotation = annotation)

draw.HBA("rTTAGGGx29", wt.length = 1811, center = SITE, window = 300, 
	annotation = annotation, ylabels = ylabels, ylim = ylim)
draw.occup("rTTAGGGx29", wt.length = 1811, center = SITE, window = 300, annotation = annotation)

draw.HBA("rTTAGGGx2", wt.length = 1811, center = SITE, window = 300, 
	annotation = annotation, ylabels = ylabels, ylim = ylim)
draw.occup("rTTAGGGx2", wt.length = 1811, center = SITE, window = 300, annotation = annotation)

draw.HBA("rTTAGGGx4", wt.length = 1811, center = SITE, window = 300, 
	annotation = annotation, ylabels = ylabels, ylim = ylim)
draw.occup("rTTAGGGx4", wt.length = 1811, center = SITE, window = 300, annotation = annotation)

draw.HBA("rTTAGGGx6", wt.length = 1811, center = SITE, window = 300, 
	annotation = annotation, ylabels = ylabels, ylim = ylim)
draw.occup("rTTAGGGx6", wt.length = 1811, center = SITE, window = 300, annotation = annotation)

draw.HBA("rTTAGGGx12", wt.length = 1811, center = SITE, window = 300, 
	annotation = annotation, ylabels = ylabels, ylim = ylim)
draw.occup("rTTAGGGx12", wt.length = 1811, center = SITE, window = 300, annotation = annotation)

draw.HBA("rTGTAGGx6", wt.length = 1811, center = SITE, window = 300, 
	annotation = annotation, ylabels = ylabels, ylim = ylim)
draw.occup("rTGTAGGx6", wt.length = 1811, center = SITE, window = 300, annotation = annotation)

draw.HBA("rTGTAGGx12", wt.length = 1811, center = SITE, window = 300, 
	annotation = annotation, ylabels = ylabels, ylim = ylim)
draw.occup("rTGTAGGx12", wt.length = 1811, center = SITE, window = 300, annotation = annotation)

draw.HBA("rTGTGAGx6", wt.length = 1811, center = SITE, window = 300, 
	annotation = annotation, ylabels = ylabels, ylim = ylim)
draw.occup("rTGTGAGx6", wt.length = 1811, center = SITE, window = 300, annotation = annotation)

draw.HBA("rTGTGAGx12", wt.length = 1811, center = SITE, window = 300, 
	annotation = annotation, ylabels = ylabels, ylim = ylim)
draw.occup("rTGTGAGx12", wt.length = 1811, center = SITE, window = 300, annotation = annotation)


#####################################################
## Figure 5A
## source: 20190703_nuCpos_widom601.txt
#####################################################

library(rasterVis)
load(file = "../RData/Widom_lHBA.RData")
levelplot(x = gridNew, useRaster = TRUE, at = seq(-9, 5, 1),
	colorkey = list(at = seq(-10, 6, 0.5), col = COL), col.regions = COL, 
	scales = list(x = list(at = seq(7, 127, 20), labels = seq(80, 200, 20)), 
			y = list(at = seq(gridNum/2, nrow(gridNew), gridNum), 
			labels = LETTERS[13:1])), 
	xlab = "Dyad position (bp)", 
	ylab = "Segments")


#####################################################
## Figure 5B
## source: 20190723_Widom_localHBA_plots.txt
#####################################################

load(file = "../RData/Widom_lHBA.RData")
pos <- 154
ori601df[pos-73,]

main <- paste("pos=", pos, sep = "")
barplot(as.numeric(ori601df[pos-73, 3:15]), 
	xlab = "Segments", ylab = "local HBA", ylim = c(-4, 4), 
	main = main, yaxt = "n", names.arg = LETTERS[1:13])
axis(2, lwd = 2)
abline(h = 0)

# Figure 5B
x <- 1:13
plot(x = 1:13, y = numeric(length = 13), type = "n", xaxt = "n", yaxt = "n", 
	ylab = "local HBA", ylim = c(-10, 5), xlim = c(0.5, 13.5), 
	main = "Positions: 150-158", bty = "n")

for(i in 1:13){
	l <- i + 2
	lines(x = seq(i-0.4, i+0.4, 0.1), y = as.numeric(ori601df[(150-73):(158-73), l]))
}

for(i in 1:4){
	pos <- i + 149	# 150:153
	points(x = x-0.5+0.1*i, y = as.numeric(ori601df[pos-73, 3:15]), 
	pch = 20, col = "black")
}
for(i in 6:9){
	pos <- i + 149	# 155:158
	points(x = x-0.5+0.1*i, y = as.numeric(ori601df[pos-73, 3:15]), 
	pch = 20, col = "black")
}
pos <- 154
	points(x = 1:13, y = as.numeric(ori601df[pos-73, 3:15]), 
	pch = 20, col = "red")
axis(1, lwd = 2, at = 1:13, labels = LETTERS[1:13], tick = FALSE)
axis(2, lwd = 2)

for(i in seq(0.5, 13.5, 1)){
	lines(x = c(i, i), y = c(-10, 5))
}


#####################################################
## Figure S7
## source: 20190725_localHBA_plots_1.txt
#####################################################

source(file = "../scripts/lHBA.rotation.R")
# load(file = "../RData/lHBA_Widom601.RData")
# lHBA.rotation(testSeq, seqName, species = "sc", ylim = c(-10, 5))
load(file = "../RData/lHBA_NCP147.RData")
lHBA.rotation(testSeq, seqName, species = "sc", ylim = c(-10, 5))
load(file = "../RData/lHBA_Widom601LR.RData")
lHBA.rotation(testSeq, seqName, species = "sc", ylim = c(-10, 5))
load(file = "../RData/lHBA_Widom601L.RData")
lHBA.rotation(testSeq, seqName, species = "sc", ylim = c(-10, 5))
load(file = "../RData/lHBA_Widom601R.RData")
lHBA.rotation(testSeq, seqName, species = "sc", ylim = c(-10, 5))


###################################################
## Figures 6A and 6B
## 20190722_nuCpos_ura4.txt
###################################################

# Figure 6A
library(rasterVis)
COL = BuRdTheme()$regions$col

load(file = "../RData/lHBA_ura4_WT.RData")
levelplot(x = gridNew, useRaster = TRUE, at = seq(-9, 5, 1),
	colorkey = list(at = seq(-10, 6, 0.5), col = COL), col.regions = COL, 
	scales = list(x = list(at = seq(1-73, 2000, 250), labels = seq(-500, 1500, 250)), 	
			y = list(at = seq(gridNum/2, nrow(gridNew), gridNum), 
			labels = LETTERS[13:1])), 
	xlab = "Dyad position (bp)", 
	ylab = "Segments")

levelplot(x = gridNew, useRaster = TRUE, at = seq(-9, 5, 1),
	colorkey = list(at = seq(-10, 6, 0.5), col = COL), col.regions = COL, 
	scales = list(x = list(at = dyads+500-73, labels = dyads), 	
			y = list(at = seq(gridNum/2, nrow(gridNew), gridNum), 
			labels = LETTERS[13:1])), 
	xlab = "Dyad position (bp)", 
	ylab = "Segments")


load(file = "../RData/lHBA_ura4_Dyad.RData")
levelplot(x = gridNew, useRaster = TRUE, at = seq(-9, 5, 1),
	colorkey = list(at = seq(-10, 6, 0.5), col = COL), col.regions = COL, 
	# scales = list(x = list(at = seq(1, 2000, 250), labels = seq(-500, 1500, 250)), 	
	scales = list(x = list(at = seq(1-73, 2000, 250), labels = seq(-500, 1500, 250)), 	
			y = list(at = seq(gridNum/2, nrow(gridNew), gridNum), 
			labels = LETTERS[13:1])), 
	xlab = "Dyad position (bp)", 
	ylab = "Segments")

levelplot(x = gridNew, useRaster = TRUE, at = seq(-9, 5, 1),
	colorkey = list(at = seq(-10, 6, 0.5), col = COL), col.regions = COL, 
	# scales = list(x = list(at = dyads+500, labels = dyads), 	
	scales = list(x = list(at = dyads+500-73, labels = dyads), 	
			y = list(at = seq(gridNum/2, nrow(gridNew), gridNum), 
			labels = LETTERS[13:1])), 
	xlab = "Dyad position (bp)", 
	ylab = "Segments")


load(file = "../RData/lHBA_ura4_Linker.RData")
levelplot(x = gridNew, useRaster = TRUE, at = seq(-9, 5, 1),
	colorkey = list(at = seq(-10, 6, 0.5), col = COL), col.regions = COL, 
	# scales = list(x = list(at = seq(1, 2000, 250), labels = seq(-500, 1500, 250)), 	
	scales = list(x = list(at = seq(1-73, 2000, 250), labels = seq(-500, 1500, 250)), 	
			y = list(at = seq(gridNum/2, nrow(gridNew), gridNum), 
			labels = LETTERS[13:1])), 
	xlab = "Dyad position (bp)", 
	ylab = "Segments")

levelplot(x = gridNew, useRaster = TRUE, at = seq(-9, 5, 1),
	colorkey = list(at = seq(-10, 6, 0.5), col = COL), col.regions = COL, 
	# scales = list(x = list(at = dyads+500, labels = dyads), 	
	scales = list(x = list(at = dyads+500-73, labels = dyads), 	
			y = list(at = seq(gridNum/2, nrow(gridNew), gridNum), 
			labels = LETTERS[13:1])), 
	xlab = "Dyad position (bp)", 
	ylab = "Segments")


load(file = "../RData/lHBA_ura4_Int.RData")
levelplot(x = gridNew, useRaster = TRUE, at = seq(-9, 5, 1),
	colorkey = list(at = seq(-10, 6, 0.5), col = COL), col.regions = COL, 
	# scales = list(x = list(at = seq(1, 2000, 250), labels = seq(-500, 1500, 250)), 	
	scales = list(x = list(at = seq(1-73, 2000, 250), labels = seq(-500, 1500, 250)), 	
			y = list(at = seq(gridNum/2, nrow(gridNew), gridNum), 
			labels = LETTERS[13:1])), 
	xlab = "Dyad position (bp)", 
	ylab = "Segments")

levelplot(x = gridNew, useRaster = TRUE, at = seq(-9, 5, 1),
	colorkey = list(at = seq(-10, 6, 0.5), col = COL), col.regions = COL, 
	# scales = list(x = list(at = dyads+500, labels = dyads), 	
	scales = list(x = list(at = dyads+500-73, labels = dyads), 	
			y = list(at = seq(gridNum/2, nrow(gridNew), gridNum), 
			labels = LETTERS[13:1])), 
	xlab = "Dyad position (bp)", 
	ylab = "Segments")

## Figure 6B
load(file = "../RData/ura4_seqs.RData")
par(mfrow = c(3, 2))
for(i in 1:3){
	NUC <- paste("P", i, sep = "")
	df.empty <- data.frame(relpos = -15:15, HBA = numeric(length = 31), 
			lHBA_A = numeric(length = 31), lHBA_B = numeric(length = 31), 
			lHBA_C = numeric(length = 31), lHBA_D = numeric(length = 31), 
			lHBA_E = numeric(length = 31), lHBA_F = numeric(length = 31), 
			lHBA_G = numeric(length = 31), lHBA_H = numeric(length = 31), 
			lHBA_I = numeric(length = 31), lHBA_J = numeric(length = 31), 
			lHBA_K = numeric(length = 31), lHBA_L = numeric(length = 31), 
			lHBA_M = numeric(length = 31))
	df.WT <- df.empty 
	df.Dyad <- df.empty 
	df.Linker <- df.empty 
	df.Int <- df.empty 

	pos <- dyads[i]
	for(j in -15:15){
		relpos <- pos + j
		seq.WT <- substring(seqUra4_0.5k.WT, first = 500+relpos-73, last = 500+relpos+73)
		seq.Dyad <- substring(seqUra4_0.5k.Dyad, first = 500+relpos-73, last = 500+relpos+73)
		seq.Linker <- substring(seqUra4_0.5k.Linker, first = 500+relpos-73, last = 500+relpos+73)
		seq.Int <- substring(seqUra4_0.5k.Int, first = 500+relpos-73, last = 500+relpos+73)
		df.WT$HBA[j+16] <- HBA(seq.WT, species = "sp", silent = TRUE)
		df.Dyad$HBA[j+16] <- HBA(seq.Dyad, species = "sp", silent = TRUE)
		df.Linker$HBA[j+16] <- HBA(seq.Linker, species = "sp", silent = TRUE)
		df.Int$HBA[j+16] <- HBA(seq.Int, species = "sp", silent = TRUE)
		df.WT[j+16,3:15] <- localHBA(seq.WT, species = "sp", silent = TRUE)
		df.Dyad[j+16,3:15] <- localHBA(seq.Dyad, species = "sp", silent = TRUE)
		df.Linker[j+16,3:15] <- localHBA(seq.Linker, species = "sp", silent = TRUE)
		df.Int[j+16,3:15] <- localHBA(seq.Int, species = "sp", silent = TRUE)
	}
	
	plot(x = -15:15, y = df.WT$HBA, xlim = c(-15, 15), ylim = c(-15, 15), 
		xlab = "Dyad position (bp)", ylab = "HBA", type = "n", 
		xaxt = "n", yaxt = "n", )
	lines(x = -15:15, y = df.Dyad$HBA, lwd = 2, col = "green4")
	lines(x = -15:15, y = df.Linker$HBA, lwd = 2, col = "orange")
	lines(x = -15:15, y = df.Int$HBA, lwd = 2, col = "black")
	lines(x = -15:15, y = df.WT$HBA, lwd = 2, col = "gray")
	title(paste(NUC, ", center=", pos, sep = ""))
	box(lwd = 2)
	axis(side = 1, lwd.ticks = 2, at = seq(-15, 15, 5))
	axis(side = 2, lwd.ticks = 2, at = seq(-10, 10, 10))
	
	j <- 0
	plot(x = 1:13, as.numeric(df.WT[j+16,3:15]), type = "n", 
		xaxt = "n", yaxt = "n", ylim = c(-10, 7.5), xlab = "Segments", ylab = "local HBA")
		lines(x = 1:13, as.numeric(df.Dyad[j+16,3:15]), lwd = 2, col = "green4") # Dyad
		lines(x = 1:13, as.numeric(df.Linker[j+16,3:15]), lwd = 2, col = "orange")  # Linker
		lines(x = 1:13, as.numeric(df.Int[j+16,3:15]), lwd = 2, col = "black") # Int
		lines(x = 1:13, as.numeric(df.WT[j+16,3:15]), lwd = 2, col = "gray")
		title(paste("+", i, ", center=", pos, ", relpos=", j, sep = ""))
		axis(1, at = 1:13, labels = LETTERS[1:13], lwd = 2)
		axis(2, lwd = 2)
		box(lwd = 2)
}


##################################################
## Figure S8A
## source: 20190305_TwoLetterFreq_1.txt
##################################################

library(rasterVis)
COL <- BuRdTheme()$regions$col

load(file = "../RData/dinuc.RData")
Sp.two.inv <- Sp.two
for(i in 1:ncol(Sp.two)) Sp.two.inv[,ncol(Sp.two)-i+1] <- Sp.two[,i]
x <- seq(1, nrow(Sp.two.inv), 1)
y <- seq(1, ncol(Sp.two.inv), 1)
Sc.grid <- expand.grid(x=x, y=y)
Sc.grid$z <- as.vector(log10(Sp.two.inv))

Sc.two.inv <- Sc.two
for(i in 1:ncol(Sc.two)) Sc.two.inv[,ncol(Sc.two)-i+1] <- Sc.two[,i]
x <- seq(1, nrow(Sc.two.inv), 1)
y <- seq(1, ncol(Sc.two.inv), 1)
Sp.grid <- expand.grid(x=x, y=y)
Sp.grid$z <- as.vector(log10(Sp.two.inv))

Mm.two.inv <- Mm.two	
for(i in 1:ncol(Mm.two)) Mm.two.inv[,ncol(Mm.two)-i+1] <- Mm.two[,i]
x <- seq(1, nrow(Mm.two.inv), 1)
y <- seq(1, ncol(Mm.two.inv), 1)
Mm.grid <- expand.grid(x=x, y=y)
Mm.grid$z <- as.vector(log10(Mm.two.inv))


max.value <- max(Sc.grid$z, Sp.grid$z, Mm.grid$z, na.rm = TRUE)
min.value <- min(Sc.grid$z, Sp.grid$z, Mm.grid$z, na.rm = TRUE)
my.at <- seq(min.value, max.value, length.out = length(COL) -1)
my.ckey <- list(at = my.at, col = COL)
ylabels <- c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", 
		"GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")
ylabels[1:length(ylabels)] <- ylabels[length(ylabels):1]

levelplot(z~x*y, Sp.grid, useRaster = FALSE, at = my.at, colorkey = my.ckey, col.regions = COL, 
	scales = list(x = list(at = c(1, seq(20, 140, 20)), labels = TRUE), 
			y = list(at = seq(1, 16, 1) - 8 + 8), labels = ylabels), 
	xlab = "Start Position of 2-letter in Nucleosome", 
	ylab = "2-letter", main = "2-letter frequency in Sp nucleosome (log10)")

levelplot(z~x*y, Sc.grid, useRaster = FALSE, at = my.at, colorkey = my.ckey, col.regions = COL, 
	scales = list(x = list(at = c(1, seq(20, 140, 20)), labels = TRUE), 
			y = list(at = seq(1, 16, 1) - 8 + 8), labels = ylabels), 
	xlab = "Start Position of 2-letter in Nucleosome", 
	ylab = "2-letter", main = "2-letter frequency in Sc nucleosome (log10)")

levelplot(z~x*y, Mm.grid, useRaster = FALSE, at = my.at, colorkey = my.ckey, col.regions = COL, 
	scales = list(x = list(at = c(1, seq(20, 140, 20)), labels = TRUE), 
			y = list(at = seq(1, 16, 1) - 8 + 8), labels = ylabels), 
	xlab = "Start Position of 2-letter in Nucleosome", 
	ylab = "2-letter", main = "2-letter frequency in Mm nucleosome (log10)")


##################################################
## Figure S8B
## source: 20190830_SpS_contents.txt
##################################################

load(file = "../RData/dinuc.RData")
source(file = "../scripts/plotDinuc.R")
Sc.two.df <- data.frame(Sc.two)
Sp.two.df <- data.frame(Sp.two)
Mm.two.df <- data.frame(Mm.two)
plotDinuc(species = "sc")
plotDinuc(species = "sp")
plotDinuc(species = "mm")


##################################################
## Figure 7B
## source: 20190305_TwoLetterFreq_1.txt
##################################################

load(file = "../RData/dinucNorm.RData")
source(file = "../scripts/plotDinucNorm.R")

par(mfrow = c(4, 1))
dinucs <- c("CG")
for(i in 1:length(dinucs)) plotDinucNorm(dinuc = dinucs[i])


##################################################
## Figure 7A
## source: 20190305_TwoLetterFreq_1.txt
##################################################

library(rasterVis)
COL <- BuRdTheme()$regions$col
load(file = "../RData/dinucNorm2.RData")
max.value <- max(Sc.grid$z, Sp.grid$z, Mm.grid$z, na.rm = TRUE)
min.value <- min(Sc.grid$z, Sp.grid$z, Mm.grid$z, na.rm = TRUE)
my.at <- seq(min.value, max.value, length.out = length(COL) -1)
my.ckey <- list(at = my.at, col = COL)
ylabels <- c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", 
		"GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")
ylabels[1:length(ylabels)] <- ylabels[length(ylabels):1]
	
levelplot(z~x*y, Sp.grid, useRaster = FALSE, at = my.at, colorkey = my.ckey, col.regions = COL, 
	scales = list(x = list(at = c(1, seq(20, 140, 20)), labels = TRUE), 
			y = list(at = seq(1, 16, 1) - 8 + 8), labels = ylabels), 
	xlab = "Start Position of 2-letter in Nucleosome", 
	ylab = "2-letter", main = "Normalized 2-letter frequency in Sp")

levelplot(z~x*y, Sc.grid, useRaster = FALSE, at = my.at, colorkey = my.ckey, col.regions = COL, 
	scales = list(x = list(at = c(1, seq(20, 140, 20)), labels = TRUE), 
			y = list(at = seq(1, 16, 1) - 8 + 8), labels = ylabels), 
	xlab = "Start Position of 2-letter in Nucleosome", 
	ylab = "2-letter", main = "Normalized 2-letter frequency in Sc")

levelplot(z~x*y, Mm.grid, useRaster = FALSE, at = my.at, colorkey = my.ckey, col.regions = COL, 
	scales = list(x = list(at = c(1, seq(20, 140, 20)), labels = TRUE), 
			y = list(at = seq(1, 16, 1) - 8 + 8), labels = ylabels), 
	xlab = "Start Position of 2-letter in Nucleosome", 
	ylab = "2-letter", main = "Normalized 2-letter frequency in Mm")



