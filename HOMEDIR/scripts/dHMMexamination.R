setwd("RData")

prefixes <- c("NuPoPScSc", "NuPoPScSp", "nuCposScSc", "nuCposScSp", 
	"NuPoPSpSc", "NuPoPSpSp", "nuCposSpSc", "nuCposSpSp")

for(i in 1:length(prefixes)){
	prefix <- prefixes[i]
	filename <- paste(prefix, "_summary.RData", sep = "")
	load(file = filename)
}

Rate.Redundant <- data.frame(Prefix = prefixes, W0 = numeric(length(prefixes)), 
			W1 = numeric(length(prefixes)), W2 = numeric(length(prefixes)), 
			W3 = numeric(length(prefixes)), W4 = numeric(length(prefixes)), 
			W5 = numeric(length(prefixes)), W6 = numeric(length(prefixes)), 
			W7 = numeric(length(prefixes)), W8 = numeric(length(prefixes)), 
			W9 = numeric(length(prefixes)))
for(i in 1:length(prefixes)){
	prefix <- prefixes[i]
	object.name <- paste(prefix, ".summary", sep = "")
	Rate.Redundant[i, 2:11] <- get(object.name)[1,(12+seq(1, 55, 6))]
}

Rate.Unique <- data.frame(Prefix = prefixes, W0 = numeric(length(prefixes)), 
			W1 = numeric(length(prefixes)), W2 = numeric(length(prefixes)), 
			W3 = numeric(length(prefixes)), W4 = numeric(length(prefixes)), 
			W5 = numeric(length(prefixes)), W6 = numeric(length(prefixes)), 
			W7 = numeric(length(prefixes)), W8 = numeric(length(prefixes)), 
			W9 = numeric(length(prefixes)))
for(i in 1:length(prefixes)){
	prefix <- prefixes[i]
	object.name <- paste(prefix, ".summary", sep = "")
	Rate.Unique[i, 2:11] <- get(object.name)[1,(11+seq(1, 55, 6))]
}

AUC.summary <- data.frame(Prefix = prefixes, Package = character(length(prefixes)), 
			Target = character(length(prefixes)), Model = character(length(prefixes)), 
			AUC.UniqueW0 = numeric(length(prefixes)), 
			AUC.RedundantW0 = numeric(length(prefixes)), stringsAsFactors = FALSE)
for(i in 1:length(prefixes)){
	prefix <- prefixes[i]
	object.name <- paste(prefix, ".summary", sep = "")
	for(j in 1:4)	AUC.summary[i, j] <- as.character(get(object.name)[1,j])
	AUC.summary[i, 5:6] <- get(object.name)[1,14:15]
}

# > AUC.summary
#       Prefix Package Target Model AUC.UniqueW0 AUC.RedundantW0
# 1  NuPoPScSc   NuPoP     Sc    Sc    0.5318071       0.5263059
# 2  NuPoPScSp   NuPoP     Sc    Sp    0.5328873       0.5277299
# 3 nuCposScSc  nuCpos     Sc    Sc    0.7506980       0.6745802
# 4 nuCposScSp  nuCpos     Sc    Sp    0.6623388       0.6313702
# 5  NuPoPSpSc   NuPoP     Sp    Sc    0.5255519       0.5149591
# 6  NuPoPSpSp   NuPoP     Sp    Sp    0.5258837       0.5152085
# 7 nuCposSpSc  nuCpos     Sp    Sc    0.6668910       0.6255247
# 8 nuCposSpSp  nuCpos     Sp    Sp    0.7432545       0.6674207

num.Summary <- data.frame(species = c("sc", "sc", "sc", "sc", "sp", "sp", "sp", "sp"), 
			model = c("NuPoP_sc", "NuPoP_sp", "nuCpos_sc", "nuCpos_sp", 
				"NuPoP_sc", "NuPoP_sp", "nuCpos_sc", "nuCpos_sp"), 
			AllPred = integer(8), Viterbi = integer(8), 
			NonRedVit = integer(8), RedVit = integer(8))
	prefix <- prefixes[i]
	object.name <- paste(prefix, ".summary", sep = "")

num.Summary <- data.frame(matrix(nrow = length(prefixes), ncol = 13))
colnames(num.Summary)[1:9] <- colnames(nuCposScSc.summary)[1:9]
colnames(num.Summary)[10:13] <- colnames(nuCposScSc.summary)[16:19]
for(i in 1:length(prefixes)){
	prefix <- prefixes[i]
	object.name <- paste(prefix, ".summary", sep = "")
	for(j in 1:4)	num.Summary[i, j] <- as.character(get(object.name)[1,j])
	num.Summary[i,5:9] <- get(object.name)[1,5:9]
	num.Summary[i,10:13] <- get(object.name)[1,16:19]
}

num.Summary$RedOnlyW1 <- num.Summary$RedundantW1 - num.Summary$UniqueW1
num.Summary$OthersW1 <- num.Summary$Viterbi - num.Summary$UniqueW1 - num.Summary$RedOnlyW1

# > num.Summary
#       Prefix Package Target Model Genome.size Num.Unique Num.Redundant Predicted Viterbi
# 1  NuPoPScSc   NuPoP     Sc    Sc    12071326      67548        344709   2576532   52945
# 2  NuPoPScSp   NuPoP     Sc    Sp    12071326      67548        344709   2720517   58626
# 3 nuCposScSc  nuCpos     Sc    Sc    12071326      67548        344709   2890137   65293
# 4 nuCposScSp  nuCpos     Sc    Sp    12071326      67548        344709   2933607   60683
# 5  NuPoPSpSc   NuPoP     Sp    Sc    12571820      75828        425653   2493027   47876
# 6  NuPoPSpSp   NuPoP     Sp    Sp    12571820      75828        425653   2689223   53657
# 7 nuCposSpSc  nuCpos     Sp    Sc    12571820      75828        425653   3265303   68851
# 8 nuCposSpSp  nuCpos     Sp    Sp    12571820      75828        425653   2888281   67671
#   UniqueW1 RedundantW1 Rate.UniqueW1 Rate.RedundantW1 RedOnlyW1 OthersW1
# 1     1673        6663    0.03159883        0.1258476      4990    46282
# 2     1833        7462    0.03126599        0.1272814      5629    51164
# 3    15327       33241    0.23474186        0.5091051     17914    32052
# 4     7238       22261    0.11927558        0.3668408     15023    38422
# 5     1326        5186    0.02769655        0.1083215      3860    42690
# 6     1454        5915    0.02709805        0.1102372      4461    47742
# 7     9356       25558    0.13588764        0.3712074     16202    43293
# 8    17247       36689    0.25486545        0.5421673     19442    30982

prefixes <- c("NuPoPScSc", "NuPoPScSp", "nuCposScSc", "nuCposScSp", 
	"NuPoPSpSc", "NuPoPSpSp", "nuCposSpSc", "nuCposSpSp")

for(i in 1:length(prefixes)){
	prefix <- prefixes[i]
	if(length(grep("ScSc", prefix)) == 1 )	Target <- "Sc"
	if(length(grep("ScSp", prefix)) == 1)	Target <- "Sc"
	if(length(grep("SpSc", prefix)) == 1)	Target <- "Sp"
	if(length(grep("SpSp", prefix)) == 1)	Target <- "Sp"
	if(Target == "Sc")	file.name <- paste("nature11142_s2_", prefix, ".RData", sep = "")
	if(Target == "Sp")	file.name <- paste("sd01_", prefix, ".RData", sep = "")
	load(file = file.name)
	if(Target == "Sc")	file.name <- paste("nature11142_s3_", prefix, ".RData", sep = "")
	if(Target == "Sp")	file.name <- paste("sd02_", prefix, ".RData", sep = "")
	load(file = file.name)
}

get.count400 <- function(in.dist){
	out <- data.frame(pos = 0:401, count = integer(length = 402), prob = integer(length = 402))
	for(i in 1:401){
		out$count[i] <- length(in.dist[abs(in.dist) == (i-1)])
	}
	out$count[402] <- length(in.dist) - sum(out$count)
	out$prob <- out$count / sum(out$count)
	return(out)
}

nature11142_s3.NuPoPScSc.PredVsChem.plot400 <- get.count400(nature11142_s3.NuPoPScSc[["PredVsChem"]]$dist)
nature11142_s3.nuCposScSc.PredVsChem.plot400 <- get.count400(nature11142_s3.nuCposScSc[["PredVsChem"]]$dist)

save(Rate.Redundant, file = "Rate_Redundant_yeasts.RData")
save(Rate.Unique, file = "Rate_Unique_yeasts.RData")
save(AUC.summary, file = "AUC_summary_yeasts.RData")
save(num.Summary, file = "num_summary_yeasts.RData")
save(nature11142_s3.NuPoPScSc.PredVsChem.plot400, file = "nature11142_s3.NuPoPScSc.PredVsChem.plot400.RData")
save(nature11142_s3.nuCposScSc.PredVsChem.plot400, file = "nature11142_s3.nuCposScSc.PredVsChem.plot400.RData")


prefixes <- c("NuPoPMmMm", "nuCposMmMm")
mm.chr <- "chr19"

for(i in 1:length(prefixes)){
	prefix <- prefixes[i]
	filename <- paste(prefix, "_", mm.chr, "_summary.RData", sep = "")
	load(file = filename)
}

Rate.Redundant <- data.frame(Prefix = prefixes, W0 = numeric(length(prefixes)), 
			W1 = numeric(length(prefixes)), W2 = numeric(length(prefixes)), 
			W3 = numeric(length(prefixes)), W4 = numeric(length(prefixes)), 
			W5 = numeric(length(prefixes)), W6 = numeric(length(prefixes)), 
			W7 = numeric(length(prefixes)), W8 = numeric(length(prefixes)), 
			W9 = numeric(length(prefixes)))
for(i in 1:length(prefixes)){
	prefix <- prefixes[i]
	object.name <- paste(prefix, ".", mm.chr, ".summary", sep = "")
	Rate.Redundant[i, 2:11] <- get(object.name)[1,(12+seq(1, 55, 6))]
}

Rate.Unique <- data.frame(Prefix = prefixes, W0 = numeric(length(prefixes)), 
			W1 = numeric(length(prefixes)), W2 = numeric(length(prefixes)), 
			W3 = numeric(length(prefixes)), W4 = numeric(length(prefixes)), 
			W5 = numeric(length(prefixes)), W6 = numeric(length(prefixes)), 
			W7 = numeric(length(prefixes)), W8 = numeric(length(prefixes)), 
			W9 = numeric(length(prefixes)))
for(i in 1:length(prefixes)){
	prefix <- prefixes[i]
	object.name <- paste(prefix, ".", mm.chr, ".summary", sep = "")
	Rate.Unique[i, 2:11] <- get(object.name)[1,(11+seq(1, 55, 6))]
}

AUC.summary <- data.frame(Prefix = prefixes, Package = character(length(prefixes)), 
			Target = character(length(prefixes)), Model = character(length(prefixes)), 
			AUC.UniqueW0 = numeric(length(prefixes)), 
			AUC.RedundantW0 = numeric(length(prefixes)), stringsAsFactors = FALSE)
for(i in 1:length(prefixes)){
	prefix <- prefixes[i]
	object.name <- paste(prefix, ".", mm.chr, ".summary", sep = "")
	for(j in 1:4)	AUC.summary[i, j] <- as.character(get(object.name)[1,j])
	AUC.summary[i, 5:6] <- get(object.name)[1,14:15]
}

# > AUC.summary
#       Prefix Package Target Model AUC.UniqueW0 AUC.RedundantW0
# 1  NuPoPMmMm   NuPoP     Mm    Mm    0.5106179       0.5085794
# 2 nuCposMmMm  nuCpos     Mm    Mm    0.7107004       0.6456698

num.Summary <- data.frame(species = c("mm", "mm"), 
			model = c("NuPoP_mm", "nuCpos_mm"), 
			AllPred = integer(2), Viterbi = integer(2), 
			NonRedVit = integer(2), RedVit = integer(2))

num.Summary <- data.frame(matrix(nrow = length(prefixes), ncol = 13))
colnames(num.Summary)[1:9] <- colnames(nuCposMmMm.chr19.summary)[1:9]
colnames(num.Summary)[10:13] <- colnames(nuCposMmMm.chr19.summary)[16:19]
for(i in 1:length(prefixes)){
	prefix <- prefixes[i]
	object.name <- paste(prefix, ".", mm.chr, ".summary", sep = "")
	for(j in 1:4)	num.Summary[i, j] <- as.character(get(object.name)[1,j])
	num.Summary[i,5:9] <- get(object.name)[1,5:9]
	num.Summary[i,10:13] <- get(object.name)[1,16:19]
}

num.Summary$RedOnlyW1 <- num.Summary$RedundantW1 - num.Summary$UniqueW1
num.Summary$OthersW1 <- num.Summary$Viterbi - num.Summary$UniqueW1 - num.Summary$RedOnlyW1

# > num.Summary
#       Prefix Package Target Model Genome.size Num.Unique Num.Redundant Predicted Viterbi
# 1  NuPoPMmMm   NuPoP     Mm    Mm    61342430     249210       2044747  11488236  266308
# 2 nuCposMmMm  nuCpos     Mm    Mm    61342430     249210       2044747  23005219  268503
#   UniqueW1 RedundantW1 Rate.UniqueW1 Rate.RedundantW1 RedOnlyW1 OthersW1
# 1     3823       29426    0.01435556        0.1104961     25603   236882
# 2    26803      110239    0.09982384        0.4105690     83436   158264

save(Rate.Redundant, file = "Rate_Redundant_mouse.RData")
save(Rate.Unique, file = "Rate_Unique_mouse.RData")
save(AUC.summary, file = "AUC_summary_mouse.RData")
save(num.Summary, file = "num_summary_mouse.RData")
