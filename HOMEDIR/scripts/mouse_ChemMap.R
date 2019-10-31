setwd("mm_chemicalmap")

# > list.files()
# [1] "chr19_Chemical_NCPscore.txt" "chr19_unique.map_95pc.txt"  

mm9_chr19_Red30 <- read.table(file = "chr19_Chemical_NCPscore.txt")[, c(2,3)]
mm9_chr19_Uni <- read.table(file = "chr19_unique.map_95pc.txt")[, 2]
mm9_chr19_Red <- mm9_chr19_Red30[,1]

# > str(mm9_chr19_Red)
#  int [1:19965481] 3000003 3000004 3000008 3000012 3000014 3000015 3000019 3000022 3000024 3000025 ...
# > str(mm9_chr19_Uni)
#  int [1:249210] 3000090 3000476 3000718 3001092 3002180 3002624 3003038 3003385 3004079 3004274 ...

# > nrow(mm9_chr19_Red30)
# [1] 19965481
# > length(mm9_chr19_Red)
# [1] 19965481
# > length(mm9_chr19_Uni)
# [1] 249210

TOTAL <- 61342431	# chr19 total length
PERbp <- 30		# nucleosomes per PERbp bp

mm9_chr19_Red30 <- mm9_chr19_Red30[order(mm9_chr19_Red30[,2], decreasing = TRUE),]

# > head(mm9_chr19_Red30, n = 3)
#                V2     V3
# 19854638 61000660 35.047
# 19854645 61000672 26.513
# 14731392 45724748 25.277

# > tail(mm9_chr19_Red30, n = 3)
#                V2    V3
# 19964247 61339444 0.001
# 19964766 61340701 0.001
# 19965271 61341864 0.001

# > nrow(mm9_chr19_Red30)
# [1] 19965481

mm9_chr19_Red30 <- mm9_chr19_Red30[1:(TOTAL/PERbp),]

# > nrow(mm9_chr19_Red30)
# [1] 2044747

mm9_chr19_Red30 <- mm9_chr19_Red30[order(mm9_chr19_Red30[,1], decreasing = FALSE),]

# > head(mm9_chr19_Red30, n = 3)
#          V2    V3
# 33  3000090 0.426
# 232 3000636 0.490
# 261 3000718 0.497

# > tail(mm9_chr19_Red30, n = 3)
#                V2    V3
# 19965454 61342356 0.516
# 19965457 61342363 0.530
# 19965458 61342368 0.588

mm9_chr19_Red30 <- mm9_chr19_Red30[,1]

# > length(mm9_chr19_Red30)
# [1] 2044747


setwd("../RData")
save(mm9_chr19_Red30, file = "mm9_chr19_Red30.RData")
save(mm9_chr19_Red, file = "mm9_chr19_Red.RData")
save(mm9_chr19_Uni, file = "mm9_chr19_Uni.RData")
