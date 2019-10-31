wd <- getwd()
setwd("../RData")

load(file = "mm9_chr19_Uni.RData")
# load(file = "mm9_chr19_Red.RData")
load(file = "mm9_chr19_Red30.RData")
setwd(wd)

# > length(mm9_chr19_Red)
# [1] 19965481

# > length(mm9_chr19_Red30)
# [1] 2044747
# > length(mm9_chr19_Uni)
# [1] 249210
