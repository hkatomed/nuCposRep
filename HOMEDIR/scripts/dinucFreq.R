library(Biostrings)

setwd("RData")
load(file = "nature11142_s2_147.RData")
load(file = "sd01_147.RData")
load(file = "chem_mm9_NucDNA.RData")

nature11142_s2.147 <- DNAStringSet(nature11142_s2.147$seq)
sd01.147 <- DNAStringSet(sd01.147$seq)

Nindex <- letterFrequency(sd01.147, letters = "N") > 0

nature11142_s2.147.noNA <- nature11142_s2.147
sd01.147.noNA <- sd01.147[which(Nindex == FALSE)]

Sc.two <- matrix(ncol = 4^2, nrow = 147-1)
for(i in 1:(147-1)){
	Sc.two[i,] <- oligonucleotideFrequency(
		subseq(nature11142_s2.147.noNA, start = i, end = i + 1), 
		width = 2, as.prob = TRUE, simplify.as = "collapsed", as.array = FALSE)
}

Sp.two <- matrix(ncol = 4^2, nrow = 147-1)
for(i in 1:(147-1)){
	Sp.two[i,] <- oligonucleotideFrequency(
		subseq(sd01.147.noNA, start = i, end = i + 1), 
		width = 2, as.prob = TRUE, simplify.as = "collapsed", as.array = FALSE)
}

Mm.two <- matrix(ncol = 4^2, nrow = 147-1)
for(i in 1:(147-1)){
	Mm.two[i,] <- oligonucleotideFrequency(
		subseq(chem.mm9.NucDNA, start = i, end = i + 1), 
		width = 2, as.prob = TRUE, simplify.as = "collapsed", as.array = FALSE)
}

colnames(Sc.two) <- c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", 
			"GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")

colnames(Sp.two) <- c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", 
			"GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")

colnames(Mm.two) <- c("AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", 
			"GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT")

save(list = c("Sc.two", "Sp.two", "Mm.two"), file = "dinuc.RData")


load(file = "sc_genome_2011.RData")
load(file = "sp_genome.RData")
load(file = "mm9_unmasked.RData")

Sc.two.genome <- oligonucleotideFrequency(sc.genome.2011, 
		width = 2, as.prob = TRUE, simplify.as = "collapsed", as.array = FALSE)

Sp.two.genome <- oligonucleotideFrequency(sp.genome, 
		width = 2, as.prob = TRUE, simplify.as = "collapsed", as.array = FALSE)

Mm.two.genome <- oligonucleotideFrequency(mm9.unmasked, 
		width = 2, as.prob = TRUE, simplify.as = "collapsed", as.array = FALSE)

Sc.twoNorm <- Sc.two
for(i in 1:146)	Sc.twoNorm[i,1:16] <- Sc.twoNorm[i,1:16] / Sc.two.genome

Sp.twoNorm <- Sp.two
for(i in 1:146)	Sp.twoNorm[i,1:16] <- Sp.twoNorm[i,1:16] / Sp.two.genome

Mm.twoNorm <- Mm.two
for(i in 1:146)	Mm.twoNorm[i,1:16] <- Mm.twoNorm[i,1:16] / Mm.two.genome

Sc.twoNorm.df <- data.frame(Sc.twoNorm)
Sp.twoNorm.df <- data.frame(Sp.twoNorm)
Mm.twoNorm.df <- data.frame(Mm.twoNorm)

save(list = c("Sc.twoNorm.df", "Sp.twoNorm.df", "Mm.twoNorm.df"), file = "dinucNorm.RData")


Sp.twoNorm.inv <- Sp.twoNorm
for(i in 1:ncol(Sp.twoNorm)) Sp.twoNorm.inv[,ncol(Sp.twoNorm)-i+1] <- Sp.twoNorm[,i]
x <- seq(1, nrow(Sp.twoNorm.inv), 1)
y <- seq(1, ncol(Sp.twoNorm.inv), 1)
Sp.grid <- expand.grid(x=x, y=y)
Sp.grid$z <- as.vector(Sp.twoNorm.inv)

Sc.twoNorm.inv <- Sc.twoNorm
for(i in 1:ncol(Sc.twoNorm)) Sc.twoNorm.inv[,ncol(Sc.twoNorm)-i+1] <- Sc.twoNorm[,i]
x <- seq(1, nrow(Sc.twoNorm.inv), 1)
y <- seq(1, ncol(Sc.twoNorm.inv), 1)
Sc.grid <- expand.grid(x=x, y=y)
Sc.grid$z <- as.vector(Sc.twoNorm.inv)

Mm.twoNorm.inv <- Mm.twoNorm
for(i in 1:ncol(Mm.twoNorm)) Mm.twoNorm.inv[,ncol(Mm.twoNorm)-i+1] <- Mm.twoNorm[,i]
x <- seq(1, nrow(Mm.twoNorm.inv), 1)
y <- seq(1, ncol(Mm.twoNorm.inv), 1)
Mm.grid <- expand.grid(x=x, y=y)
Mm.grid$z <- as.vector(Mm.twoNorm.inv)

save(list = c("Sp.grid", "Sc.grid", "Mm.grid"), file = "dinucNorm2.RData")


