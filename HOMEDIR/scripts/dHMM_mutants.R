
library(Biostrings)
library(nuCpos)

setwd("RData")
load(file = "sc_genome_2011.RData")
sc.genome <- sc.genome.2011

seqBAR1 <- sc.genome[["chr9"]][(322342-262-3000):(322342-262+3000)]
seqBAR1[(3000+262-157):(3000+262-157+4)] <- DNAString("GTACC")
annotation = data.frame(name = c("alpha2", "BAR1", "TATA"), color = c("purple", "orange", "blue"), 
			left = c(3001, 3263, 3129), right = c(3026, 5026, 3135), 
			stringsAsFactors = FALSE)

seqA20 <- "AAAAAAAAAAAAAAAAAAAAGTAC"
seqA30 <- "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTAC"
seqCG4 <- "CGCGCGCGGTAC"
seqCG5 <- "CGCGCGCGCGGTAC"
seqCTG12 <- "TGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGTAC"
seqSac5 <- "ACGAGCTCGAGCTCGAGCTCGAGCTCGAGCTCGAGTAC"
seqSac6 <- "ACGAGCTCGAGCTCGAGCTCGAGCTCGAGCTCGAGCTCGAGTAC"

ylim.HBA <- c(-25, 15)
plot.window <- 501

wt <- mutNuCpos(seqBAR1, site = 3109, ins= "", species="sc", prob.dyad = TRUE, 
		smoothHBA = FALSE, plot.window = plot.window, ylim.HBA = ylim.HBA, 
		annotation = annotation, 
		full = TRUE)

A20 <- mutNuCpos(seqBAR1, site = 3109, ins= seqA20, species="sc", prob.dyad = TRUE, 
		smoothHBA = FALSE, plot.window = plot.window, ylim.HBA = ylim.HBA, 
		annotation = annotation, 
		full = TRUE)

A30 <- mutNuCpos(seqBAR1, site = 3109, ins= seqA30, species="sc", prob.dyad = TRUE, 
		smoothHBA = FALSE, plot.window = plot.window, ylim.HBA = ylim.HBA, 
		annotation = annotation, 
		full = TRUE)

CG4 <- mutNuCpos(seqBAR1, site = 3109, ins= seqCG4, species="sc", prob.dyad = TRUE, 
		smoothHBA = FALSE, plot.window = plot.window, ylim.HBA = ylim.HBA, 
		annotation = annotation, 
		full = TRUE)

CG5 <- mutNuCpos(seqBAR1, site = 3109, ins= seqCG5, species="sc", prob.dyad = TRUE, 
		smoothHBA = FALSE, plot.window = plot.window, ylim.HBA = ylim.HBA, 
		annotation = annotation, 
		full = TRUE)

CTG12 <- mutNuCpos(seqBAR1, site = 3109, ins= seqCTG12, species="sc", prob.dyad = TRUE, 
		smoothHBA = FALSE, plot.window = plot.window, ylim.HBA = ylim.HBA, 
		annotation = annotation, 
		full = TRUE)

Sac5 <- mutNuCpos(seqBAR1, site = 3109, ins= seqSac5, species="sc", prob.dyad = TRUE, 
		smoothHBA = FALSE, plot.window = plot.window, ylim.HBA = ylim.HBA, 
		annotation = annotation, 
		full = TRUE)

Sac6 <- mutNuCpos(seqBAR1, site = 3109, ins= seqSac6, species="sc", prob.dyad = TRUE, 
		smoothHBA = FALSE, plot.window = plot.window, ylim.HBA = ylim.HBA, 
		annotation = annotation, 
		full = TRUE)

objects <- c("wt", "A20", "A30", "CG4", "CG5", "CTG12", "Sac5", "Sac6", "annotation")
save(list = objects, file = "mutNuCpos_BAR1.RData")


TALS <- paste(scan(file =
	system.file("extdata", "TALS.fasta", package = "nuCpos"),
	what = character(), skip = 1), sep = "", collapse = "")

TTAGGGx2 <- paste(scan(file =
	system.file("extdata", "TTAGGGx2.fasta", package = "nuCpos"),
	what = character(), skip = 1), sep = "", collapse = "")

TTAGGGx4 <- paste(scan(file =
	system.file("extdata", "TTAGGGx4.fasta", package = "nuCpos"),
	what = character(), skip = 1), sep = "", collapse = "")

TTAGGGx6 <- paste(scan(file =
	system.file("extdata", "TTAGGGx6.fasta", package = "nuCpos"),
	what = character(), skip = 1), sep = "", collapse = "")

TTAGGGx12 <- paste(scan(file =
	system.file("extdata", "TTAGGGx12.fasta", package = "nuCpos"),
	what = character(), skip = 1), sep = "", collapse = "")

TTAGGGx29 <- paste(scan(file =
	system.file("extdata", "TTAGGGx29.fasta", package = "nuCpos"),
	what = character(), skip = 1), sep = "", collapse = "")

TGTAGGx6 <- paste(scan(file =
	system.file("extdata", "TGTAGGx6.fasta", package = "nuCpos"),
	what = character(), skip = 1), sep = "", collapse = "")

TGTAGGx12 <- paste(scan(file =
	system.file("extdata", "TGTAGGx12.fasta", package = "nuCpos"),
	what = character(), skip = 1), sep = "", collapse = "")

TGTGAGx6 <- paste(scan(file =
	system.file("extdata", "TGTGAGx6.fasta", package = "nuCpos"),
	what = character(), skip = 1), sep = "", collapse = "")

TGTGAGx12 <- paste(scan(file =
	system.file("extdata", "TGTGAGx12.fasta", package = "nuCpos"),
	what = character(), skip = 1), sep = "", collapse = "")

SITE <- 1464
annotation = data.frame(name = "alpha2", color = "purple",
	left = 1534, right = 1559, stringsAsFactors = FALSE)

ylim.HBA <- c(-20, 10)
plot.window <- 501

wt <- mutNuCpos(TALS, site = SITE, ins= "", species="sc", prob.dyad = TRUE, 
		smoothHBA = FALSE, plot.window = plot.window, ylim.HBA = ylim.HBA, 
		annotation = annotation, full = TRUE)

rTTAGGGx29 <- mutNuCpos(TALS, site = SITE, ins= TTAGGGx29, species="sc", prob.dyad = TRUE, 
		smoothHBA = FALSE, plot.window = plot.window, ylim.HBA = ylim.HBA, 
		annotation = annotation, full = TRUE)

rTTAGGGx2 <- mutNuCpos(TALS, site = SITE, ins= TTAGGGx2, species="sc", prob.dyad = TRUE, 
		smoothHBA = FALSE, plot.window = plot.window, ylim.HBA = ylim.HBA, 
		annotation = annotation, full = TRUE)

rTTAGGGx4 <- mutNuCpos(TALS, site = SITE, ins= TTAGGGx4, species="sc", prob.dyad = TRUE, 
		smoothHBA = FALSE, plot.window = plot.window, ylim.HBA = ylim.HBA, 
		annotation = annotation, full = TRUE)

rTTAGGGx6 <- mutNuCpos(TALS, site = SITE, ins= TTAGGGx6, species="sc", prob.dyad = TRUE, 
		smoothHBA = FALSE, plot.window = plot.window, ylim.HBA = ylim.HBA, 
		annotation = annotation, full = TRUE)

rTTAGGGx12 <- mutNuCpos(TALS, site = SITE, ins= TTAGGGx12, species="sc", prob.dyad = TRUE, 
		smoothHBA = FALSE, plot.window = plot.window, ylim.HBA = ylim.HBA, 
		annotation = annotation, full = TRUE)

rTGTAGGx6 <- mutNuCpos(TALS, site = SITE, ins= TGTAGGx6, species="sc", prob.dyad = TRUE, 
		smoothHBA = FALSE, plot.window = plot.window, ylim.HBA = ylim.HBA, 
		annotation = annotation, full = TRUE)

rTGTAGGx12 <- mutNuCpos(TALS, site = SITE, ins= TGTAGGx12, species="sc", prob.dyad = TRUE, 
		smoothHBA = FALSE, plot.window = plot.window, ylim.HBA = ylim.HBA, 
		annotation = annotation, full = TRUE)

rTGTGAGx6 <- mutNuCpos(TALS, site = SITE, ins= TGTGAGx6, species="sc", prob.dyad = TRUE, 
		smoothHBA = FALSE, plot.window = plot.window, ylim.HBA = ylim.HBA, 
		annotation = annotation, full = TRUE)

rTGTGAGx12 <- mutNuCpos(TALS, site = SITE, ins= TGTGAGx12, species="sc", prob.dyad = TRUE, 
		smoothHBA = FALSE, plot.window = plot.window, ylim.HBA = ylim.HBA, 
		annotation = annotation, full = TRUE)


objects <- c("wt", "rTTAGGGx29", "rTTAGGGx2", "rTTAGGGx4", 
	"rTTAGGGx6", "rTTAGGGx12", "rTGTAGGx6", "rTGTAGGx12", 
	"rTGTGAGx6", "rTGTGAGx12", "annotation", "SITE")
save(list = objects, file = "mutNuCpos_TALS.RData")
