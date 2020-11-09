library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)
library(gplots)
library(lattice)

states <- import("/home/jw18713/Quicktest/BG3Kc-state11-BG3.bed")

BG3_gro <- get(load("/home/jw18713/Project1/data/model_predictions/BG3_FuzzyGR.Rda")) #Imported as BG3_gro
BG3_enhancers <- get(load("~/Archive_Year1/grown_predicted_enhancers_BG3.Rda"))
BG3_starr <- get(load("~/Archive_Year1/starr.Rda"))

BG3_starr_and_predicted <- reduce(c(BG3_enhancers[BG3_enhancers %over% BG3_starr],
                                BG3_starr[BG3_starr %over% BG3_enhancers]))
BG3_fuzzy <- BG3_enhancers[!BG3_enhancers %over% BG3_starr]
BG3_starr_specific <- BG3_starr[!BG3_starr %over% BG3_enhancers]


BG3_both <- gro[gro %over% BG3_starr_and_predicted]
BG3_fuzzy_only <- gro[gro %over% BG3_fuzzy]
BG3_STARR_only <- gro[gro %over% BG3_starr_specific]
BG3_neither <- gro[!gro %over% BG3_both & !gro %over% BG3_fuzzy_only & !gro %over% BG3_STARR_only]

BG3_gro$group <- 0
BG3_gro$group[which(BG3_gro %over% BG3_starr_and_predicted)] <- "Both"
BG3_gro$group[which(BG3_gro %over% BG3_fuzzy)] <- "XAI"
BG3_gro$group[which(BG3_gro %over% BG3_starr_specific)] <- "STARR"
BG3_gro$group[which(!BG3_gro %over% BG3_both & !BG3_gro %over% BG3_fuzzy_only & !BG3_gro %over% BG3_STARR_only)] <- "Neither"

BG3_gro$group <- as.factor(BG3_gro$group)

BG3_overlaps <- findOverlaps(BG3_gro, states)
BG3_covered <- BG3_gro[queryHits(BG3_overlaps)]

BG3_covered$state <- as.factor(states$name[subjectHits(BG3_overlaps)])
BG3_covered$rgb <- states$itemRgb[subjectHits(BG3_overlaps)]

BG3_results_matrix <- matrix(0, 5, 11)

# BG3 cells
BG3_results_matrix[1,1] <- length(which(BG3_covered$group == "Both" & BG3_covered$state == 1))/length(which(BG3_covered$group == "Both"))
BG3_results_matrix[1,2] <- length(which(BG3_covered$group == "Both" & BG3_covered$state == 2))/length(which(BG3_covered$group == "Both"))
BG3_results_matrix[1,3] <- length(which(BG3_covered$group == "Both" & BG3_covered$state == 3))/length(which(BG3_covered$group == "Both"))
BG3_results_matrix[1,4] <- length(which(BG3_covered$group == "Both" & BG3_covered$state == 4))/length(which(BG3_covered$group == "Both"))
BG3_results_matrix[1,5] <- length(which(BG3_covered$group == "Both" & BG3_covered$state == 5))/length(which(BG3_covered$group == "Both"))
BG3_results_matrix[1,6] <- length(which(BG3_covered$group == "Both" & BG3_covered$state == 6))/length(which(BG3_covered$group == "Both"))
BG3_results_matrix[1,7] <- length(which(BG3_covered$group == "Both" & BG3_covered$state == 7))/length(which(BG3_covered$group == "Both"))
BG3_results_matrix[1,8] <- length(which(BG3_covered$group == "Both" & BG3_covered$state == 8))/length(which(BG3_covered$group == "Both"))
BG3_results_matrix[1,9] <- length(which(BG3_covered$group == "Both" & BG3_covered$state == 9))/length(which(BG3_covered$group == "Both"))
BG3_results_matrix[1,10] <- length(which(BG3_covered$group == "Both" & BG3_covered$state == 10))/length(which(BG3_covered$group == "Both"))
BG3_results_matrix[1,11] <- length(which(BG3_covered$group == "Both" & BG3_covered$state == 11))/length(which(BG3_covered$group == "Both"))

BG3_results_matrix[2,1] <- length(which(BG3_covered$group == "XAI" & BG3_covered$state == 1))/length(which(BG3_covered$group == "XAI"))
BG3_results_matrix[2,2] <- length(which(BG3_covered$group == "XAI" & BG3_covered$state == 2))/length(which(BG3_covered$group == "XAI"))
BG3_results_matrix[2,3] <- length(which(BG3_covered$group == "XAI" & BG3_covered$state == 3))/length(which(BG3_covered$group == "XAI"))
BG3_results_matrix[2,4] <- length(which(BG3_covered$group == "XAI" & BG3_covered$state == 4))/length(which(BG3_covered$group == "XAI"))
BG3_results_matrix[2,5] <- length(which(BG3_covered$group == "XAI" & BG3_covered$state == 5))/length(which(BG3_covered$group == "XAI"))
BG3_results_matrix[2,6] <- length(which(BG3_covered$group == "XAI" & BG3_covered$state == 6))/length(which(BG3_covered$group == "XAI"))
BG3_results_matrix[2,7] <- length(which(BG3_covered$group == "XAI" & BG3_covered$state == 7))/length(which(BG3_covered$group == "XAI"))
BG3_results_matrix[2,8] <- length(which(BG3_covered$group == "XAI" & BG3_covered$state == 8))/length(which(BG3_covered$group == "XAI"))
BG3_results_matrix[2,9] <- length(which(BG3_covered$group == "XAI" & BG3_covered$state == 9))/length(which(BG3_covered$group == "XAI"))
BG3_results_matrix[2,10] <- length(which(BG3_covered$group == "XAI" & BG3_covered$state == 10))/length(which(BG3_covered$group == "XAI"))
BG3_results_matrix[2,11] <- length(which(BG3_covered$group == "XAI" & BG3_covered$state == 11))/length(which(BG3_covered$group == "XAI"))

BG3_results_matrix[3,1] <- length(which(BG3_covered$group == "STARR" & BG3_covered$state == 1))/length(which(BG3_covered$group == "STARR"))
BG3_results_matrix[3,2] <- length(which(BG3_covered$group == "STARR" & BG3_covered$state == 2))/length(which(BG3_covered$group == "STARR"))
BG3_results_matrix[3,3] <- length(which(BG3_covered$group == "STARR" & BG3_covered$state == 3))/length(which(BG3_covered$group == "STARR"))
BG3_results_matrix[3,4] <- length(which(BG3_covered$group == "STARR" & BG3_covered$state == 4))/length(which(BG3_covered$group == "STARR"))
BG3_results_matrix[3,5] <- length(which(BG3_covered$group == "STARR" & BG3_covered$state == 5))/length(which(BG3_covered$group == "STARR"))
BG3_results_matrix[3,6] <- length(which(BG3_covered$group == "STARR" & BG3_covered$state == 6))/length(which(BG3_covered$group == "STARR"))
BG3_results_matrix[3,7] <- length(which(BG3_covered$group == "STARR" & BG3_covered$state == 7))/length(which(BG3_covered$group == "STARR"))
BG3_results_matrix[3,8] <- length(which(BG3_covered$group == "STARR" & BG3_covered$state == 8))/length(which(BG3_covered$group == "STARR"))
BG3_results_matrix[3,9] <- length(which(BG3_covered$group == "STARR" & BG3_covered$state == 9))/length(which(BG3_covered$group == "STARR"))
BG3_results_matrix[3,10] <- length(which(BG3_covered$group == "STARR" & BG3_covered$state == 10))/length(which(BG3_covered$group == "STARR"))
BG3_results_matrix[3,11] <- length(which(BG3_covered$group == "STARR" & BG3_covered$state == 11))/length(which(BG3_covered$group == "STARR"))

BG3_results_matrix[4,1] <- length(which(BG3_covered$group == "Neither" & BG3_covered$state == 1))/length(which(BG3_covered$group == "Neither"))
BG3_results_matrix[4,2] <- length(which(BG3_covered$group == "Neither" & BG3_covered$state == 2))/length(which(BG3_covered$group == "Neither"))
BG3_results_matrix[4,3] <- length(which(BG3_covered$group == "Neither" & BG3_covered$state == 3))/length(which(BG3_covered$group == "Neither"))
BG3_results_matrix[4,4] <- length(which(BG3_covered$group == "Neither" & BG3_covered$state == 4))/length(which(BG3_covered$group == "Neither"))
BG3_results_matrix[4,5] <- length(which(BG3_covered$group == "Neither" & BG3_covered$state == 5))/length(which(BG3_covered$group == "Neither"))
BG3_results_matrix[4,6] <- length(which(BG3_covered$group == "Neither" & BG3_covered$state == 6))/length(which(BG3_covered$group == "Neither"))
BG3_results_matrix[4,7] <- length(which(BG3_covered$group == "Neither" & BG3_covered$state == 7))/length(which(BG3_covered$group == "Neither"))
BG3_results_matrix[4,8] <- length(which(BG3_covered$group == "Neither" & BG3_covered$state == 8))/length(which(BG3_covered$group == "Neither"))
BG3_results_matrix[4,9] <- length(which(BG3_covered$group == "Neither" & BG3_covered$state == 9))/length(which(BG3_covered$group == "Neither"))
BG3_results_matrix[4,10] <- length(which(BG3_covered$group == "Neither" & BG3_covered$state == 10))/length(which(BG3_covered$group == "Neither"))
BG3_results_matrix[4,11] <- length(which(BG3_covered$group == "Neither" & BG3_covered$state == 11))/length(which(BG3_covered$group == "Neither"))

BG3_results_matrix[5,1] <- length(which(BG3_covered$state == 1))/length(BG3_covered)
BG3_results_matrix[5,2] <- length(which(BG3_covered$state == 2))/length(BG3_covered)
BG3_results_matrix[5,3] <- length(which(BG3_covered$state == 3))/length(BG3_covered)
BG3_results_matrix[5,4] <- length(which(BG3_covered$state == 4))/length(BG3_covered)
BG3_results_matrix[5,5] <- length(which(BG3_covered$state == 5))/length(BG3_covered)
BG3_results_matrix[5,6] <- length(which(BG3_covered$state == 6))/length(BG3_covered)
BG3_results_matrix[5,7] <- length(which(BG3_covered$state == 7))/length(BG3_covered)
BG3_results_matrix[5,8] <- length(which(BG3_covered$state == 8))/length(BG3_covered)
BG3_results_matrix[5,9] <- length(which(BG3_covered$state == 9))/length(BG3_covered)
BG3_results_matrix[5,10] <- length(which(BG3_covered$state == 10))/length(BG3_covered)
BG3_results_matrix[5,11] <- length(which(BG3_covered$state == 11))/length(BG3_covered)



S2_gro <- get(load("/home/jw18713/Project1/data/model_predictions/S2_FuzzyGR.Rda")) #Imported as S2_gro
S2_enhancers <- get(load("~/Archive_Year1/grown_predicted_enhancers_S2.Rda"))
S2_starr <- get(load("~/Archive_Year1/starr_S2.Rda"))


S2_starr_and_predicted <- reduce(c(S2_enhancers[S2_enhancers %over% S2_starr],
                                S2_starr[S2_starr %over% S2_enhancers]))
S2_fuzzy <- S2_enhancers[!S2_enhancers %over% S2_starr]
S2_starr_specific <- S2_starr[!S2_starr %over% S2_enhancers]


S2_both <- gro[gro %over% S2_starr_and_predicted]
S2_fuzzy_only <- gro[gro %over% S2_fuzzy]
S2_STARR_only <- gro[gro %over% S2_starr_specific]
S2_neither <- gro[!gro %over% S2_both & !gro %over% S2_fuzzy_only & !gro %over% S2_STARR_only]

S2_gro$group <- 0
S2_gro$group[which(S2_gro %over% S2_starr_and_predicted)] <- "Both"
S2_gro$group[which(S2_gro %over% S2_fuzzy)] <- "XAI"
S2_gro$group[which(S2_gro %over% S2_starr_specific)] <- "STARR"
S2_gro$group[which(!S2_gro %over% S2_both & !S2_gro %over% S2_fuzzy_only & !S2_gro %over% S2_STARR_only)] <- "Neither"

S2_gro$group <- as.factor(S2_gro$group)

S2_overlaps <- findOverlaps(S2_gro, states)
S2_covered <- S2_gro[queryHits(S2_overlaps)]

S2_covered$state <- as.factor(states$name[subjectHits(S2_overlaps)])
S2_covered$rgb <- states$itemRgb[subjectHits(S2_overlaps)]

S2_results_matrix <- matrix(0, 5, 11)

# S2 cells
S2_results_matrix[1,1] <- length(which(S2_covered$group == "Both" & S2_covered$state == 1))/length(which(S2_covered$group == "Both"))
S2_results_matrix[1,2] <- length(which(S2_covered$group == "Both" & S2_covered$state == 2))/length(which(S2_covered$group == "Both"))
S2_results_matrix[1,3] <- length(which(S2_covered$group == "Both" & S2_covered$state == 3))/length(which(S2_covered$group == "Both"))
S2_results_matrix[1,4] <- length(which(S2_covered$group == "Both" & S2_covered$state == 4))/length(which(S2_covered$group == "Both"))
S2_results_matrix[1,5] <- length(which(S2_covered$group == "Both" & S2_covered$state == 5))/length(which(S2_covered$group == "Both"))
S2_results_matrix[1,6] <- length(which(S2_covered$group == "Both" & S2_covered$state == 6))/length(which(S2_covered$group == "Both"))
S2_results_matrix[1,7] <- length(which(S2_covered$group == "Both" & S2_covered$state == 7))/length(which(S2_covered$group == "Both"))
S2_results_matrix[1,8] <- length(which(S2_covered$group == "Both" & S2_covered$state == 8))/length(which(S2_covered$group == "Both"))
S2_results_matrix[1,9] <- length(which(S2_covered$group == "Both" & S2_covered$state == 9))/length(which(S2_covered$group == "Both"))
S2_results_matrix[1,10] <- length(which(S2_covered$group == "Both" & S2_covered$state == 10))/length(which(S2_covered$group == "Both"))
S2_results_matrix[1,11] <- length(which(S2_covered$group == "Both" & S2_covered$state == 11))/length(which(S2_covered$group == "Both"))

S2_results_matrix[2,1] <- length(which(S2_covered$group == "XAI" & S2_covered$state == 1))/length(which(S2_covered$group == "XAI"))
S2_results_matrix[2,2] <- length(which(S2_covered$group == "XAI" & S2_covered$state == 2))/length(which(S2_covered$group == "XAI"))
S2_results_matrix[2,3] <- length(which(S2_covered$group == "XAI" & S2_covered$state == 3))/length(which(S2_covered$group == "XAI"))
S2_results_matrix[2,4] <- length(which(S2_covered$group == "XAI" & S2_covered$state == 4))/length(which(S2_covered$group == "XAI"))
S2_results_matrix[2,5] <- length(which(S2_covered$group == "XAI" & S2_covered$state == 5))/length(which(S2_covered$group == "XAI"))
S2_results_matrix[2,6] <- length(which(S2_covered$group == "XAI" & S2_covered$state == 6))/length(which(S2_covered$group == "XAI"))
S2_results_matrix[2,7] <- length(which(S2_covered$group == "XAI" & S2_covered$state == 7))/length(which(S2_covered$group == "XAI"))
S2_results_matrix[2,8] <- length(which(S2_covered$group == "XAI" & S2_covered$state == 8))/length(which(S2_covered$group == "XAI"))
S2_results_matrix[2,9] <- length(which(S2_covered$group == "XAI" & S2_covered$state == 9))/length(which(S2_covered$group == "XAI"))
S2_results_matrix[2,10] <- length(which(S2_covered$group == "XAI" & S2_covered$state == 10))/length(which(S2_covered$group == "XAI"))
S2_results_matrix[2,11] <- length(which(S2_covered$group == "XAI" & S2_covered$state == 11))/length(which(S2_covered$group == "XAI"))

S2_results_matrix[3,1] <- length(which(S2_covered$group == "STARR" & S2_covered$state == 1))/length(which(S2_covered$group == "STARR"))
S2_results_matrix[3,2] <- length(which(S2_covered$group == "STARR" & S2_covered$state == 2))/length(which(S2_covered$group == "STARR"))
S2_results_matrix[3,3] <- length(which(S2_covered$group == "STARR" & S2_covered$state == 3))/length(which(S2_covered$group == "STARR"))
S2_results_matrix[3,4] <- length(which(S2_covered$group == "STARR" & S2_covered$state == 4))/length(which(S2_covered$group == "STARR"))
S2_results_matrix[3,5] <- length(which(S2_covered$group == "STARR" & S2_covered$state == 5))/length(which(S2_covered$group == "STARR"))
S2_results_matrix[3,6] <- length(which(S2_covered$group == "STARR" & S2_covered$state == 6))/length(which(S2_covered$group == "STARR"))
S2_results_matrix[3,7] <- length(which(S2_covered$group == "STARR" & S2_covered$state == 7))/length(which(S2_covered$group == "STARR"))
S2_results_matrix[3,8] <- length(which(S2_covered$group == "STARR" & S2_covered$state == 8))/length(which(S2_covered$group == "STARR"))
S2_results_matrix[3,9] <- length(which(S2_covered$group == "STARR" & S2_covered$state == 9))/length(which(S2_covered$group == "STARR"))
S2_results_matrix[3,10] <- length(which(S2_covered$group == "STARR" & S2_covered$state == 10))/length(which(S2_covered$group == "STARR"))
S2_results_matrix[3,11] <- length(which(S2_covered$group == "STARR" & S2_covered$state == 11))/length(which(S2_covered$group == "STARR"))

S2_results_matrix[4,1] <- length(which(S2_covered$group == "Neither" & S2_covered$state == 1))/length(which(S2_covered$group == "Neither"))
S2_results_matrix[4,2] <- length(which(S2_covered$group == "Neither" & S2_covered$state == 2))/length(which(S2_covered$group == "Neither"))
S2_results_matrix[4,3] <- length(which(S2_covered$group == "Neither" & S2_covered$state == 3))/length(which(S2_covered$group == "Neither"))
S2_results_matrix[4,4] <- length(which(S2_covered$group == "Neither" & S2_covered$state == 4))/length(which(S2_covered$group == "Neither"))
S2_results_matrix[4,5] <- length(which(S2_covered$group == "Neither" & S2_covered$state == 5))/length(which(S2_covered$group == "Neither"))
S2_results_matrix[4,6] <- length(which(S2_covered$group == "Neither" & S2_covered$state == 6))/length(which(S2_covered$group == "Neither"))
S2_results_matrix[4,7] <- length(which(S2_covered$group == "Neither" & S2_covered$state == 7))/length(which(S2_covered$group == "Neither"))
S2_results_matrix[4,8] <- length(which(S2_covered$group == "Neither" & S2_covered$state == 8))/length(which(S2_covered$group == "Neither"))
S2_results_matrix[4,9] <- length(which(S2_covered$group == "Neither" & S2_covered$state == 9))/length(which(S2_covered$group == "Neither"))
S2_results_matrix[4,10] <- length(which(S2_covered$group == "Neither" & S2_covered$state == 10))/length(which(S2_covered$group == "Neither"))
S2_results_matrix[4,11] <- length(which(S2_covered$group == "Neither" & S2_covered$state == 11))/length(which(S2_covered$group == "Neither"))

S2_results_matrix[5,1] <- length(which(S2_covered$state == 1))/length(S2_covered)
S2_results_matrix[5,2] <- length(which(S2_covered$state == 2))/length(S2_covered)
S2_results_matrix[5,3] <- length(which(S2_covered$state == 3))/length(S2_covered)
S2_results_matrix[5,4] <- length(which(S2_covered$state == 4))/length(S2_covered)
S2_results_matrix[5,5] <- length(which(S2_covered$state == 5))/length(S2_covered)
S2_results_matrix[5,6] <- length(which(S2_covered$state == 6))/length(S2_covered)
S2_results_matrix[5,7] <- length(which(S2_covered$state == 7))/length(S2_covered)
S2_results_matrix[5,8] <- length(which(S2_covered$state == 8))/length(S2_covered)
S2_results_matrix[5,9] <- length(which(S2_covered$state == 9))/length(S2_covered)
S2_results_matrix[5,10] <- length(which(S2_covered$state == 10))/length(S2_covered)
S2_results_matrix[5,11] <- length(which(S2_covered$state == 11))/length(S2_covered)






names<- c("TSS proximal",
                "Elongation",
                "Enhancer",
                "Active TSS",
                "Active intron",
                "Weak intron",
                "Competent",
                "Polycomb",
                "Heterochromatin",
                "Heterochr. in eu.",
                "Basal")

cs <- rep(0, 11)
for (i in seq_along(cs)){
  colpicker <- unique(BG3_covered$rgb[BG3_covered$state == i])
  cs[i] <- colpicker
}

BG3_plot_matrix <- t(BG3_results_matrix)
colnames(BG3_plot_matrix) <- c("Both", "XAI Only", "STARR-seq\nOnly", "Neither", "Whole Genome")

S2_plot_matrix <- t(S2_results_matrix)
colnames(S2_plot_matrix) <- c("Both", "XAI Only", "STARR-seq\nOnly", "Neither", "Whole Genome")

pdf("/home/jw18713/Project1/Paper_Plots/Lenka_compare.pdf",
  width = 12, height = 16, pointsize = 14)
par(mar=c(4,4,4,12.5), mfrow = c(2,1), cex = 1.2)
par(xpd = NA)
barplot(BG3_plot_matrix, col=cs, xlab="Predicted By", ylab="Percentage of BG3_results in Region", main="Chromatin States of\n Predicted Enhancers in BG3 Cells", yaxt = "none")
axis(2, seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10), "%"), las=2)
legend(6.5,0.75, names, fill=cs, bty = "n")

barplot(S2_plot_matrix, col=cs, xlab="Predicted By", ylab="Percentage of BG3_results in Region", main="Chromatin States of\nPredicted Enhancers in S2 Cells", yaxt = "none")
axis(2, seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10), "%"), las=2)
legend(6.5,0.75, names, fill=cs, bty = "n")

dev.off()

BG3_hmap <- BG3_results_matrix

BG3_hmap_aves <- BG3_hmap

BG3_hmap_aves[1,] <- BG3_hmap[1,]/sum(BG3_hmap[1,])
BG3_hmap_aves[2,] <- BG3_hmap[2,]/sum(BG3_hmap[2,])
BG3_hmap_aves[3,] <- BG3_hmap[3,]/sum(BG3_hmap[3,])
BG3_hmap_aves[4,] <- BG3_hmap[4,]/sum(BG3_hmap[4,])
BG3_hmap_aves[5,] <- BG3_hmap[5,]/sum(BG3_hmap[5,])

BG3_hmap_comp <- BG3_hmap_aves[1:4,]
BG3_hmap_comp[1,] <- log2(BG3_hmap_comp[1,]/BG3_hmap_aves[5,])
BG3_hmap_comp[2,] <- log2(BG3_hmap_comp[2,]/BG3_hmap_aves[5,])
BG3_hmap_comp[3,] <- log2(BG3_hmap_comp[3,]/BG3_hmap_aves[5,])
BG3_hmap_comp[4,] <- log2(BG3_hmap_comp[4,]/BG3_hmap_aves[5,])


S2_hmap <- S2_results_matrix

S2_hmap_aves <- S2_hmap

S2_hmap_aves[1,] <- S2_hmap[1,]/sum(S2_hmap[1,])
S2_hmap_aves[2,] <- S2_hmap[2,]/sum(S2_hmap[2,])
S2_hmap_aves[3,] <- S2_hmap[3,]/sum(S2_hmap[3,])
S2_hmap_aves[4,] <- S2_hmap[4,]/sum(S2_hmap[4,])
S2_hmap_aves[5,] <- S2_hmap[5,]/sum(S2_hmap[5,])

S2_hmap_comp <- S2_hmap_aves[1:4,]
S2_hmap_comp[1,] <- log2(S2_hmap_comp[1,]/S2_hmap_aves[5,])
S2_hmap_comp[2,] <- log2(S2_hmap_comp[2,]/S2_hmap_aves[5,])
S2_hmap_comp[3,] <- log2(S2_hmap_comp[3,]/S2_hmap_aves[5,])
S2_hmap_comp[4,] <- log2(S2_hmap_comp[4,]/S2_hmap_aves[5,])

colnames(S2_hmap_comp) <- names
colnames(BG3_hmap_comp) <- names
rownames(S2_hmap_comp) <- c("Both", "XAI Only", "STARR-seq\nOnly", "Neither")
rownames(BG3_hmap_comp) <- c("Both", "XAI Only", "STARR-seq\nOnly", "Neither")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
custom_at <- seq(-5.5, 5.5, by = (11/30))
cols_contrast <- rev(c(colorRampPalette(c(cbbPalette[7], "white", cbbPalette[6]))(30)))



pdf(file = "/home/jw18713/Project1/Paper_Plots/Lenka_compare_hmap_BG3.pdf")

  levelplot(t(BG3_hmap_comp),
    at = custom_at,
    main = "Predicted Region Annotation Comparisons BG3",
    xlab = "log2 Enrichment Difference",
    ylab = "Predicted By",
    scales=list(x=list(rot=90)),
    col.regions = cols_contrast)

dev.off()

pdf(file = "/home/jw18713/Project1/Paper_Plots/Lenka_compare_hmap_S2.pdf")

    levelplot(t(S2_hmap_comp),
      at = custom_at,
      main = "Predicted Region Annotation Comparisons S2",
      xlab = "log2 Enrichment Difference",
      ylab = "Predicted By",
      scales=list(x=list(rot=90)),
      col.regions = cols_contrast)
dev.off()
