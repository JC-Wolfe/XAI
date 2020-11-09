library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)
library(gplots)
library(lattice)

load("/home/jw18713/Project1/data/model_predictions/S2_FuzzyGR.Rda") # gro
ann_genome <- get(load(file = "/home/jw18713/annotated_GRange_FI_done.Rda"))
S2_enhancers <- get(load("grown_predicted_enhancers_S2.Rda"))
S2_starr <- get(load("starr_S2.Rda"))

S2_starr_and_predicted <- reduce(c(S2_enhancers[S2_enhancers %over% S2_starr],
                                S2_starr[S2_starr %over% S2_enhancers]))
S2_fuzzy <- S2_enhancers[!S2_enhancers %over% S2_starr]
S2_starr_specific <- S2_starr[!S2_starr %over% S2_enhancers]


S2_both <- gro[gro %over% S2_starr_and_predicted]
S2_fuzzy_only <- gro[gro %over% S2_fuzzy]
S2_STARR_only <- gro[gro %over% S2_starr_specific]
S2_neither <- gro[!gro %over% S2_both & !gro %over% S2_fuzzy_only & !gro %over% S2_STARR_only]


gro$group <- 0
gro$group[which(gro %over% S2_starr_and_predicted)] <- "Both"
gro$group[which(gro %over% S2_fuzzy)] <- "XAI"
gro$group[which(gro %over% S2_starr_specific)] <- "STARR"
gro$group[which(!gro %over% S2_both & !gro %over% S2_fuzzy_only & !gro %over% S2_STARR_only)] <- "Neither"

gro$group <- as.factor(gro$group)

gro$annotations <- 0
overlaps <- findOverlaps(gro, ann_genome)
scores <- rep(0, length(gro))
scores[as.vector(queryHits(overlaps))] <- ann_genome$annotations[as.vector(subjectHits(overlaps))]
gro$annotations <- scores

results_matrix <- matrix(0, 5, 7)
colnames(results_matrix) <- c("Intergenic", "Promoter", "First Intron", "Other Introns", "Exon", "5' UTR", "3' UTR")
rownames(results_matrix) <- c("Both", "XAI Only", "STARR-seq Only", "Neither", "Whole Genome")


results_matrix[1,1] <- length(which(gro$group == "Both" & gro$annotations == "intergenic"))
results_matrix[1,2] <- length(which(gro$group == "Both" & gro$annotations == "promoter"))
results_matrix[1,3] <- length(which(gro$group == "Both" & gro$annotations == "first_intron"))
results_matrix[1,4] <- length(which(gro$group == "Both" & gro$annotations == "intron"))
results_matrix[1,5] <- length(which(gro$group == "Both" & gro$annotations == "exon"))
results_matrix[1,6] <- length(which(gro$group == "Both" & gro$annotations == "five_prime_UTR"))
results_matrix[1,7] <- length(which(gro$group == "Both" & gro$annotations == "three_prime_UTR"))

results_matrix[2,1] <- length(which(gro$group == "XAI" & gro$annotations == "intergenic"))
results_matrix[2,2] <- length(which(gro$group == "XAI" & gro$annotations == "promoter"))
results_matrix[2,3] <- length(which(gro$group == "XAI" & gro$annotations == "first_intron"))
results_matrix[2,4] <- length(which(gro$group == "XAI" & gro$annotations == "intron"))
results_matrix[2,5] <- length(which(gro$group == "XAI" & gro$annotations == "exon"))
results_matrix[2,6] <- length(which(gro$group == "XAI" & gro$annotations == "five_prime_UTR"))
results_matrix[2,7] <- length(which(gro$group == "XAI" & gro$annotations == "three_prime_UTR"))

results_matrix[3,1] <- length(which(gro$group == "STARR" & gro$annotations == "intergenic"))
results_matrix[3,2] <- length(which(gro$group == "STARR" & gro$annotations == "promoter"))
results_matrix[3,3] <- length(which(gro$group == "STARR" & gro$annotations == "first_intron"))
results_matrix[3,4] <- length(which(gro$group == "STARR" & gro$annotations == "intron"))
results_matrix[3,5] <- length(which(gro$group == "STARR" & gro$annotations == "exon"))
results_matrix[3,6] <- length(which(gro$group == "STARR" & gro$annotations == "five_prime_UTR"))
results_matrix[3,7] <- length(which(gro$group == "STARR" & gro$annotations == "three_prime_UTR"))

results_matrix[4,1] <- length(which(gro$group == "Neither" & gro$annotations == "intergenic"))
results_matrix[4,2] <- length(which(gro$group == "Neither" & gro$annotations == "promoter"))
results_matrix[4,3] <- length(which(gro$group == "Neither" & gro$annotations == "first_intron"))
results_matrix[4,4] <- length(which(gro$group == "Neither" & gro$annotations == "intron"))
results_matrix[4,5] <- length(which(gro$group == "Neither" & gro$annotations == "exon"))
results_matrix[4,6] <- length(which(gro$group == "Neither" & gro$annotations == "five_prime_UTR"))
results_matrix[4,7] <- length(which(gro$group == "Neither" & gro$annotations == "three_prime_UTR"))

results_matrix[5,1] <- length(which(ann_build$annotations == "intergenic"))
results_matrix[5,2] <- length(which(ann_build$annotations == "promoter"))
results_matrix[5,3] <- length(which(ann_build$annotations == "first_intron"))
results_matrix[5,4] <- length(which(ann_build$annotations == "intron"))
results_matrix[5,5] <- length(which(ann_build$annotations == "exon"))
results_matrix[5,6] <- length(which(ann_build$annotations == "five_prime_UTR"))
results_matrix[5,7] <- length(which(ann_build$annotations == "three_prime_UTR"))

hmap <- results_matrix
results_matrix <- t(results_matrix)
results_matrix_aves <- results_matrix

results_matrix_aves[,1] <- results_matrix[,1]/sum(results_matrix[,1])
results_matrix_aves[,2] <- results_matrix[,2]/sum(results_matrix[,2])
results_matrix_aves[,3] <- results_matrix[,3]/sum(results_matrix[,3])
results_matrix_aves[,4] <- results_matrix[,4]/sum(results_matrix[,4])
results_matrix_aves[,5] <- results_matrix[,5]/sum(results_matrix[,5])



cs <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pdf("/home/jw18713/Project1/Paper_Plots/Figure2/Figure2D.pdf")
par(mar=c(4,4,4,10.5))
par(xpd = NA)
barplot(results_matrix_aves, col=cs, xlab="Result Classification", ylab="Percentage of Results in Region", main="S2 Annotation Comparisons", yaxt = "none")
axis(2, seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10), "%"), las=2)
legend(6.25,0.75,rownames(results_matrix_aves), fill=cs, bty = "n")
dev.off()



hmap_aves <- hmap


hmap_aves[1,] <- hmap[1,]/sum(hmap[1,])
hmap_aves[2,] <- hmap[2,]/sum(hmap[2,])
hmap_aves[3,] <- hmap[3,]/sum(hmap[3,])
hmap_aves[4,] <- hmap[4,]/sum(hmap[4,])
hmap_aves[5,] <- hmap[5,]/sum(hmap[5,])

hmap_comp <- hmap_aves[1:4,]
hmap_comp[1,] <- log2(hmap_comp[1,]/hmap_aves[5,])
hmap_comp[2,] <- log2(hmap_comp[2,]/hmap_aves[5,])
hmap_comp[3,] <- log2(hmap_comp[3,]/hmap_aves[5,])
hmap_comp[4,] <- log2(hmap_comp[4,]/hmap_aves[5,])



cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
custom_at <- seq(-2.5, 2.5, by = (5/30))
cols_contrast <- rev(c(colorRampPalette(c(cbbPalette[7], "white", cbbPalette[6]))(30)))



pdf(file = "/home/jw18713/Project1/Paper_Plots/Figure2/Figure2D_Var.pdf")
  levelplot(t(hmap_comp),
    at = custom_at,
    main = "Predicted Region Annotation Comparisons S2",
    xlab = "log2 Enrichment Difference",
    ylab = "Predicted By",
    scales=list(x=list(rot=90)),
    col.regions = cols_contrast)
dev.off()
