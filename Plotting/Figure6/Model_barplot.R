library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

Arch <- c("Accuracy" = 0.73949, "Precision" = 0.63316, "Recall" = 0.77066)

setwd("~/dm6_promoter_correction/Figure6_stuff")
# Plotting Arch
pdf(file = "Performance_Barplot.pdf",
width = 6, height = 6, pointsize = 14)
par(mar=c(4,3,4,5.5), cex = 1.2)
par(xpd = NA)
    barplot(Arch, col = cbbPalette[2:4], beside = F, ylim = c(0,1),
            main = "Common vs Putative\nConfusion Matrix Statistics")
dev.off()
