library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


S2_cmat <- matrix(0,1,3)
colnames(S2_cmat) <- c("Accuracy", "Precision", "Recall")
rownames(S2_cmat) <- c("XAI")

S2_cmat[1,1] <- 0.53812
S2_cmat[1,2] <- 0.55283
S2_cmat[1,3] <- 0.55668

# Plotting BG3
pdf(file = "/home/jw18713/Project1/Paper_Plots/Figure1/Dev_Hk_Bar.pdf",
width = 5.5, height = 6, pointsize = 14)
par(mar=c(4,3,4,4.5), cex = 1.2)
par(xpd = NA)
    barplot(S2_cmat, col = cbbPalette[2:4], beside = T, ylim = c(0,1),
            main = "Developmental vs Housekeeping\nConfusion Matrix Statistics"
            )
    legend(10,0.75,rownames(BG3_cmat), fill=cbbPalette[2:3])
dev.off()
