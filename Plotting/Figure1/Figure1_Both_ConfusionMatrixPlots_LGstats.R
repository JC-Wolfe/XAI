library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# BG3 matrix
BG3_cmat <- matrix(0,2,3)
colnames(BG3_cmat) <- c("Accuracy", "Precision", "Recall")
rownames(BG3_cmat) <- c("XAI", "NN")

BG3_cmat[1,1] <- 0.83237
BG3_cmat[1,2] <- 0.52308
BG3_cmat[1,3] <- 0.76792

BG3_cmat[2,1] <- 0.85409
BG3_cmat[2,2] <- 0.52851
BG3_cmat[2,3] <- 0.79842

# S2 Matrix

S2_cmat <- matrix(0,2,3)
colnames(S2_cmat) <- c("Accuracy", "Precision", "Recall")
rownames(S2_cmat) <- c("XAI", "NN")

S2_cmat[1,1] <- 0.84778
S2_cmat[1,2] <- 0.52101
S2_cmat[1,3] <- 0.74141

S2_cmat[2,1] <- 0.87596
S2_cmat[2,2] <- 0.52618
S2_cmat[2,3] <- 0.75430

# Plotting BG3
pdf(file = "/home/jw18713/Project1/Paper_Plots/Figure1/Figure1B(LG).pdf",
width = 12, height = 6, pointsize = 14)
par(mar=c(4,3,4,5.5), mfrow=c(1,2), cex = 1.2)
par(xpd = NA)
    barplot(BG3_cmat, col = cbbPalette[2:3], beside = T, ylim = c(0,1),
            main = "BG3 Confusion Matrix Statistics"
            )
    barplot(S2_cmat, col = cbbPalette[2:3], beside = T, ylim = c(0,1),
            main = "S2 Confusion Matrix Statistics"
            )
    legend(10,0.75,rownames(BG3_cmat), fill=cbbPalette[2:3])
dev.off()
