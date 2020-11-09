library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)
library(VennDiagram)

S2_enhancers <- get(load("grown_predicted_enhancers_S2.Rda"))
BG3_enhancers <- get(load("grown_predicted_enhancers_BG3.Rda"))
S2_starr <- get(load("starr_S2.Rda"))
BG3_starr <- get(load("starr.Rda"))
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pdf("/home/jw18713/Project1/Paper_Plots/Figure2/BG3_Enhancer_Widths.pdf")
hist(log2(width(BG3_enhancers)), axes=F, main="BG3 Predicted Enhancer Widths", xlab="Width (bp)", col=cbbPalette[3], cex.main=2, cex.lab=1.5, xlim = c(3,14), ylim = c(0,3000))
axis(1, at=c(log2(10), log2(50), log2(150), log2(400), log2(1000), log2(4000), log2(16000)),
    labels=c(10, 50, 150, 400, 1000, 4000, 16000))
axis(2, at=seq(0,3000, by = 500), labels = seq(0,3000, by = 500))
abline(v = log2(50), col="red", lty=2, lwd=2)
abline(v = log2(1000), col="red", lty=2, lwd=2)
text(4.25, 2925, "Fragments")
text(7.75, 2925, "Enhancers")
text(12, 2925, "Super\nEnhancers")
dev.off()

pdf("/home/jw18713/Project1/Paper_Plots/Figure2/S2_Enhancer_Widths.pdf")
hist(log2(width(S2_Enhancers)), axes=F, main="S2 Predicted Enhancer Widths", xlab="Width (bp)", col=cbbPalette[3], cex.main=2, cex.lab=1.5, xlim = c(3,14), ylim = c(0,3000))
axis(1, at=c(log2(10), log2(50), log2(150), log2(400), log2(1000), log2(4000), log2(16000)),
    labels=c(10, 50, 150, 400, 1000, 4000, 16000))
axis(2, at=seq(0,3000, by = 500), labels = seq(0,3000, by = 500))
abline(v = log2(50), col="red", lty=2, lwd=2)
abline(v = log2(1000), col="red", lty=2, lwd=2)
text(4.25, 2925, "Fragments")
text(7.75, 2925, "Enhancers")
text(12, 2925, "Super\nEnhancers")
dev.off()


BG3_fuzzy_only <- BG3_enhancers[!BG3_enhancers %over% BG3_starr]
BG3_shared <- BG3_enhancers[BG3_enhancers %over% BG3_starr]

S2_fuzzy_only <- S2_enhancers[!S2_enhancers %over% S2_starr]
S2_shared <- S2_enhancers[S2_enhancers %over% S2_starr]

# Novel Enhancers

pdf("/home/jw18713/Project1/Paper_Plots/Figure2/Figure2E_Novel.pdf")
hist(log2(width(BG3_fuzzy_only)), axes=F, main="BG3 XAI Only Enhancer Widths", xlab="Width (bp)", col=cbbPalette[3], cex.main=2, cex.lab=1.5, xlim = c(3,14), ylim = c(0,3000), breaks = 20)
axis(1, at=c(log2(10), log2(50), log2(150), log2(400), log2(1000), log2(4000), log2(16000)),
    labels=c(10, 50, 150, 400, 1000, 4000, 16000))
axis(2, at=seq(0,3000, by = 500), labels = seq(0,3000, by = 500))
abline(v = log2(50), col="red", lty=2, lwd=2)
abline(v = log2(1000), col="red", lty=2, lwd=2)
text(4.25, 2925, "Fragments")
text(7.75, 2925, "Enhancers")
text(12, 2925, "Super\nEnhancers")
dev.off()

pdf("/home/jw18713/Project1/Paper_Plots/Figure2/Figure2F_Novel.pdf")
hist(log2(width(S2_fuzzy_only)), axes=F, main="S2 XAI Only Enhancer Widths", xlab="Width (bp)", col=cbbPalette[3], cex.main=2, cex.lab=1.5, xlim = c(3,14), ylim = c(0,3000), breaks = 20)
axis(1, at=c(log2(10), log2(50), log2(150), log2(400), log2(1000), log2(4000), log2(16000)),
    labels=c(10, 50, 150, 400, 1000, 4000, 16000))
axis(2, at=seq(0,3000, by = 500), labels = seq(0,3000, by = 500))
abline(v = log2(50), col="red", lty=2, lwd=2)
abline(v = log2(1000), col="red", lty=2, lwd=2)
text(4.25, 2925, "Fragments")
text(7.75, 2925, "Enhancers")
text(12, 2925, "Super\nEnhancers")
dev.off()

# Shared Enhancers

pdf("/home/jw18713/Project1/Paper_Plots/Figure2/Figure2E_Shared.pdf")
hist(log2(width(BG3_shared)), axes=F, main="BG3 Both Enhancer Widths", xlab="Width (bp)", col=cbbPalette[3], cex.main=2, cex.lab=1.5, xlim = c(3,14), ylim = c(0,400), breaks = 20)
axis(1, at=c(log2(10), log2(50), log2(150), log2(400), log2(1000), log2(4000), log2(16000)),
    labels=c(10, 50, 150, 400, 1000, 4000, 16000))
axis(2, at=seq(0,400, by = 50), labels = seq(0,400, by = 50))
abline(v = log2(50), col="red", lty=2, lwd=2)
abline(v = log2(1000), col="red", lty=2, lwd=2)
text(4.25, 390, "Fragments")
text(7.75, 390, "Enhancers")
text(12, 390, "Super\nEnhancers")
dev.off()

pdf("/home/jw18713/Project1/Paper_Plots/Figure2/Figure2F_Shared.pdf")
hist(log2(width(S2_shared)), axes=F, main="S2 Both Enhancer Widths", xlab="Width (bp)", col=cbbPalette[3], cex.main=2, cex.lab=1.5, xlim = c(3,14), ylim = c(0,400), breaks = 20)
axis(1, at=c(log2(10), log2(50), log2(150), log2(400), log2(1000), log2(4000), log2(16000)),
    labels=c(10, 50, 150, 400, 1000, 4000, 16000))
axis(2, at=seq(0,400, by = 50), labels = seq(0,400, by = 50))
abline(v = log2(50), col="red", lty=2, lwd=2)
abline(v = log2(1000), col="red", lty=2, lwd=2)
text(4.25, 390, "Fragments")
text(7.75, 390, "Enhancers")
text(12, 390, "Super\nEnhancers")
dev.off()
