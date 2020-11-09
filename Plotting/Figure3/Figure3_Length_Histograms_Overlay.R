
library(GenomicRanges)
library(rtracklayer)

load("Archive_Year1/BG3_novel_proximal.Rda")
load("Archive_Year1/BG3_novel_distal.Rda")
load("Archive_Year1/S2_novel_proximal.Rda")
load("Archive_Year1/S2_novel_distal.Rda")

pdf(file = "/home/jw18713/Project1/Paper_Plots/Figure3/Figure3C&D.pdf",
width = 18, height = 8, pointsize = 14)
par(mar=c(5,5.5,5,6), mfrow=c(1,2), cex = 1.2)
par(xpd = NA)

      # BG3 proximal
      hist(log2(width(BG3_novel_proximal)), axes=F, main="BG3 Novel Enhancer Widths",
        xlab="Width (bp)", col="#D55E0080", cex.main=2, cex.lab=1.5, xlim = c(3,14), ylim = c(0,1500),
        breaks = 20, ylab=NA)
          axis(1, at=c(log2(10), log2(50), log2(150), log2(400), log2(1000), log2(4000), log2(16000)),
              labels=c(10, 50, 150, 400, 1000, 4000, 16000))
          axis(2, at=seq(0,1500, by = 250), labels = seq(0,1500, by = 250), las=2)
          lines(x=c(log2(50),log2(50)), y=c(0,1500), col="#D55E00", lty=2, lwd=2)
          lines(x=c(log2(1000),log2(1000)), y=c(0,1500), col="#D55E00", lty=2, lwd=2)
          text(4.25, 1425, "Fragments")
          text((log2(1000) - log2(50))/2 + log2(50), 1425, "Enhancers")
          text((14 - log2(1000))/2 + log2(1000), 1425, "Super\nEnhancers")

      #BG3 distal only
      hist(log2(width(BG3_novel_distal)), col="#0072B280", add= T, breaks = 20)

      legend(14.5,1000,c("Distal\nOnly", "Proximal"), fill=c("#0072B280", "#D55E0080"), bty="n")

      #S2 proximal
      hist(log2(width(S2_novel_proximal)), axes=F, main="S2 Novel Enhancer Widths",
        xlab="Width (bp)", col="#D55E0080", cex.main=2, cex.lab=1.5, xlim = c(3,14), ylim = c(0,1500),
        breaks = 20, ylab=NA)
          axis(1, at=c(log2(10), log2(50), log2(150), log2(400), log2(1000), log2(4000), log2(16000)),
              labels=c(10, 50, 150, 400, 1000, 4000, 16000))
          axis(2, at=seq(0,1500, by = 250), labels = seq(0,1500, by = 250), las=2)
          lines(x=c(log2(50),log2(50)), y=c(0,1500), col="#D55E00", lty=2, lwd=2)
          lines(x=c(log2(1000),log2(1000)), y=c(0,1500), col="#D55E00", lty=2, lwd=2)
          text(4.25, 1425, "Fragments")
          text((log2(1000) - log2(50))/2 + log2(50), 1425, "Enhancers")
          text((14 - log2(1000))/2 + log2(1000), 1425, "Super\nEnhancers")

      #S2 distal only
      hist(log2(width(S2_novel_distal)), col="#0072B280", add = T, breaks = 20)
      legend(14.5,1000,c("Distal\nOnly", "Proximal"), fill=c("#0072B280", "#D55E0080"), bty="n")

dev.off()
