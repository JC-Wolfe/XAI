library(GenomicRanges)
library(rtracklayer)

load("~/dm6_promoter_correction/Expression_lists/S2_background_proximal.Rda")
load("~/dm6_promoter_correction/Expression_lists/S2_background_distal_only.Rda")
load("~/dm6_promoter_correction/Expression_lists/BG3_background_proximal.Rda")
load("~/dm6_promoter_correction/Expression_lists/BG3_background_distal_only.Rda")
load("~/dm6_promoter_correction/Expression_lists/S2_putative_proximal.Rda")
load("~/dm6_promoter_correction/Expression_lists/S2_putative_distal_only.Rda")
load("~/dm6_promoter_correction/Expression_lists/BG3_putative_proximal.Rda")
load("~/dm6_promoter_correction/Expression_lists/BG3_putative_distal_only.Rda")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7")

# Ok, so this is my histogram then
pdf(file = "~/dm6_promoter_correction/Redone_plots/Figure3D.pdf",
width = 14, height = 14, pointsize = 14)
par(mar=c(4.9,4.5,5.5,9), mfrow=c(2,2), cex = 1.2)
par(xpd = NA)

hist(log2(BG3_putative_distal_only$max_Expression+1),
  main="BG3 Putative Enhancers", xlab="maximum expression of target genes (FPKM)",
  col="#0072B280", cex.main=2, cex.lab=1, axes=F,
  ylim = c(0,3000), breaks = 16, ylab = NA)
axis(1, at=seq(0,16,by=4), labels=c(0,2^seq(4,16, by = 4)))
axis(2, at=seq(0,3000, by = 1000), labels = seq(0,3000, by = 1000), las=2)

hist(log2(BG3_putative_proximal$max_Expression+1), col="#D55E0080",
add=T, breaks = 16)

# S2 stuff

hist(log2(S2_putative_distal_only$max_Expression+1),
  main="S2 Putative Enhancers",
  xlab="maximum expression of target genes (FPKM)", col="#0072B280", cex.main=2, cex.lab=1, axes=F,
  ylim = c(0,3000), breaks = 16, ylab = NA)
axis(1, at=seq(0,16,by=4), labels=c(0,2^seq(4,16, by = 4)))
axis(2, at=seq(0,3000, by = 1000), labels = seq(0,3000, by = 1000), las=2)

hist(log2(S2_putative_proximal$max_Expression+1), col="#D55E0080",
  add=T, breaks = 16)

legend(16.5,1875,c("Distal only", "Proximal"), fill=c("#0072B280", "#D55E0080"),
  bty='n')

#--------------------------Background Expression-------------------------------#
#------------------------------------------------------------------------------#


hist(log2(BG3_background_distal_only$max_Expression+1),
  main="BG3 whole genome",
  xlab="maximum expression of target genes (FPKM)", col="#56B4E980", cex.main=2, cex.lab=1, axes=F,
  ylim = c(0,60000), breaks = 16, ylab = NA)
axis(1, at=seq(0,16,by=4), labels=c(0,2^seq(4,16, by = 4)))
axis(2, at=seq(0,60000, by = 15000), labels = seq(0,60000, by = 15000), las=2)

hist(log2(BG3_background_proximal$max_Expression+1), col="#E69F0080",
add=T, breaks = 16)

# S2 stuff

hist(log2(S2_background_distal_only$max_Expression+1),
  main="S2 whole genome",
  xlab="maximum expression of target genes (FPKM)", col="#56B4E980", cex.main=2, cex.lab=1, axes=F,
  ylim = c(0,60000), breaks = 16, ylab = NA)
axis(1, at=seq(0,16,by=4), labels=c(0,2^seq(4,16, by = 4)))
axis(2, at=seq(0,60000, by = 15000), labels = seq(0,60000, by = 15000), las=2)

hist(log2(S2_background_proximal$max_Expression+1), col="#E69F0080",
  add=T, breaks = 16)

legend(16.5,(37500),c("Distal only", "Proximal"), fill=c("#56B4E980", "#E69F0080"),
  bty='n')


dev.off()

wilcox.test(log2(BG3_putative_distal_only$max_Expression+1),
  log2(BG3_putative_proximal$max_Expression+1))

wilcox.test(log2(S2_putative_distal_only$max_Expression+1),
  log2(S2_putative_proximal$max_Expression+1))

wilcox.test(log2(BG3_background_distal_only$max_Expression+1),
  log2(BG3_background_proximal$max_Expression+1))

wilcox.test(log2(S2_background_distal_only$max_Expression+1),
  log2(S2_background_proximal$max_Expression+1))
