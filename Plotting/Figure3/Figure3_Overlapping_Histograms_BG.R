library(GenomicRanges)
library(rtracklayer)

load("/home/jw18713/rerun_glists/putative/BG3_put_dist.Rda")
load("/home/jw18713/rerun_glists/putative/BG3_put_prox.Rda")
load("/home/jw18713/rerun_glists/bg/BG3_bg_prox.Rda")
load("/home/jw18713/rerun_glists/bg/BG3_bg_dist.Rda")


load("/home/jw18713/rerun_glists/putative/S2_put_dist.Rda")
load("/home/jw18713/rerun_glists/putative/S2_put_prox.Rda")
load("/home/jw18713/rerun_glists/bg/S2_bg_prox.Rda")
load("/home/jw18713/rerun_glists/bg/S2_bg_dist.Rda")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Proximal BG3 (BG3P)

# Now which of these don't have any expresion?
no_expression_BG3P <- BG3_put_prox[sapply(BG3_put_prox, length) == 0]

# And which of these are expressed?
enhancers_contacted_BG3P <- BG3_put_prox[sapply(BG3_put_prox, length) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$FPKM, na.rm = T))
}

max_to_plot_BG3P <- sapply(enhancers_contacted_BG3P, maxExp)

meanExp <- function(gr){
  return(mean(gr$FPKM, na.rm = T))
}

means_plot_BG3P <- sapply(enhancers_contacted_BG3P, meanExp)

no_expression_max_BG3P <- c(no_expression_BG3P, enhancers_contacted_BG3P[max_to_plot_BG3P == 0])
expression_max_BG3P <- max_to_plot_BG3P[max_to_plot_BG3P > 0]
expression_means_BG3P <- means_plot_BG3P[means_plot_BG3P > 0]

# DISTAL ONLY BG3 (BG3D)

# Now which of these don't have any expresion?
no_expression_BG3D <- BG3_put_dist[sapply(BG3_put_dist, length) == 0]

# And which of these are expressed?
enhancers_contacted_BG3D <- BG3_put_dist[sapply(BG3_put_dist, length) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$FPKM, na.rm = T))
}

max_to_plot_BG3D <- sapply(enhancers_contacted_BG3D, maxExp)

meanExp <- function(gr){
  return(mean(gr$FPKM, na.rm = T))
}

means_plot_BG3D <- sapply(enhancers_contacted_BG3D, meanExp)

no_expression_max_BG3D <- c(no_expression_BG3D, enhancers_contacted_BG3D[max_to_plot_BG3D == 0])
expression_max_BG3D <- max_to_plot_BG3D[max_to_plot_BG3D > 0]
expression_means_BG3D <- means_plot_BG3D[means_plot_BG3D > 0]

# DISTAL ONLY BG3 background (BG3_bg_do)

# Now which of these don't have any expresion?
no_expression_BG3_bg_do <- BG3_bg_dist[sapply(BG3_bg_dist, length) == 0]

# And which of these are expressed?
enhancers_contacted_BG3_bg_do <- BG3_bg_dist[sapply(BG3_bg_dist, length) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$FPKM, na.rm = T))
}

max_to_plot_BG3_bg_do <- sapply(enhancers_contacted_BG3_bg_do, maxExp)

meanExp <- function(gr){
  return(mean(gr$FPKM, na.rm = T))
}

means_plot_BG3_bg_do <- sapply(enhancers_contacted_BG3_bg_do, meanExp)

no_expression_max_BG3_bg_do <- c(no_expression_BG3_bg_do, enhancers_contacted_BG3_bg_do[max_to_plot_BG3_bg_do == 0])
expression_max_BG3_bg_do <- max_to_plot_BG3_bg_do[max_to_plot_BG3_bg_do > 0]
expression_means_BG3_bg_do <- means_plot_BG3_bg_do[means_plot_BG3_bg_do > 0]

# BG3_background_proximal (BG3_bg_p)

# Now which of these don't have any expresion?
no_expression_BG3_bg_p <- BG3_bg_prox[sapply(BG3_bg_prox, length) == 0]

# And which of these are expressed?
enhancers_contacted_BG3_bg_p <- BG3_bg_prox[sapply(BG3_bg_prox, length) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$FPKM, na.rm = T))
}

max_to_plot_BG3_bg_p <- sapply(enhancers_contacted_BG3_bg_p, maxExp)

meanExp <- function(gr){
  return(mean(gr$FPKM, na.rm = T))
}

means_plot_BG3_bg_p <- sapply(enhancers_contacted_BG3_bg_p, meanExp)

no_expression_max_BG3_bg_p <- c(no_expression_BG3_bg_p, enhancers_contacted_BG3_bg_p[max_to_plot_BG3_bg_p == 0])
expression_max_BG3_bg_p <- max_to_plot_BG3_bg_p[max_to_plot_BG3_bg_p > 0]
expression_means_BG3_bg_p <- means_plot_BG3_bg_p[means_plot_BG3_bg_p > 0]


# Proximal S2 (S2P)

# Now which of these don't have any expresion?
no_expression_S2P <- S2_put_prox[sapply(S2_put_prox, length) == 0]

# And which of these are expressed?
enhancers_contacted_S2P <- S2_put_prox[sapply(S2_put_prox, length) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$FPKM, na.rm = T))
}

max_to_plot_S2P <- sapply(enhancers_contacted_S2P, maxExp)

meanExp <- function(gr){
  return(mean(gr$FPKM, na.rm = T))
}

means_plot_S2P <- sapply(enhancers_contacted_S2P, meanExp)

no_expression_max_S2P <- c(no_expression_S2P, enhancers_contacted_S2P[max_to_plot_S2P == 0])
expression_max_S2P <- max_to_plot_S2P[max_to_plot_S2P > 0]
expression_means_S2P <- means_plot_S2P[means_plot_S2P > 0]


# DISTAL ONLY S2 (S2D)

# Now which of these don't have any expresion?
no_expression_S2D <- S2_put_dist[sapply(S2_put_dist, length) == 0]

# And which of these are expressed?
enhancers_contacted_S2D <- S2_put_dist[sapply(S2_put_dist, length) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$FPKM, na.rm = T))
}

max_to_plot_S2D <- sapply(enhancers_contacted_S2D, maxExp)

meanExp <- function(gr){
  return(mean(gr$FPKM, na.rm = T))
}

means_plot_S2D <- sapply(enhancers_contacted_S2D, meanExp)

no_expression_max_S2D <- c(no_expression_S2D, enhancers_contacted_S2D[max_to_plot_S2D == 0])
expression_max_S2D <- max_to_plot_S2D[max_to_plot_S2D > 0]
expression_means_S2D <- means_plot_S2D[means_plot_S2D > 0]

# S2 background distal only (S2_bg_do)

# Now which of these don't have any expresion?
no_expression_S2_bg_do <- S2_bg_dist[sapply(S2_bg_dist, length) == 0]

# And which of these are expressed?
enhancers_contacted_S2_bg_do <- S2_bg_dist[sapply(S2_bg_dist, length) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$FPKM, na.rm = T))
}

max_to_plot_S2_bg_do <- sapply(enhancers_contacted_S2_bg_do, maxExp)

meanExp <- function(gr){
  return(mean(gr$FPKM, na.rm = T))
}

means_plot_S2_bg_do <- sapply(enhancers_contacted_S2_bg_do, meanExp)

no_expression_max_S2_bg_do <- c(no_expression_S2_bg_do, enhancers_contacted_S2_bg_do[max_to_plot_S2_bg_do == 0])
expression_max_S2_bg_do <- max_to_plot_S2_bg_do[max_to_plot_S2_bg_do > 0]
expression_means_S2_bg_do <- means_plot_S2_bg_do[means_plot_S2_bg_do > 0]

# DISTAL ONLY S2 (S2_bg_p)

# Now which of these don't have any expresion?
no_expression_S2_bg_p <- S2_bg_prox[sapply(S2_bg_prox, length) == 0]

# And which of these are expressed?
enhancers_contacted_S2_bg_p <- S2_bg_prox[sapply(S2_bg_prox, length) > 0]

# Max Expression Data
maxExp <- function(gr){
  return(max(gr$FPKM, na.rm = T))
}

max_to_plot_S2_bg_p <- sapply(enhancers_contacted_S2_bg_p, maxExp)

meanExp <- function(gr){
  return(mean(gr$FPKM, na.rm = T))
}

means_plot_S2_bg_p <- sapply(enhancers_contacted_S2_bg_p, meanExp)

no_expression_max_S2_bg_p <- c(no_expression_S2_bg_p, enhancers_contacted_S2_bg_p[max_to_plot_S2_bg_p == 0])
expression_max_S2_bg_p <- max_to_plot_S2_bg_p[max_to_plot_S2_bg_p > 0]
expression_means_S2_bg_p <- means_plot_S2_bg_p[means_plot_S2_bg_p > 0]

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Ok, so this is my histogram then
pdf(file = "/home/jw18713/Project1/Paper_Plots/Figure3/Figure3D_max.pdf",
width = 14, height = 12, pointsize = 14)
par(mar=c(4.9,5.5,5.5,8), mfrow=c(2,2), cex = 1.2)
par(xpd = NA)

hist(log2(expression_max_BG3D+1), main="BG3 Enhancers\nMax Contacted Expression",
  xlab="Expression FPKM", col="#0072B280", cex.main=2, cex.lab=1.5, axes=F,
  ylim = c(0,3000), breaks = 16, ylab = NA)
axis(1, at=seq(0,16,by=4), labels=c(0,2^seq(4,16, by = 4)))
axis(2, at=seq(0,3000, by = 1000), labels = seq(0,3000, by = 1000), las=2)

hist(log2(expression_max_BG3P+1), col="#D55E0080", add=T, breaks = 16)

legend(16.5,2000,c("Distal", "Proximal"), fill=c("#0072B280", "#D55E0080"), bty='n')

# S2 stuff

hist(log2(expression_max_S2D+1), main="S2 Enhancers\nMax Contacted Expression",
  xlab="Expression FPKM", col="#0072B280", cex.main=2, cex.lab=1.5, axes=F,
  ylim = c(0,3000), breaks = 16, ylab = NA)
axis(1, at=seq(0,16,by=4), labels=c(0,2^seq(4,16, by = 4)))
axis(2, at=seq(0,3000, by = 1000), labels = seq(0,3000, by = 1000), las=2)

hist(log2(expression_max_S2P+1), col="#D55E0080", add=T, breaks = 16)

legend(16.5,2000,c("Distal", "Proximal"), fill=c("#0072B280", "#D55E0080"), bty='n')

# And background maximums
# Ok, so this is my histogram then

hist(log2(expression_max_BG3_bg_do+1), main="BG3 Background",
  xlab="Expression FPKM", col="#56B4E980", cex.main=2, cex.lab=1.5, axes=F,
  ylim = c(0,60000), breaks = 16, ylab = NA)
axis(1, at=seq(0,16,by=4), labels=c(0,2^seq(4,16, by = 4)))
axis(2, at=seq(0,60000, by = 15000), labels = seq(0,60000, by = 15000), las=2)

hist(log2(expression_max_BG3_bg_p+1), col="#E69F0080", add=T, breaks = 16)

legend(16.5,2000,c("Distal", "Proximal"), fill=c("#56B4E980", "#E69F0080"), bty='n')

hist(log2(expression_max_S2_bg_do+1), main="S2 Background",
  xlab="Expression FPKM", col="#56B4E980", cex.main=2, cex.lab=1.5, axes=F,
  ylim = c(0,60000), breaks = 16, ylab = NA)
axis(1, at=seq(0,16,by=4), labels=c(0,2^seq(4,16, by = 4)))
axis(2, at=seq(0,60000, by = 15000), labels = seq(0,60000, by = 15000), las=2)

hist(log2(expression_max_S2_bg_p+1), col="#E69F0080", add=T, breaks = 16)

legend(16.5,40000,c("Distal", "Proximal"), fill=c("#56B4E980", "#E69F0080"), bty='n')

dev.off()




# Ok, so this is my histogram then
pdf(file = "/home/jw18713/Project1/Paper_Plots/Figure3/Figure3_mean.pdf",
width = 14, height = 12, pointsize = 14)
par(mar=c(4.9,5.5,5.5,8), mfrow=c(2,2), cex = 1.2)
par(xpd = NA)

hist(log2(expression_means_BG3D+1), main="BG3 Enhancers\nMean Contacted Expression",
  xlab="Expression FPKM", col="#0072B280", cex.main=2, cex.lab=1.5, axes=F,
  ylim = c(0,3000), breaks = seq(0,16,by=1), ylab = NA)
axis(1, at=seq(0,16,by=4), labels=c(0,2^seq(4,16, by = 4)))
axis(2, at=seq(0,3000, by = 1000), labels = seq(0,3000, by = 1000), las=2)

hist(log2(expression_means_BG3P+1), col="#D55E0080", add=T, breaks = seq(0,16,by=1))

legend(16.5,2000,c("Distal", "Proximal"), fill=c("#0072B280", "#D55E0080"), bty='n')

# S2 stuff

hist(log2(expression_means_S2D+1), main="S2 Enhancers\nMean Contacted Expression",
  xlab="Expression FPKM", col="#0072B280", cex.main=2, cex.lab=1.5, axes=F,
  ylim = c(0,3000), breaks = seq(0,16,by=1), ylab = NA)
axis(1, at=seq(0,16,by=4), labels=c(0,2^seq(4,16, by = 4)))
axis(2, at=seq(0,3000, by = 1000), labels = seq(0,3000, by = 1000), las=2)

hist(log2(expression_means_S2P+1), col="#D55E0080", add=T, breaks = seq(0,16,by=1))

legend(16.5,2000,c("Distal", "Proximal"), fill=c("#0072B280", "#D55E0080"), bty='n')

# And background means
# Ok, so this is my histogram then

hist(log2(expression_means_BG3_bg_do+1), main="BG3 Background",
  xlab="Expression FPKM", col="#56B4E980", cex.main=2, cex.lab=1.5, axes=F,
  ylim = c(0,60000), breaks = seq(0,16,by=1), ylab = NA)
axis(1, at=seq(0,16,by=4), labels=c(0,2^seq(4,16, by = 4)))
axis(2, at=seq(0,60000, by = 15000), labels = seq(0,60000, by = 15000), las=2)

hist(log2(expression_means_BG3_bg_p+1), col="#E69F0080", add=T, breaks = seq(0,16,by=1))

legend(16.5,40000,c("Distal", "Proximal"), fill=c("#56B4E980", "#E69F0080"), bty='n')

hist(log2(expression_means_S2_bg_do+1), main="S2 Background",
  xlab="Expression FPKM", col="#56B4E980", cex.main=2, cex.lab=1.5, axes=F,
  ylim = c(0,60000), breaks = seq(0,16,by=1), ylab = NA)
axis(1, at=seq(0,16,by=4), labels=c(0,2^seq(4,16, by = 4)))
axis(2, at=seq(0,60000, by = 15000), labels = seq(0,60000, by = 15000), las=2)

hist(log2(expression_means_S2_bg_p+1), col="#E69F0080", add=T, breaks = seq(0,16,by=1))

legend(16.5,40000,c("Distal", "Proximal"), fill=c("#56B4E980", "#E69F0080"), bty='n')

dev.off()

# Stats
wilcox.test(log2(expression_max_BG3_bg_p+1), log2(expression_max_BG3_bg_do+1))
wilcox.test(log2(expression_max_S2_bg_p+1), log2(expression_max_S2_bg_do+1))
wilcox.test(log2(expression_max_BG3D+1), log2(expression_max_BG3P+1))
wilcox.test(log2(expression_max_S2D+1), log2(expression_max_S2P+1))
