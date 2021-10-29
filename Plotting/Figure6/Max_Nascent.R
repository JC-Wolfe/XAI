
library(genomation)
library(GenomicRanges)
library(rtracklayer)

tracks <- read.csv("~/Review_Paper_Plots/Mila_boxplots/Enhancer_csv_new.csv")
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

common <- tracks[tracks$IsCommon == 1,]
putative <- tracks[tracks$IsCommon == 0,]

background <- read.csv("~/Review_Paper_Plots/Mila_boxplots/Background_new.csv")

stats_matrix <- matrix(0, 3, dimnames = list(
  c("Putative/Common", "Putative/Background", "Common/Background"),
  c("X3NTseq_max")
))

setwd("~/dm6_promoter_correction/Figure6_stuff")

#---------------------------------Max Nascent----------------------------------#
#------------------------------------------------------------------------------#

pdf("Max_X3NTseq_boxplot.pdf", width = 8, height = 6, pointsize = 14)
par(mar = c(4,5.5,4,8), xpd = T)

boxplot(putative$X3NTseq_max, common$X3NTseq_max, background$X3NTseq_max,
  names = c("Putative", "Common", "Background"),
  ylim = c(0, 10000),
  outline = F, col = cbbPalette[2:4], main = "3'NT-seq")
dev.off()

  # Stats code
stats_matrix[1,1] <- wilcox.test(putative$X3NTseq_max,
  common$X3NTseq_max)$p.value
stats_matrix[2,1] <- wilcox.test(putative$X3NTseq_max,
  background$X3NTseq_max)$p.value
stats_matrix[3,1] <- wilcox.test(common$X3NTseq_max,
  background$X3NTseq_max)$p.value

write.csv(stats_matrix, file = "Nascent_stats_matrix.csv")
