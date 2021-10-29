library(pROC)


setwd("~/dm6_promoter_correction/Figure6_stuff")

Arch <- read.csv("Fig6_Model_Results.tsv", sep = "\t")

v1 <- roc(Arch$IsCommon, Arch$Conf..Perc.1)

setwd("~/dm6_promoter_correction/Figure6_stuff")
pdf("B&W_roc.pdf")
plot(v1, main = "Common vs Putative ROC")
dev.off()
