library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)
library(VennDiagram)

S2_enhancers <- get(load("grown_predicted_enhancers_S2.Rda"))
S2_starr <- get(load("starr_S2.Rda"))
BG3_enhancers <- get(load("grown_predicted_enhancers_BG3.Rda"))
BG3_starr <- get(load("starr.Rda"))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

generatePairwiseVenn <- function(set1, set2, categories, cols=c("#0072B2", "#D55E00"), cat.pos = c(-20, 20), cat.dist = c(0.07, 0.07), human = TRUE){

  v <- draw.pairwise.venn(area1 = length(set1), area2 = length(set2), cross.area = sum(set1%over%set2),
                          category = categories, col = "transparent", fill = cols,
                          alpha = 0.6, label.col = rep("black", 3), cex = 1.2,
                          cat.col =  cols, cat.cex = 1.4,
                          cat.pos = cat.pos,
                          cat.dist = cat.dist,
                          margin = 0.2,
                          euler.d =TRUE, scaled = T
  )
  if(human){
    for(i in 5:7){
      v[[i]]$label  <- as.vector(sciNotation(as.numeric(v[[i]]$label)))
    }
  }

  return(v)
}

pdf("/home/jw18713/Project1/Paper_Plots/Figure1/BG3_Fuzzy_vs_STARR_Venn.pdf")
plot.new()
v1 <- generatePairwiseVenn(BG3_enhancers, BG3_starr, c("XAi", "STARR-seq"),
                           cols=c(cbbPalette[6], cbbPalette[7]),
                           cat.pos = c(-30, 150),
                           cat.dist = c(0.075, 0.075), human = F)
title(main = "BG3")
dev.off()

pdf("/home/jw18713/Project1/Paper_Plots/Figure1/S2_Fuzzy_vs_STARR_Venn.pdf")
plot.new()
v1 <- generatePairwiseVenn(S2_enhancers, S2_starr, c("XAi", "STARR-seq"),
                           cols=c(cbbPalette[6], cbbPalette[7]),
                           cat.pos = c(-30, 150),
                           cat.dist = c(0.075, 0.075), human = F)
title(main = "S2")
dev.off()
