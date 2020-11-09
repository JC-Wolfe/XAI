install.packages("PRROC")
library(PRROC)
# Load the required packages
require(ggplot2)
require(ggseqlogo)

read_LG_output <- function(x){
  out_df <- read.delim(x, sep = "\t", stringsAsFactors = F)
  out_gr <- GRanges(seqnames = out_df$seqnames,
    ranges = IRanges(start = out_df$start, end = out_df$end),
    strand = out_df$strand)
  out_gr$Developmental <- out_df$Developmental
  out_gr$Predicted_Class <- out_df$Predicted.Class
  out_gr$Conf_Perc_1 <- out_df$Conf..Perc.1
  out_gr$Conf_Perc_0 <- out_df$Conf..Perc.0
  out_gr$Errors <- out_df$Errors
  return(out_gr)
}

DvHK_XAI <- read_LG_output("dev_hk_results.tsv")


fg1 <- DvHK_XAI$Conf_Perc_1[DvHK_XAI$Developmental == 1]
bg1 <- DvHK_XAI$Conf_Perc_1[DvHK_XAI$Developmental == 0]

# DvHK ROC Curve
DvHK_roc <- roc.curve(scores.class0 = fg1, scores.class1 = bg1, curve = T)

# DvHK PR Curve
DvHK_pr <- pr.curve(scores.class0 = fg1, scores.class1 = bg1, curve = T)



png("/home/jw18713/Project1/Paper_Plots/FigureS1/DvHk_roc.png", width=500, height=500)
plot(DvHK_roc, main = "Developmental vs Housekeeping ROC")
dev.off()
png("/home/jw18713/Project1/Paper_Plots/FigureS1/DvHk_pr.png", width=500, height=500)
plot(DvHK_pr, main = "Developmental vs Housekeeping PRC")
dev.off()


# VennDiagram package
if(!require("VennDiagram", character.only = TRUE)){
  install.packages("VennDiagram")
}
library(VennDiagram)

# install package gridExtra if not already installed
if(!require("gridExtra", character.only = TRUE)){
  install.packages("gridExtra")
}
library(gridExtra)

# install package grid if not already installed
if(!require("grid", character.only = TRUE)){
  install.packages("grid")
}
library(grid)

# install package gridBase if not already installed
if(!require("gridBase", character.only = TRUE)){
  install.packages("gridBase")
}
library(gridBase)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

generateTripleVenn <- function(set1, set2, set3, categories, cols=c("#0072B2", "#D55E00", "#CC79A7"), cat.pos = c(-20, 0, 20), cat.dist = c(0.07, 0.07, 0.07), human = TRUE){

  v <- draw.triple.venn(area1 = length(set1), area2 = length(set2), area3 = length(set3),
                          n12 = sum(set1%over%set2),
                          n23 = sum(set2%over%set3),
                          n13 = sum(set1%over%set3),
                          n123 = sum(set1%over%set2 & set1%over%set3),
                          category = categories, col = "transparent", fill = cols,
                          alpha = 0.6, label.col = rep("black", 7), cex = 1.2,
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



S2_enhancers_neither <- get(load("~/Archive_Year1/S2_novel_neither.Rda"))
S2_SE_neither <- S2_enhancers_neither[width(S2_enhancers_neither) >= 1000]
S2_proximal_enhancers <- get(load("~/Archive_Year1/S2_novel_proximal.Rda"))
S2_proximal_SE <- S2_proximal_enhancers[width(S2_proximal_enhancers) >= 1000]
S2_distal_enhancers <- get(load("~/Archive_Year1/S2_novel_distal_all.Rda"))
S2_distal_SE <- S2_distal_enhancers[width(S2_distal_enhancers) >= 1000]

dev <- reduce(DvHK_XAI[DvHK_XAI$Developmental == 1])
hk <- reduce(DvHK_XAI[DvHK_XAI$Developmental == 0])
SE <- c(S2_distal_SE, S2_proximal_SE, S2_SE_neither)

v1 <- generateTripleVenn(dev, hk, SE, c("Developmental\nEnhancers", "Housekeeping\nEnhancers", "Super Enhancers"),
                           cols=c(cbbPalette[6], cbbPalette[7], cbbPalette[8]),
                           cat.pos = c(-60, 60, 180),
                           cat.dist = c(0.125, 0.125, 0.05),
                           human = F)

pdf("~/Project1/Paper_Plots/FigureS1/Dev_Hk_SE_Venn.pdf", width=5, height=5,pointsize = 10);
par(mar=c(0, 0, 0, 0)+0.1)
generateTripleVenn(dev, hk, SE, c("Developmental\nEnhancers", "Housekeeping\nEnhancers", "Super Enhancers"),
                           cols=c(cbbPalette[6], cbbPalette[7], cbbPalette[8]),
                           cat.pos = c(-150, 0, 150),
                           cat.dist = c(0.2, 0.2, 0.2),
                           human = F)
dev.off()
