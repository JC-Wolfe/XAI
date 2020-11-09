install.packages("PRROC")
library(PRROC)
# Load the required packages
require(ggplot2)
require(ggseqlogo)

BG3_XAI <- get(load("/home/jw18713/Project1/data/model_predictions/BG3_FuzzyGR.Rda"))
S2_XAI <- get(load("/home/jw18713/Project1/data/model_predictions/S2_FuzzyGR.Rda"))

fg1 <- BG3_XAI$Conf..Perc.1[BG3_XAI$STARR_seq_binary == 1]
bg1 <- BG3_XAI$Conf..Perc.1[BG3_XAI$STARR_seq_binary == 0]

# BG3 ROC Curve
BG3_roc <- roc.curve(scores.class0 = fg1, scores.class1 = bg1, curve = T)

# BG3 PR Curve
BG3_pr <- pr.curve(scores.class0 = fg1, scores.class1 = bg1, curve = T)

fg2 <- S2_XAI$Conf..Perc.1[S2_XAI$STARR_seq_binary == 1]
bg2 <- S2_XAI$Conf..Perc.1[S2_XAI$STARR_seq_binary == 0]

# S2 ROC Curve
S2_roc <- roc.curve(scores.class0 = fg2, scores.class1 = bg2, curve = T)

# S2 PR Curve
S2_pr <- pr.curve(scores.class0 = fg2, scores.class1 = bg2, curve = T)



png("/home/jw18713/Project1/Paper_Plots/FigureS1/Supp_Figure1A.png", width=500, height=500)
plot(BG3_roc, main = "BG3 ROC")
dev.off()
png("/home/jw18713/Project1/Paper_Plots/FigureS1/Supp_Figure1B.png", width=500, height=500)
plot(BG3_pr, main = "BG3 PRC")
dev.off()
png("/home/jw18713/Project1/Paper_Plots/FigureS1/Supp_Figure1C.png", width=500, height=500)
plot(S2_roc, main = "S2 ROC")
dev.off()
png("/home/jw18713/Project1/Paper_Plots/FigureS1/Supp_Figure1D.png", width=500, height=500)
plot(S2_pr, main = "S2 PRC")
dev.off()
