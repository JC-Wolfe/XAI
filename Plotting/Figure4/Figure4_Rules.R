library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)

rules <- read.csv("~/Project1/New_rules_csv.csv", stringsAsFactors = F)

rules_df <- data.frame(rules[,1:8])

names(rules_df) <- c("ID", "M1", "L1", "M2", "L2", "M3", "L3", "Classification")

load("Positive_ordered_marks.Rda") # marks_pos_scored
load("Negative_ordered_marks.Rda") # marks_neg_scored

pos_df <- rules_df[which(rules_df$Classification==1),]
pos_plot <- pos_df
neg_df <- rules_df[which(rules_df$Classification==0),]
neg_plot <- neg_df

mark_number <- function(x, marks){
  x <- which(marks==x)
}

# Positive rules

pos_plot$M1 <- sapply(pos_df$M1, mark_number, marks_pos_scored)
pos_plot$L1 <- as.character(pos_df$L1)
pos_plot$M2 <- sapply(pos_df$M2, mark_number, marks_pos_scored)
pos_plot$L2 <- as.character(pos_df$L2)
pos_plot$M3 <- sapply(pos_df$M3, mark_number, marks_pos_scored)
pos_plot$L3 <- as.character(pos_df$L3)
pos_plot$Classification <- pos_df$Classification

pos_plot$L1[pos_plot$L1=="high"] <- "forestgreen"
pos_plot$L1[pos_plot$L1=="medium"] <- "orange"
pos_plot$L1[pos_plot$L1=="low"] <- "red"

pos_plot$L2[pos_plot$L2=="high"] <- "forestgreen"
pos_plot$L2[pos_plot$L2=="medium"] <- "orange"
pos_plot$L2[pos_plot$L2=="low"] <- "red"

pos_plot$L3[pos_plot$L3=="high"] <- "forestgreen"
pos_plot$L3[pos_plot$L3=="medium"] <- "orange"
pos_plot$L3[pos_plot$L3=="low"] <- "red"

pos_plot$L2[pos_plot$L2 == ""] <- NA
pos_plot$L3[pos_plot$L3 == ""] <- NA

# Negative rules

neg_plot$M1 <- sapply(neg_df$M1, mark_number, marks_neg_scored)
neg_plot$L1 <- as.character(neg_df$L1)
neg_plot$M2 <- sapply(neg_df$M2, mark_number, marks_neg_scored)
neg_plot$L2 <- as.character(neg_df$L2)
neg_plot$M3 <- sapply(neg_df$M3, mark_number, marks_neg_scored)
neg_plot$L3 <- as.character(neg_df$L3)
neg_plot$Classification <- neg_df$Classification

neg_plot$L1[neg_plot$L1=="high"] <- "forestgreen"
neg_plot$L1[neg_plot$L1=="medium"] <- "orange"
neg_plot$L1[neg_plot$L1=="low"] <- "red"

neg_plot$L2[neg_plot$L2=="high"] <- "forestgreen"
neg_plot$L2[neg_plot$L2=="medium"] <- "orange"
neg_plot$L2[neg_plot$L2=="low"] <- "red"

neg_plot$L3[neg_plot$L3=="high"] <- "forestgreen"
neg_plot$L3[neg_plot$L3=="medium"] <- "orange"
neg_plot$L3[neg_plot$L3=="low"] <- "red"

neg_plot$L2[neg_plot$L2 == ""] <- NA
neg_plot$L3[neg_plot$L3 == ""] <- NA




# All of this is plotting stuff


pdf("/home/jw18713/Project1/Paper_Plots/Figure4/Additional_Rules_Circles.pdf",
  width = 10, height = 15, pointsize = 14)
par(mar=c(3,5.5,8.5,8), mfrow=c(2,1), cex = 1.2)
par(xpd = NA)
plot(0,type="n",axes=F,xlab="",ylab="",xlim=c(1,26),ylim=c(1,length(pos_plot$M1)))

for (i in seq(1, 26)){
  lines(x=c(i,i), y=c(1,length(pos_plot[,1])), col="grey", lty=1, lwd=1)
  #abline(v = i, col="grey", lty=1, lwd=1)
}

for (i in seq(1, length(pos_plot[,1]))){
  lines(x=c(1,26), y=c(i,i), col="grey", lty=2, lwd=0.5)
  #abline(h = i, col="grey", lty=2, lwd=0.5)
}

for (i in seq_along(pos_plot$M1)){
  symbols(pos_plot$M1[i], i, circles = 0.45, add = T, inches = F, bg = pos_plot$L1[i])
}

for (i in seq_along(pos_plot$M2)){
  symbols(pos_plot$M2[i], i, circles = 0.45, add = T, inches = F, bg = pos_plot$L2[i])
}

for (i in seq_along(pos_plot$M3)){
  symbols(pos_plot$M3[i], i, circles = 0.45, add = T, inches = F, bg = pos_plot$L3[i])
}

axis(3, at=seq(1,26, by = 1), labels=F)
text(seq_along(marks_pos_scored),y = length(pos_plot[,1]) + 1, srt = 315, adj = 1,labels = marks_pos_scored, xpd = TRUE,cex=0.75)
axis(2, at=seq(1,length(pos_plot$ID), by = 1), labels=pos_plot$ID, las=2)

title(ylab = "Rule ID", cex.lab = 1, line = 2)
title(main = "BG3 Enhancer Rules", cex.main = 2, line = 5)

legend(27, length(pos_plot[,1])*0.75, legend = c("High", "Medium", "Low"),
       fill = c("forestgreen", "orange", "red"),
       title = "Fuzzy Membership",
       bty="n")

plot(0,type="n",axes=F,xlab="",ylab="",xlim=c(1,26),ylim=c(1,length(neg_plot$M1)))

for (i in seq(1, 26)){
  lines(x=c(i,i), y=c(1,length(neg_plot[,1])), col="grey", lty=1, lwd=1)
  #abline(v = i, col="grey", lty=1, lwd=1)
}

for (i in seq(1, length(neg_plot[,1]))){
  lines(x=c(1,26), y=c(i,i), col="grey", lty=2, lwd=0.5)
  #abline(h = i, col="grey", lty=2, lwd=0.5)
}

for (i in seq_along(neg_plot$M1)){
  symbols(neg_plot$M1[i], i, circles = 0.45, add = T, inches = F, bg = neg_plot$L1[i])
}

for (i in seq_along(neg_plot$M2)){
  symbols(neg_plot$M2[i], i, circles = 0.45, add = T, inches = F, bg = neg_plot$L2[i])
}

for (i in seq_along(neg_plot$M3)){
  symbols(neg_plot$M3[i], i, circles = 0.45, add = T, inches = F, bg = neg_plot$L3[i])
}

axis(3, at=seq(1,26, by = 1), labels=F)
text(seq_along(marks_neg_scored),y =length(neg_plot[,1]) + 0.55, srt = 315, adj = 1,labels = marks_neg_scored, xpd = TRUE,cex=0.75)
axis(2, at=seq(1,length(neg_plot$ID), by = 1), labels=neg_plot$ID, las=2)

title(ylab = "Rule ID", cex.lab = 1, line = 2)
title(main = "BG3 Non-Enhancer Rules", cex.main = 2, line = 5)

legend(27, length(neg_plot[,1])*0.75, legend = c("High", "Medium", "Low"),
       fill = c("forestgreen", "orange", "red"),
       title = "Fuzzy Membership",
       bty="n")

dev.off()
