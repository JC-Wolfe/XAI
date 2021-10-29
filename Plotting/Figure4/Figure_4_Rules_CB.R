library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

rules <- read.csv("~/Project1/BG3_FL_Rules.csv",
  stringsAsFactors = F)

rules_df <- data.frame(rules$Rule.Id,
            sapply(strsplit(rules$Antecedent.1, " "),"[[",1),
            sapply(strsplit(rules$Antecedent.1, " "),"[[",3))

names(rules_df) <- c("Id", "M1", "C1")

rules_df$M2 <- NA
rules_df$C2 <- NA

for (i in seq_along(rules$Antecedent.2)){
  splitrule <- strsplit(rules$Antecedent.2[i], " ")
  if (length(splitrule[[1]]) > 1){
    rules_df$M2[i] <- splitrule[[1]][1]
    rules_df$C2[i] <- splitrule[[1]][3]
  }
}

rules_df$M3 <- NA
rules_df$C3 <- NA

for (i in seq_along(rules$Antecedent.3)){
  splitrule <- strsplit(rules$Antecedent.3[i], " ")
  if (length(splitrule[[1]]) > 1){
    rules_df$M3[i] <- splitrule[[1]][1]
    rules_df$C3[i] <- splitrule[[1]][3]
  }
}

rules_df$Dominance <- rules$Dominance
rules_df$Result <- rules$Result

results <- get(load("/home/jw18713/Project1/data/model_predictions/BG3_FuzzyGR.Rda")) #This is specific to this model, this is bad
marks <- names(mcols(results)[1:26])

plotting_df <- rules_df

mark_number <- function(x, marks){
  x <- which(marks==x)
}

#plotting_df$M1 <- sapply(rules_df$M1, mark_number)
plotting_df$C1 <- as.character(rules_df$C1)
#plotting_df$M2 <- sapply(rules_df$M2, mark_number)
plotting_df$C2 <- as.character(rules_df$C2)
#plotting_df$M3 <- sapply(rules_df$M3, mark_number)
plotting_df$C3 <- as.character(rules_df$C3)
plotting_df$Dominance <- rules_df$Dominance
plotting_df$Result <- rules_df$Result

plotting_df$C1[plotting_df$C1=="high"] <- cbbPalette[4]
plotting_df$C1[plotting_df$C1=="medium"] <- cbbPalette[5]
plotting_df$C1[plotting_df$C1=="low"] <- cbbPalette[7]

plotting_df$C2[plotting_df$C2=="high"] <- cbbPalette[4]
plotting_df$C2[plotting_df$C2=="medium"] <- cbbPalette[5]
plotting_df$C2[plotting_df$C2=="low"] <- cbbPalette[7]

plotting_df$C3[plotting_df$C3=="high"] <- cbbPalette[4]
plotting_df$C3[plotting_df$C3=="medium"] <- cbbPalette[5]
plotting_df$C3[plotting_df$C3=="low"] <- cbbPalette[7]

pos_df <- plotting_df[which(plotting_df$Result==1),]

score_buffer <- rep(0, length(marks))
for (i in seq_along(score_buffer)){

  M1_shared <- pos_df[which(pos_df$M1==marks[i]),]
  M1_score <- 1 * length(which(M1_shared$C1==cbbPalette[4]))+
              -1 * length(which(M1_shared$C1==cbbPalette[7]))

  M2_shared <- pos_df[which(pos_df$M2==marks[i]),]
  M2_score <- 1 * length(which(M2_shared$C2==cbbPalette[4]))+
              -1 * length(which(M2_shared$C2==cbbPalette[7]))

  M3_shared <- pos_df[which(pos_df$M3==marks[i]),]
  M3_score <- 1 * length(which(M3_shared$C3==cbbPalette[4]))+
              -1 * length(which(M3_shared$C3==cbbPalette[7]))

  score_buffer[i] <- M1_score + M2_score + M3_score
}

names(score_buffer) <- marks
score_buffer <- sort(score_buffer, decreasing = T)
marks_pos_scored <- names(score_buffer)

pos_df$M1 <- sapply(pos_df$M1, mark_number, marks = marks_pos_scored)
pos_df$M2 <- sapply(pos_df$M2, mark_number, marks = marks_pos_scored)
pos_df$M3 <- sapply(pos_df$M3, mark_number, marks = marks_pos_scored)

pos_df <- pos_df[order(pos_df$Dominance, decreasing = F),]


neg_df <- plotting_df[which(plotting_df$Result==0),]

neg_score_buffer <- rep(0, length(marks))
for (i in seq_along(neg_score_buffer)){

  M1_shared <- neg_df[which(neg_df$M1==marks[i]),]
  M1_score <- 1 * length(which(M1_shared$C1==cbbPalette[4]))+
              -1 * length(which(M1_shared$C1==cbbPalette[7]))

  M2_shared <- neg_df[which(neg_df$M2==marks[i]),]
  M2_score <- 1 * length(which(M2_shared$C2==cbbPalette[4]))+
              -1 * length(which(M2_shared$C2==cbbPalette[7]))

  M3_shared <- neg_df[which(neg_df$M3==marks[i]),]
  M3_score <- 1 * length(which(M3_shared$C3==cbbPalette[4]))+
              -1 * length(which(M3_shared$C3==cbbPalette[7]))

  neg_score_buffer[i] <- M1_score + M2_score + M3_score
}

names(neg_score_buffer) <- marks
neg_score_buffer <- sort(neg_score_buffer, decreasing = T)
marks_neg_scored <- names(neg_score_buffer)

neg_df$M1 <- sapply(neg_df$M1, mark_number, marks = marks_neg_scored)
neg_df$M2 <- sapply(neg_df$M2, mark_number, marks = marks_neg_scored)
neg_df$M3 <- sapply(neg_df$M3, mark_number, marks = marks_neg_scored)

neg_df <- neg_df[order(neg_df$Dominance, decreasing = F),]
save(marks_neg_scored, file = "Negative_ordered_marks.Rda")

# And plotting both

pdf("/home/jw18713/Project1/Paper_Plots/Figure4/BG3_Enhancer_Rules.pdf",
  width = 10, height = 10, pointsize = 14)
par(mar=c(3,5.5,8.5,8), cex = 1.2)
par(xpd = NA)
plot(0,type="n",axes=F,xlab="",ylab="",xlim=c(1,26),ylim=c(1,length(pos_df$M1)))

for (i in seq(1, 26)){
  lines(x=c(i,i), y=c(1,length(pos_df[,1])), col="grey", lty=1, lwd=1)
  #abline(v = i, col="grey", lty=1, lwd=1)
}

for (i in seq(1, length(pos_df[,1]))){
  lines(x=c(1,26), y=c(i,i), col="grey", lty=2, lwd=0.5)
  #abline(h = i, col="grey", lty=2, lwd=0.5)
}

for (i in seq_along(pos_df$M1)){
  symbols(pos_df$M1[i], i, circles = 0.45, add = T, inches = F, bg = pos_df$C1[i])
}

for (i in seq_along(pos_df$M2)){
  symbols(pos_df$M2[i], i, circles = 0.45, add = T, inches = F, bg = pos_df$C2[i])
}

for (i in seq_along(pos_df$M3)){
  symbols(pos_df$M3[i], i, circles = 0.45, add = T, inches = F, bg = pos_df$C3[i])
}

axis(3, at=seq(1,26, by = 1), labels=F)
text(seq_along(marks_pos_scored),y = length(pos_df[,1]) + 1.75, srt = 315, adj = 1,labels = marks_pos_scored, xpd = TRUE,cex=0.75)
axis(2, at=seq(1,length(pos_df$Id), by = 1), labels=pos_df$Id, las=2)

title(ylab = "Rule ID", cex.lab = 1, line = 2)
title(main = "Enhancer Rules", cex.main = 2, line = 5)

legend(27, length(pos_df[,1])*0.75, legend = c("High", "Medium", "Low"),
       fill = c(cbbPalette[4], cbbPalette[5], cbbPalette[7]),
       bty="n")
dev.off()


pdf("/home/jw18713/Project1/Paper_Plots/Figure4/BG3_Non-Enhancer_Rules.pdf",
  width = 10, height = 10, pointsize = 14)
par(mar=c(3,5.5,8.5,8), cex = 1.2)
par(xpd = NA)
plot(0,type="n",axes=F,xlab="",ylab="",xlim=c(1,26),ylim=c(1,length(neg_df$M1)))

for (i in seq(1, 26)){
  lines(x=c(i,i), y=c(1,length(neg_df[,1])), col="grey", lty=1, lwd=1)
  #abline(v = i, col="grey", lty=1, lwd=1)
}

for (i in seq(1, length(neg_df[,1]))){
  lines(x=c(1,26), y=c(i,i), col="grey", lty=2, lwd=0.5)
  #abline(h = i, col="grey", lty=2, lwd=0.5)
}

for (i in seq_along(neg_df$M1)){
  symbols(neg_df$M1[i], i, circles = 0.45, add = T, inches = F, bg = neg_df$C1[i])
}

for (i in seq_along(neg_df$M2)){
  symbols(neg_df$M2[i], i, circles = 0.45, add = T, inches = F, bg = neg_df$C2[i])
}

for (i in seq_along(neg_df$M3)){
  symbols(neg_df$M3[i], i, circles = 0.45, add = T, inches = F, bg = neg_df$C3[i])
}

axis(3, at=seq(1,26, by = 1), labels=F)
text(seq_along(marks_neg_scored),y =length(neg_df[,1]) + 1.75, srt = 315, adj = 1,labels = marks_neg_scored, xpd = TRUE,cex=0.75)
axis(2, at=seq(1,length(neg_df$Id), by = 1), labels=neg_df$Id, las=2)

title(ylab = "Rule ID", cex.lab = 1, line = 2)
title(main = "Non-Enhancer Rules", cex.main = 2, line = 5)

legend(27, length(neg_df[,1])*0.75, legend = c("High", "Medium", "Low"),
       fill = c(cbbPalette[4], cbbPalette[5], cbbPalette[7]),
       bty="n")

dev.off()


# CUSTOM RULES

rules <- read.csv("~/Project1/New_rules_csv.csv", stringsAsFactors = F)

rules_df <- data.frame(rules[,1:8])

names(rules_df) <- c("ID", "M1", "L1", "M2", "L2", "M3", "L3", "Classification")



pos_plotting_df <- rules_df

pos_plotting_df$M1 <- sapply(rules_df$M1, mark_number, marks_pos_scored)
pos_plotting_df$L1 <- as.character(rules_df$L1)
pos_plotting_df$M2 <- sapply(rules_df$M2, mark_number, marks_pos_scored)
pos_plotting_df$L2 <- as.character(rules_df$L2)
pos_plotting_df$M3 <- sapply(rules_df$M3, mark_number, marks_pos_scored)
pos_plotting_df$L3 <- as.character(rules_df$L3)
pos_plotting_df$Classification <- rules_df$Classification

pos_plotting_df$L1[pos_plotting_df$L1=="high"] <- cbbPalette[4]
pos_plotting_df$L1[pos_plotting_df$L1=="medium"] <- cbbPalette[5]
pos_plotting_df$L1[pos_plotting_df$L1=="low"] <- cbbPalette[7]

pos_plotting_df$L2[pos_plotting_df$L2=="high"] <- cbbPalette[4]
pos_plotting_df$L2[pos_plotting_df$L2=="medium"] <- cbbPalette[5]
pos_plotting_df$L2[pos_plotting_df$L2=="low"] <- cbbPalette[7]

pos_plotting_df$L3[pos_plotting_df$L3=="high"] <- cbbPalette[4]
pos_plotting_df$L3[pos_plotting_df$L3=="medium"] <- cbbPalette[5]
pos_plotting_df$L3[pos_plotting_df$L3=="low"] <- cbbPalette[7]

pos_plotting_df$L2[pos_plotting_df$L2 == ""] <- NA
pos_plotting_df$L3[pos_plotting_df$L3 == ""] <- NA

pos_df <- pos_plotting_df[which(pos_plotting_df$Classification==1),]


neg_plotting_df <- rules_df

neg_plotting_df$M1 <- sapply(rules_df$M1, mark_number, marks_neg_scored)
neg_plotting_df$L1 <- as.character(rules_df$L1)
neg_plotting_df$M2 <- sapply(rules_df$M2, mark_number, marks_neg_scored)
neg_plotting_df$L2 <- as.character(rules_df$L2)
neg_plotting_df$M3 <- sapply(rules_df$M3, mark_number, marks_neg_scored)
neg_plotting_df$L3 <- as.character(rules_df$L3)
neg_plotting_df$Classification <- rules_df$Classification

neg_plotting_df$L1[neg_plotting_df$L1=="high"] <- cbbPalette[4]
neg_plotting_df$L1[neg_plotting_df$L1=="medium"] <- cbbPalette[5]
neg_plotting_df$L1[neg_plotting_df$L1=="low"] <- cbbPalette[7]

neg_plotting_df$L2[neg_plotting_df$L2=="high"] <- cbbPalette[4]
neg_plotting_df$L2[neg_plotting_df$L2=="medium"] <- cbbPalette[5]
neg_plotting_df$L2[neg_plotting_df$L2=="low"] <- cbbPalette[7]

neg_plotting_df$L3[neg_plotting_df$L3=="high"] <- cbbPalette[4]
neg_plotting_df$L3[neg_plotting_df$L3=="medium"] <- cbbPalette[5]
neg_plotting_df$L3[neg_plotting_df$L3=="low"] <- cbbPalette[7]

neg_plotting_df$L2[neg_plotting_df$L2 == ""] <- NA
neg_plotting_df$L3[neg_plotting_df$L3 == ""] <- NA

neg_df <- neg_plotting_df[which(neg_plotting_df$Classification==0),]


pdf("/home/jw18713/Project1/Paper_Plots/Figure4/Expert_Enhancer_Rules_Circles.pdf",
  width = 10, height = 8, pointsize = 14)
par(mar=c(3,5.5,8.5,8), cex = 1.2)
par(xpd = NA)
plot(0,type="n",axes=F,xlab="",ylab="",xlim=c(1,26),ylim=c(1,length(pos_df$M1)))

for (i in seq(1, 26)){
  lines(x=c(i,i), y=c(1,length(pos_df[,1])), col="grey", lty=1, lwd=1)
  #abline(v = i, col="grey", lty=1, lwd=1)
}

for (i in seq(1, length(pos_df[,1]))){
  lines(x=c(1,26), y=c(i,i), col="grey", lty=2, lwd=0.5)
  #abline(h = i, col="grey", lty=2, lwd=0.5)
}

for (i in seq_along(pos_df$M1)){
  symbols(pos_df$M1[i], i, circles = 0.45, add = T, inches = F, bg = pos_df$L1[i])
}

for (i in seq_along(pos_df$M2)){
  symbols(pos_df$M2[i], i, circles = 0.45, add = T, inches = F, bg = pos_df$L2[i])
}

for (i in seq_along(pos_df$M3)){
  symbols(pos_df$M3[i], i, circles = 0.45, add = T, inches = F, bg = pos_df$L3[i])
}

axis(3, at=seq(1,26, by = 1), labels=F)
text(seq_along(marks_pos_scored),y = length(pos_df[,1]) + 0.75, srt = 315,
  adj = 1,labels = marks_pos_scored, xpd = TRUE,cex=0.75)
axis(2, at=seq(1,length(pos_df$ID), by = 1), labels=pos_df$ID, las=2)

title(ylab = "Rule ID", cex.lab = 1, line = 2)
title(main = "Expert Enhancer Rules", cex.main = 2, line = 5)

legend(27, length(pos_df[,1])*0.75, legend = c("High", "Medium", "Low"),
       fill = c(cbbPalette[4], cbbPalette[5], cbbPalette[7]),
       bty="n")

dev.off()


pdf("/home/jw18713/Project1/Paper_Plots/Figure4/Expert_Non-Enhancer_Rules_Circles.pdf",
  width = 10, height = 6, pointsize = 14)
par(mar=c(3,5.5,8.5,8), cex = 1.2)
par(xpd = NA)

plot(0,type="n",axes=F,xlab="",ylab="",xlim=c(1,26),ylim=c(1,length(neg_df$M1)))

for (i in seq(1, 26)){
  lines(x=c(i,i), y=c(1,length(neg_df[,1])), col="grey", lty=1, lwd=1)
  #abline(v = i, col="grey", lty=1, lwd=1)
}

for (i in seq(1, length(neg_df[,1]))){
  lines(x=c(1,26), y=c(i,i), col="grey", lty=2, lwd=0.5)
  #abline(h = i, col="grey", lty=2, lwd=0.5)
}

for (i in seq_along(neg_df$M1)){
  symbols(neg_df$M1[i], i, circles = 0.45, add = T, inches = F, bg = neg_df$L1[i])
}

for (i in seq_along(neg_df$M2)){
  symbols(neg_df$M2[i], i, circles = 0.45, add = T, inches = F, bg = neg_df$L2[i])
}

for (i in seq_along(neg_df$M3)){
  symbols(neg_df$M3[i], i, circles = 0.45, add = T, inches = F, bg = neg_df$L3[i])
}

axis(3, at=seq(1,26, by = 1), labels=F, line = 0.5)
text(seq_along(marks_neg_scored),y =length(neg_df[,1]) + 0.9, srt = 315,
  adj = 1,labels = marks_neg_scored, xpd = TRUE,cex=0.75)
axis(2, at=seq(1,length(neg_df$ID), by = 1), labels=neg_df$ID, las=2)

title(ylab = "Rule ID", cex.lab = 1, line = 2)
title(main = "Expert Non-Enhancer Rules", cex.main = 2, line = 5)

legend(27, length(neg_df[,1])*0.9, legend = c("High", "Medium", "Low"),
       fill = c(cbbPalette[4], cbbPalette[5], cbbPalette[7]),
       bty="n")

dev.off()
