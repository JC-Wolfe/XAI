library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)

rules <- read.csv("~/Quicktest/Housekeeping_Rules.csv",
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
marks <- names(mcols(results)[1:26]) #Lazily getting a list of names

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

plotting_df$C1[plotting_df$C1=="high"] <- "forestgreen"
plotting_df$C1[plotting_df$C1=="medium"] <- "orange"
plotting_df$C1[plotting_df$C1=="low"] <- "red"

plotting_df$C2[plotting_df$C2=="high"] <- "forestgreen"
plotting_df$C2[plotting_df$C2=="medium"] <- "orange"
plotting_df$C2[plotting_df$C2=="low"] <- "red"

plotting_df$C3[plotting_df$C3=="high"] <- "forestgreen"
plotting_df$C3[plotting_df$C3=="medium"] <- "orange"
plotting_df$C3[plotting_df$C3=="low"] <- "red"

pos_df <- plotting_df[which(plotting_df$Result==1),]

score_buffer <- rep(0, length(marks))
for (i in seq_along(score_buffer)){

  M1_shared <- pos_df[which(pos_df$M1==marks[i]),]
  M1_score <- 1 * length(which(M1_shared$C1=="forestgreen"))+
              -1 * length(which(M1_shared$C1=="red"))

  M2_shared <- pos_df[which(pos_df$M2==marks[i]),]
  M2_score <- 1 * length(which(M2_shared$C2=="forestgreen"))+
              -1 * length(which(M2_shared$C2=="red"))

  M3_shared <- pos_df[which(pos_df$M3==marks[i]),]
  M3_score <- 1 * length(which(M3_shared$C3=="forestgreen"))+
              -1 * length(which(M3_shared$C3=="red"))

  score_buffer[i] <- M1_score + M2_score + M3_score
}

names(score_buffer) <- marks
score_buffer <- sort(score_buffer, decreasing = T)
marks_pos_scored <- names(score_buffer)

pos_df$M1 <- sapply(pos_df$M1, mark_number, marks = marks_pos_scored)
pos_df$M2 <- sapply(pos_df$M2, mark_number, marks = marks_pos_scored)
pos_df$M3 <- sapply(pos_df$M3, mark_number, marks = marks_pos_scored)

pos_df <- pos_df[order(pos_df$Dominance, decreasing = F),]


# Negative Marks Down Here For Now Because I'm Lazy






neg_df <- plotting_df[which(plotting_df$Result==0),]

neg_score_buffer <- rep(0, length(marks))
for (i in seq_along(neg_score_buffer)){

  M1_shared <- neg_df[which(neg_df$M1==marks[i]),]
  M1_score <- 1 * length(which(M1_shared$C1=="forestgreen"))+
              -1 * length(which(M1_shared$C1=="red"))

  M2_shared <- neg_df[which(neg_df$M2==marks[i]),]
  M2_score <- 1 * length(which(M2_shared$C2=="forestgreen"))+
              -1 * length(which(M2_shared$C2=="red"))

  M3_shared <- neg_df[which(neg_df$M3==marks[i]),]
  M3_score <- 1 * length(which(M3_shared$C3=="forestgreen"))+
              -1 * length(which(M3_shared$C3=="red"))

  neg_score_buffer[i] <- M1_score + M2_score + M3_score
}

names(neg_score_buffer) <- marks
neg_score_buffer <- sort(neg_score_buffer, decreasing = T)
marks_neg_scored <- names(neg_score_buffer)

neg_df$M1 <- sapply(neg_df$M1, mark_number, marks = marks_neg_scored)
neg_df$M2 <- sapply(neg_df$M2, mark_number, marks = marks_neg_scored)
neg_df$M3 <- sapply(neg_df$M3, mark_number, marks = marks_neg_scored)

neg_df <- neg_df[order(neg_df$Dominance, decreasing = F),]

# And plotting both

pdf("/home/jw18713/Project1/Paper_Plots/Figure4/Housekeeping_rule_circles.pdf",
  width = 20, height = 10, pointsize = 14)
par(mar=c(3,5.5,8.5,8), mfrow=c(1,2), cex = 1.2)
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
title(main = "Housekeeping Enhancer Rules", cex.main = 2, line = 5)

legend(27, length(pos_df[,1])*0.75, legend = c("High", "Medium", "Low"),
       fill = c("forestgreen", "orange", "red"),
       title = "Fuzzy Membership",
       bty="n")

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
title(main = "Housekeeping Non-Enhancer Rules", cex.main = 2, line = 5)

legend(27, length(neg_df[,1])*0.75, legend = c("High", "Medium", "Low"),
       fill = c("forestgreen", "orange", "red"),
       title = "Fuzzy Membership",
       bty="n")

dev.off()

rules <- read.csv("~/Quicktest/Developmental_Rules.csv",
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
marks <- names(mcols(results)[1:26]) #Lazily getting a list of names

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

plotting_df$C1[plotting_df$C1=="high"] <- "forestgreen"
plotting_df$C1[plotting_df$C1=="medium"] <- "orange"
plotting_df$C1[plotting_df$C1=="low"] <- "red"

plotting_df$C2[plotting_df$C2=="high"] <- "forestgreen"
plotting_df$C2[plotting_df$C2=="medium"] <- "orange"
plotting_df$C2[plotting_df$C2=="low"] <- "red"

plotting_df$C3[plotting_df$C3=="high"] <- "forestgreen"
plotting_df$C3[plotting_df$C3=="medium"] <- "orange"
plotting_df$C3[plotting_df$C3=="low"] <- "red"

pos_df <- plotting_df[which(plotting_df$Result==1),]

score_buffer <- rep(0, length(marks))
for (i in seq_along(score_buffer)){

  M1_shared <- pos_df[which(pos_df$M1==marks[i]),]
  M1_score <- 1 * length(which(M1_shared$C1=="forestgreen"))+
              -1 * length(which(M1_shared$C1=="red"))

  M2_shared <- pos_df[which(pos_df$M2==marks[i]),]
  M2_score <- 1 * length(which(M2_shared$C2=="forestgreen"))+
              -1 * length(which(M2_shared$C2=="red"))

  M3_shared <- pos_df[which(pos_df$M3==marks[i]),]
  M3_score <- 1 * length(which(M3_shared$C3=="forestgreen"))+
              -1 * length(which(M3_shared$C3=="red"))

  score_buffer[i] <- M1_score + M2_score + M3_score
}

names(score_buffer) <- marks
score_buffer <- sort(score_buffer, decreasing = T)
marks_pos_scored <- names(score_buffer)

pos_df$M1 <- sapply(pos_df$M1, mark_number, marks = marks_pos_scored)
pos_df$M2 <- sapply(pos_df$M2, mark_number, marks = marks_pos_scored)
pos_df$M3 <- sapply(pos_df$M3, mark_number, marks = marks_pos_scored)

pos_df <- pos_df[order(pos_df$Dominance, decreasing = F),]


# Negative Marks Down Here For Now Because I'm Lazy






neg_df <- plotting_df[which(plotting_df$Result==0),]

neg_score_buffer <- rep(0, length(marks))
for (i in seq_along(neg_score_buffer)){

  M1_shared <- neg_df[which(neg_df$M1==marks[i]),]
  M1_score <- 1 * length(which(M1_shared$C1=="forestgreen"))+
              -1 * length(which(M1_shared$C1=="red"))

  M2_shared <- neg_df[which(neg_df$M2==marks[i]),]
  M2_score <- 1 * length(which(M2_shared$C2=="forestgreen"))+
              -1 * length(which(M2_shared$C2=="red"))

  M3_shared <- neg_df[which(neg_df$M3==marks[i]),]
  M3_score <- 1 * length(which(M3_shared$C3=="forestgreen"))+
              -1 * length(which(M3_shared$C3=="red"))

  neg_score_buffer[i] <- M1_score + M2_score + M3_score
}

names(neg_score_buffer) <- marks
neg_score_buffer <- sort(neg_score_buffer, decreasing = T)
marks_neg_scored <- names(neg_score_buffer)

neg_df$M1 <- sapply(neg_df$M1, mark_number, marks = marks_neg_scored)
neg_df$M2 <- sapply(neg_df$M2, mark_number, marks = marks_neg_scored)
neg_df$M3 <- sapply(neg_df$M3, mark_number, marks = marks_neg_scored)

neg_df <- neg_df[order(neg_df$Dominance, decreasing = F),]

# And plotting both

pdf("/home/jw18713/Project1/Paper_Plots/Figure4/Developmental_rule_circles.pdf",
  width = 20, height = 10, pointsize = 14)
par(mar=c(3,5.5,8.5,8), mfrow=c(1,2), cex = 1.2)
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
title(main = "Developmental Enhancer Rules", cex.main = 2, line = 5)

legend(27, length(pos_df[,1])*0.75, legend = c("High", "Medium", "Low"),
       fill = c("forestgreen", "orange", "red"),
       title = "Fuzzy Membership",
       bty="n")

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
title(main = "Developmental Non-Enhancer Rules", cex.main = 2, line = 5)

legend(27, length(neg_df[,1])*0.75, legend = c("High", "Medium", "Low"),
       fill = c("forestgreen", "orange", "red"),
       title = "Fuzzy Membership",
       bty="n")

dev.off()




# Do this tomorrow with set diff
