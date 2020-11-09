load("BG3_matrix.Rda")
load("S2_matrix.Rda")

BG3_matrix_neither <- matrix(0, ncol=2, nrow=4)
rownames(BG3_matrix_neither) <- c("Putative_Neither", "Common_Neither", "STARR-seq_Neither", "Background_Neither")
BG3_matrix_neither[1,2] <- 279
BG3_matrix_neither[2,2] <- 11
BG3_matrix_neither[3,2] <- 49
BG3_matrix_neither[4,2] <- 29652

S2_matrix_neither <- matrix(0, ncol=2, nrow=4)
rownames(S2_matrix_neither) <- c("Putative_Neither", "Common_Neither", "STARR-seq_Neither", "Background_Neither")
S2_matrix_neither[1,2] <- 248
S2_matrix_neither[2,2] <- 25
S2_matrix_neither[3,2] <- 31
S2_matrix_neither[4,2] <- 27461

BG3_matrix_all <- rbind(BG3_matrix, BG3_matrix_neither)
BG3_matrix_all <- BG3_matrix_all[c(1,2,9,3,4,10,5,6,11,7,8,12),]


S2_matrix_all <- rbind(S2_matrix, S2_matrix_neither)
S2_matrix_all <- S2_matrix_all[c(1,2,9,3,4,10,5,6,11,7,8,12),]

.fisherTest <-  function(x, alternative = c("two.sided", "less", "greater")){
  return(fisher.test(matrix(unlist(x),nrow=2, byrow=TRUE), alternative = alternative))
}


BG3_putative_pvalue <- .fisherTest(as.vector(t(BG3_matrix[1:2,1:2])), alternative ="two.sided")
BG3_common_pvalue <- .fisherTest(as.vector(t(BG3_matrix[3:4,1:2])), alternative ="two.sided")
BG3_starr_seq_pvalue <- .fisherTest(as.vector(t(BG3_matrix[5:6,1:2])), alternative ="two.sided")
BG3_background_seq_pvalue <- .fisherTest(as.vector(t(BG3_matrix[7:8,1:2])), alternative ="two.sided")

BG3_proximal_putative_common <- .fisherTest(as.vector(t(BG3_matrix[c(1,3),1:2])), alternative ="two.sided")
BG3_proximal_putative_starr_seq <- .fisherTest(as.vector(t(BG3_matrix[c(1,5),1:2])), alternative ="two.sided")
BG3_proximal_putative_background <- .fisherTest(as.vector(t(BG3_matrix[c(1,7),1:2])), alternative ="two.sided")


BG3_distal_putative_common <- .fisherTest(as.vector(t(BG3_matrix[c(2,4),1:2])), alternative ="two.sided")
BG3_distal_putative_starr_seq <- .fisherTest(as.vector(t(BG3_matrix[c(2,6),1:2])), alternative ="two.sided")
BG3_distal_putative_background <- .fisherTest(as.vector(t(BG3_matrix[c(2,7),1:2])), alternative ="two.sided")



S2_putative_pvalue <- .fisherTest(as.vector(t(S2_matrix[1:2,1:2])), alternative ="two.sided")
S2_common_pvalue <- .fisherTest(as.vector(t(S2_matrix[3:4,1:2])), alternative ="two.sided")
S2_starr_seq_pvalue <- .fisherTest(as.vector(t(S2_matrix[5:6,1:2])), alternative ="two.sided")
S2_background_seq_pvalue <- .fisherTest(as.vector(t(S2_matrix[7:8,1:2])), alternative ="two.sided")

S2_proximal_putative_common <- .fisherTest(as.vector(t(S2_matrix[c(1,3),1:2])), alternative ="two.sided")
S2_proximal_putative_starr_seq <- .fisherTest(as.vector(t(S2_matrix[c(1,5),1:2])), alternative ="two.sided")
S2_proximal_putative_background <- .fisherTest(as.vector(t(S2_matrix[c(1,7),1:2])), alternative ="two.sided")


S2_distal_putative_common <- .fisherTest(as.vector(t(S2_matrix[c(2,4),1:2])), alternative ="two.sided")
S2_distal_putative_starr_seq <- .fisherTest(as.vector(t(S2_matrix[c(2,6),1:2])), alternative ="two.sided")
S2_distal_putative_background <- .fisherTest(as.vector(t(S2_matrix[c(2,7),1:2])), alternative ="two.sided")



cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")



BG3_matrix_total <- matrix(0, ncol=2, nrow=4)

BG3_matrix_total[1,] <- apply(BG3_matrix_all[1:3,],2,sum)
BG3_matrix_total[2,] <- apply(BG3_matrix_all[4:6,],2,sum)
BG3_matrix_total[3,] <- apply(BG3_matrix_all[7:9,],2,sum)
BG3_matrix_total[4,] <- apply(BG3_matrix_all[10:12,],2,sum)
rownames(BG3_matrix_total) <- unlist(strsplit(rownames(BG3_matrix_all)[c(T,F,F)], "_"))[c(T,F)]
colnames(BG3_matrix_total) <- c("Expressed", "Not_Expressed")


BG3_matrix_total_putative_common <- .fisherTest(as.vector(t(BG3_matrix_total[c(1,2),1:2])), alternative ="two.sided")
BG3_matrix_total_putative_starr_seq <- .fisherTest(as.vector(t(BG3_matrix_total[c(1,3),1:2])), alternative ="two.sided")
BG3_matrix_total_putative_background <- .fisherTest(as.vector(t(BG3_matrix_total[c(1,4),1:2])), alternative ="two.sided")



S2_matrix_total_putative_common <- .fisherTest(as.vector(t(S2_matrix_total[c(1,2),1:2])), alternative ="two.sided")
S2_matrix_total_putative_starr_seq <- .fisherTest(as.vector(t(S2_matrix_total[c(1,3),1:2])), alternative ="two.sided")
S2_matrix_total_putative_background <- .fisherTest(as.vector(t(S2_matrix_total[c(1,4),1:2])), alternative ="two.sided")


BG3_matrix_total_prop <- matrix(BG3_matrix_total[,1]/(BG3_matrix_total[,1] + BG3_matrix_total[,2]), 
                              ncol=1, byrow = T)
colnames(BG3_matrix_total_prop) <- c("All")
rownames(BG3_matrix_total_prop) <- rownames(BG3_matrix_total)
BG3_matrix_total_prop <- BG3_matrix_total_prop[4:1,]
BG3_matrix_total_prop <- 100*BG3_matrix_total_prop


S2_matrix_total <- matrix(0, ncol=2, nrow=4)

S2_matrix_total[1,] <- apply(S2_matrix_all[1:3,],2,sum)
S2_matrix_total[2,] <- apply(S2_matrix_all[4:6,],2,sum)
S2_matrix_total[3,] <- apply(S2_matrix_all[7:9,],2,sum)
S2_matrix_total[4,] <- apply(S2_matrix_all[10:12,],2,sum)
rownames(S2_matrix_total) <- unlist(strsplit(rownames(S2_matrix_all)[c(T,F,F)], "_"))[c(T,F)]
colnames(S2_matrix_total) <- c("Expressed", "Not_Expressed")


S2_matrix_total_prop <- matrix(S2_matrix_total[,1]/(S2_matrix_total[,1] + S2_matrix_total[,2]), 
                               ncol=1, byrow = T)
colnames(S2_matrix_total_prop) <- c("All")
rownames(S2_matrix_total_prop) <- rownames(S2_matrix_total)
S2_matrix_total_prop <- S2_matrix_total_prop[4:1,]
S2_matrix_total_prop <- 100*S2_matrix_total_prop


BG3_matrix_prop <- matrix(BG3_matrix[,1]/(BG3_matrix[,1] + BG3_matrix[,2]), 
                          ncol=2, byrow = T)

colnames(BG3_matrix_prop) <- c("Proximal", "Distal Only")
rownames(BG3_matrix_prop) <- unlist(strsplit(rownames(BG3_matrix)[c(T,F)], "_"))[c(T,F)]
BG3_matrix_prop <- BG3_matrix_prop[4:1,]
BG3_matrix_prop <- 100*(BG3_matrix_prop)
BG3_matrix_all_prop <- cbind(BG3_matrix_prop, BG3_matrix_total_prop)
colnames(BG3_matrix_all_prop) <- c("Proximal", "Distal Only", "All")

S2_matrix_prop <- matrix(S2_matrix[,1]/(S2_matrix[,1] + S2_matrix[,2]), 
                         ncol=2, byrow = T)

colnames(S2_matrix_prop) <- c("Proximal", "Distal Only")
rownames(S2_matrix_prop) <- unlist(strsplit(rownames(S2_matrix)[c(T,F)], "_"))[c(T,F)]
S2_matrix_prop <- S2_matrix_prop[4:1,]
S2_matrix_prop <- 100*(S2_matrix_prop)
S2_matrix_all_prop <- cbind(S2_matrix_prop, S2_matrix_total_prop)
colnames(S2_matrix_all_prop) <- c("Proximal", "Distal Only", "All")



cols <- c("white", cbPalette[c(5,2,7)])


pdf(file = "Figure3_Horiz_BG_split.pdf", width = 21, height = 6, pointsize = 14)
par(mar=c(4,5.5,4,15), mfrow=c(1,2), cex = 1.2)
par(xpd = NA)

  barplot(BG3_matrix_prop, beside=T, horiz = T, col = cols,
      xlim = c(0,100), xlab = "% of enhancers that contact expressed genes", xaxt = "n", 
      main = "Enhancers in BG3 contact expressed genes")
    axis(1, seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"), las=2)
    legend(110,3.5,rev(rownames(BG3_matrix_prop)), fill=rev(cols), bty='n')

    
  barplot(S2_matrix_prop, beside=T, horiz = T, col = cols,
          xlim = c(0,100), xlab = "% of enhancers that contact expressed genes", xaxt = "n", 
          main = "Enhancers in S2 contact expressed genes")
  axis(1, seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"), las=2)
  legend(110,3.5,rev(rownames(S2_matrix_prop)), fill=rev(cols), bty='n')
  
    
dev.off()



pdf(file = "Figure3_Horiz_BG_all.pdf", width = 21, height = 6, pointsize = 14)
par(mar=c(4,5.5,4,15), mfrow=c(1,2), cex = 1.2)
par(xpd = NA)

barplot(BG3_matrix_all_prop, beside=T, horiz = T, col = cols,
        xlim = c(0,100), xlab = "% of enhancers that contact expressed genes", xaxt = "n", 
        main = "Enhancers in BG3 contact expressed genes")
axis(1, seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"), las=2)
legend(110,3.5,rev(rownames(BG3_matrix_prop)), fill=rev(cols), bty='n')


barplot(S2_matrix_all_prop, beside=T, horiz = T, col = cols,
        xlim = c(0,100), xlab = "% of enhancers that contact expressed genes", xaxt = "n", 
        main = "Enhancers in S2 contact expressed genes")
axis(1, seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"), las=2)
legend(110,3.5,rev(rownames(S2_matrix_prop)), fill=rev(cols), bty='n')


dev.off()


