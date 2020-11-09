
library(GenomicRanges)
library(rtracklayer)

new_rule_results <- read.csv("~/Project1/New_rules_csv.csv")

pdf(file = "~/Project1/Paper_Plots/Figure4/Figure4B.pdf",
width = 14, height = 20, pointsize = 14)
par(mfrow=c(3,1), cex = 1.2)
par(xpd = TRUE)

    bp1 <- barplot(new_rule_results$Average_Recall, ylim = c(-2, 2), space = 0.2,
        col = "#CC79A7", main = "Average Recall Change", yaxt = "none")
        axis(2, seq(-2, 2, 0.5), labels = paste0(seq(-2, 2, 0.5), "%"), las=2)
        lines(x=c(0,nrow(new_rule_results)+nrow(new_rule_results)*0.2),
        y=c(-2*1.1,-2*1.1), col="grey", lwd=0.5)
        text(bp1, -2*1.2, as.character(seq(1,nrow(new_rule_results))))
        text((nrow(new_rule_results)+nrow(new_rule_results)*0.2)/2, -2*1.4, "Rule ID")

    barplot(new_rule_results$Recall_1, ylim = c(-7.5, 7.5), space = 0.2,
        col = "#0072B2", main = "Enhancer Recall", yaxt = "none")
        axis(2, seq(-7.5, 7.5, 2.5), labels = paste0(seq(-7.5, 7.5, 2.5), "%"), las=2)
        lines(x=c(0,nrow(new_rule_results)+nrow(new_rule_results)*0.2),
        y=c(-7.5*1.1,-7.5*1.1), col="grey", lwd=0.5)
        text(bp1, -7.5*1.2, as.character(seq(1,nrow(new_rule_results))))
        text((nrow(new_rule_results)+nrow(new_rule_results)*0.2)/2, -7.5*1.4, "Rule ID")

    barplot(new_rule_results$Recall_0, ylim = c(-7.5, 7.5), space = 0.2,
        col = "#D55E00", main = "Non-Enhancer Recall", yaxt = "none")
        axis(2, seq(-7.5, 7.5, 2.5), labels = paste0(seq(-7.5, 7.5, 2.5), "%"), las=2)
        lines(x=c(0,nrow(new_rule_results)+nrow(new_rule_results)*0.2),
        y=c(-7.5*1.1,-7.5*1.1), col="grey", lwd=0.5)
        text(bp1, -7.5*1.2, as.character(seq(1,nrow(new_rule_results))))
        text((nrow(new_rule_results)+nrow(new_rule_results)*0.2)/2, -7.5*1.4, "Rule ID")

dev.off()
