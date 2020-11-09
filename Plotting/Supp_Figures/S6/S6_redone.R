
library(seqLogo)
library(GenomicRanges)
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(PWMEnrich)
library(PWMEnrich.Dmelanogaster.background)
library(MotifDb)
library(gridExtra)
require(ggplot2)
require(ggseqlogo)
library(genomation)
library(rtracklayer)
library(VennDiagram)


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

readTFs <- function(x){
  names <- read.table(x, stringsAsFactors = F)[,1]
  return(names)
}

# Comparing within the same cell type
BG3_put <- readTFs("/home/jw18713/Project1/Paper_Plots/PWM/BG3/BG3_putative_report_pvalue05.txt")
BG3_common <- readTFs("/home/jw18713/Project1/Paper_Plots/PWM/BG3/BG3_common_report_pvalue05.txt")
S2_put <- readTFs("/home/jw18713/Project1/Paper_Plots/PWM/S2/S2_putative_report_pvalue05.txt")
S2_common <- readTFs("/home/jw18713/Project1/Paper_Plots/PWM/S2/S2_common_report_pvalue05.txt")

# Putative specific
BG3_put2 <- BG3_put[!BG3_put %in% BG3_common]
S2_put2 <- S2_put[!S2_put %in% S2_common]

# Common specific
BG3_common2 <- BG3_common[!BG3_common %in% BG3_put]
S2_common2 <- S2_common[!S2_common %in% S2_put]

generatePairwiseVenn <- function(set1, set2, categories, cols=c("#0072B2", "#D55E00"), cat.pos = c(-20, 20), cat.dist = c(0.07, 0.07), human = TRUE){

  v <- draw.pairwise.venn(area1 = length(set1), area2 = length(set2), cross.area = sum(set1%in%set2),
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



v1 <- generatePairwiseVenn(BG3_put, BG3_common, c("BG3 putative", "BG3 common"),
                           cols=c(cbbPalette[6], cbbPalette[7]),
                           cat.pos = c(-30, 150),
                           cat.dist = c(0.075, 0.075), human = F)


v2 <- generatePairwiseVenn(S2_put, S2_common, c("S2 putative", "S2 common"),
                           cols=c(cbbPalette[6], cbbPalette[7]),
                           cat.pos = c(-30, 150),
                           cat.dist = c(0.075, 0.075), human = F)


v3 <- generatePairwiseVenn(BG3_put2, S2_put2, c("BG3 putative specific", "S2 putative specific"),
                           cols=c(cbbPalette[6], cbbPalette[7]),
                           cat.pos = c(-30, 150),
                           cat.dist = c(0.075, 0.075), human = F)


v4 <- generatePairwiseVenn(BG3_common2, S2_common2, c("BG3 common specifc", "S2 common specific"),
                           cols=c(cbbPalette[6], cbbPalette[7]),
                           cat.pos = c(-30, 150),
                           cat.dist = c(0.075, 0.075), human = F)


pdf("/home/jw18713/Project1/Paper_Plots/FigureS6/FigureS6_venn.pdf", width=8, height=8,pointsize = 10);
par(mar=c(0, 0, 0, 0)+0.1)
pushViewport(plotViewport(layout=grid.layout(2, 2),gp=gpar(cex=1.0)))
pushViewport(plotViewport(layout.pos.row=1, layout.pos.col=1))
grid.draw(v1)
grid.text(bquote(bold(""~"BG3")), y=0.93, gp = gpar(cex=1.2))
grid.text(LETTERS[1], y=0.93, x=0.07, gp = gpar(cex=1.6))
popViewport()
pushViewport(plotViewport(layout.pos.row=1, layout.pos.col=2))
grid.draw(v2)
grid.text(bquote(bold(""~"S2")), y=0.93, gp = gpar(cex=1.2))
grid.text(LETTERS[2], y=0.93, x=0.07, gp = gpar(cex=1.6))
popViewport()
pushViewport(plotViewport(layout.pos.row=2, layout.pos.col=1))
grid.draw(v3)
grid.text(bquote(bold(""~"Putative Specific")), y=0.93, gp = gpar(cex=1.2))
grid.text(LETTERS[3], y=0.93, x=0.07, gp = gpar(cex=1.6))
popViewport()
pushViewport(plotViewport(layout.pos.row=2, layout.pos.col=2))
grid.draw(v4)
grid.text(bquote(bold(""~"Common Specific")), y=0.93, gp = gpar(cex=1.2))
grid.text(LETTERS[4], y=0.93, x=0.07, gp = gpar(cex=1.6))
popViewport()
dev.off()

put_names <- BG3_put2[BG3_put2 %in% S2_put2]
common_names <- BG3_common2[BG3_common2 %in% S2_common2]

# some of the names are not the official symbols
# put_motifs[2] is dl_2, correct symbol for plot is dl
# put_motifs[3] is CG14962, correct symbol for plot is Asciz
# common_motifs[1] CNC::maf-S is CNC with maf-s, id wise I should probably leave it this way
# I'll use the original target names to get the correct motifs though

put_plot_names <- put_names
put_plot_names[2] <- "dl"
put_plot_names[3] <- "Asciz"

common_plot_names <- common_names

put_report <- get(load("/home/jw18713/BG3_put_full_report_0.5.Rda"))
common_report <- get(load("/home/jw18713/BG3_common_full_report_0.5.Rda"))

report_id_extract <- function(report, names){
  ids <- rep(0, length(names))
  for (i in seq_along(names)){
    ids[i] <- report$id[report$target == names[i]][1]
  }
  return(ids)
}

put_ids <- report_id_extract(put_report, put_names)
common_ids <- report_id_extract(common_report, common_names)

pwm_grabber <- function(x){
  pwms <- list()
  for (i in x){
    pwm <- as.list(query(MotifDb,i)[1])
    pwms <- c(pwm, pwms)
  }
  return(pwms)
}

put_pwms <- pwm_grabber(put_ids)
names(put_pwms) <- put_plot_names
common_pwms <- pwm_grabber(common_ids)
names(common_pwms) <- common_plot_names

pdf("/home/jw18713/Project1/Paper_Plots/FigureS6/put_pwms.pdf",
  width = 14, height = 6, pointsize = 14)
par(cex = 1.8)
ggseqlogo(put_pwms, ncol=3)
dev.off()

pdf("/home/jw18713/Project1/Paper_Plots/FigureS6/common_pwms.pdf",
width = 14, height = 6, pointsize = 14)
par(cex = 1.8)
ggseqlogo(common_pwms, ncol=3)
dev.off()
