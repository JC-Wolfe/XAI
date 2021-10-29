library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)
library(gplots)
library(lattice)

BG3_pred_enh <- get(load("~/dm6_promoter_correction/predicted_enhancers/grown_predicted_enhancers_BG3.Rda")) # BG3_pred_enh
BG3_gr_contacts <- get(load("~/dm6_promoter_correction/grange_contact_objects/useful_HiC_dm6_contactMap_BG3.Rda"))
BG3_gr_contacts$contact_chr <- paste0("chr", as.character(BG3_gr_contacts$contact_chr))
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
BG3_starr <- get(load("~/dm6_promoter_correction/starr_seq/starr.Rda")) # as BG3_starr
BG3_400 <- get(load("~/dm6_promoter_correction/Rda_Objects/empty_tiled_400.Rda"))
ann_genome <- get(load(file = "~/dm6_promoter_correction/dm6_annotations/dm6_annotation.Rda"))
proms <- reduce(ann_genome[ann_genome$annotations == "promoter"])
t <- 0

BG3_shared_enhancers <- BG3_pred_enh[BG3_pred_enh %over% BG3_starr]
BG3_putative_enhancers <- BG3_pred_enh[!BG3_pred_enh %over% BG3_starr]
BG3_starr_only <- BG3_starr[!BG3_starr %over% BG3_pred_enh]

contact_mapping <- function(x, contacts, promoters, distance = 5000, threshold = 0){

  over5k <- contacts[(end(contacts) + distance < contacts$contact_start
                    |!as.character(seqnames(contacts)) == contacts$contact_chr)
                    & contacts$enrichment > threshold
                    & contacts$promoter_contact == TRUE]

  nearest_promoter <- nearest(x, promoters)

  prox <- x[distance(x,promoters[nearest_promoter]) <= distance]


  # Distal overlaps
  dist <- x[x %over% over5k]
  # Proximal
  # Distal and Proximal
  d_and_p <- dist[dist %over% prox]
  # Distal only
  d_only <- dist[!dist %over% prox]
  # Proximal only
  p_only <- prox[!prox %over% dist]
  # Neither
  neither <- x[!x %over% dist & !x %over% prox]

  # GRange List Output
  out <- GRangesList(d_and_p, d_only, p_only, neither)
  names(out) <- c("Distal_&_Proximal_Overlaps",
                  "Distal_Only_Overlaps",
                  "Proximal_Only_Overlaps",
                  "No_Promoter_Overlaps")
  return(out)
}

BG3_putative <- contact_mapping(BG3_putative_enhancers, BG3_gr_contacts, proms, threshold = t)
BG3_shared <- contact_mapping(BG3_shared_enhancers, BG3_gr_contacts, proms, threshold = t)
BG3_starr_res <- contact_mapping(BG3_starr_only, BG3_gr_contacts, proms, threshold = t)
BG3_bg <- contact_mapping(BG3_400, BG3_gr_contacts, proms, threshold = t)
save(BG3_putative, file = "~/dm6_promoter_correction/Class_split_glists/BG3_putative_list.Rda")
save(BG3_shared, file = "~/dm6_promoter_correction/Class_split_glists/BG3_shared_list.Rda")
save(BG3_starr_res, file = "~/dm6_promoter_correction/Class_split_glists/BG3_starr_res_list.Rda")
save(BG3_bg, file = "~/dm6_promoter_correction/Class_split_glists/BG3_bg_list.Rda")

BG3_mat <- matrix(0,4,4)
rownames(BG3_mat) <- c("Putative", "Common", "STARR-seq", "Background")
colnames(BG3_mat) <- c("Distal & Proximal", "Distal Only", "Proximal Only", "Neither")

BG3_mat[1,1] <- length(BG3_shared[[1]])
BG3_mat[1,2] <- length(BG3_shared[[2]])
BG3_mat[1,3] <- length(BG3_shared[[3]])
BG3_mat[1,4] <- length(BG3_shared[[4]])

BG3_mat[2,1] <- length(BG3_putative[[1]])
BG3_mat[2,2] <- length(BG3_putative[[2]])
BG3_mat[2,3] <- length(BG3_putative[[3]])
BG3_mat[2,4] <- length(BG3_putative[[4]])

BG3_mat[3,1] <- length(BG3_starr_res[[1]])
BG3_mat[3,2] <- length(BG3_starr_res[[2]])
BG3_mat[3,3] <- length(BG3_starr_res[[3]])
BG3_mat[3,4] <- length(BG3_starr_res[[4]])

BG3_mat[4,1] <- length(BG3_bg[[1]])
BG3_mat[4,2] <- length(BG3_bg[[2]])
BG3_mat[4,3] <- length(BG3_bg[[3]])
BG3_mat[4,4] <- length(BG3_bg[[4]])

BG3_bar <- t(BG3_mat)
BG3_bar_aves <- BG3_bar

BG3_bar_aves[,1] <- BG3_bar[,1]/sum(BG3_bar[,1])
BG3_bar_aves[,2] <- BG3_bar[,2]/sum(BG3_bar[,2])
BG3_bar_aves[,3] <- BG3_bar[,3]/sum(BG3_bar[,3])
BG3_bar_aves[,4] <- BG3_bar[,4]/sum(BG3_bar[,4])


S2_pred_enh <- get(load("~/dm6_promoter_correction/predicted_enhancers/grown_predicted_enhancers_S2.Rda")) # S2_pred_enh
S2_gr_contacts <- get(load("~/dm6_promoter_correction/grange_contact_objects/useful_HiC_dm6_contactMap_S2.Rda"))
S2_gr_contacts$contact_chr <- paste0("chr", as.character(S2_gr_contacts$contact_chr))
S2_starr <- get(load("/home/jw18713/dm6_promoter_correction/starr_seq/starr_S2.Rda")) # as S2_starr
S2_starr <- S2_starr[seqnames(S2_starr) %in% seqlevels(S2_starr)[1:5]]
S2_400 <- get(load("~/dm6_promoter_correction/Rda_Objects/empty_tiled_400.Rda"))


S2_shared_enhancers <- S2_pred_enh[S2_pred_enh %over% S2_starr]
S2_putative_enhancers <- S2_pred_enh[!S2_pred_enh %over% S2_starr]
S2_starr_only <- S2_starr[!S2_starr %over% S2_pred_enh]



S2_putative <- contact_mapping(S2_putative_enhancers, S2_gr_contacts, proms, threshold = t)
S2_shared <- contact_mapping(S2_shared_enhancers, S2_gr_contacts, proms, threshold = t)
S2_starr_res <- contact_mapping(S2_starr_only, S2_gr_contacts, proms, threshold = t)
S2_bg <- contact_mapping(S2_400, S2_gr_contacts, proms, threshold = t)

save(S2_putative, file = "~/dm6_promoter_correction/Class_split_glists/S2_putative_list.Rda")
save(S2_shared, file = "~/dm6_promoter_correction/Class_split_glists/S2_shared_list.Rda")
save(S2_starr_res, file = "~/dm6_promoter_correction/Class_split_glists/S2_starr_res_list.Rda")
save(S2_bg, file = "~/dm6_promoter_correction/Class_split_glists/S2_bg_list.Rda")

S2_mat <- matrix(0,4,4)
rownames(S2_mat) <- c("Putative", "Common", "STARR-seq", "Background")
colnames(S2_mat) <- c("Distal & Proximal", "Distal Only", "Proximal Only", "Neither")

S2_mat[1,1] <- length(S2_shared[[1]])
S2_mat[1,2] <- length(S2_shared[[2]])
S2_mat[1,3] <- length(S2_shared[[3]])
S2_mat[1,4] <- length(S2_shared[[4]])

S2_mat[2,1] <- length(S2_putative[[1]])
S2_mat[2,2] <- length(S2_putative[[2]])
S2_mat[2,3] <- length(S2_putative[[3]])
S2_mat[2,4] <- length(S2_putative[[4]])

S2_mat[3,1] <- length(S2_starr_res[[1]])
S2_mat[3,2] <- length(S2_starr_res[[2]])
S2_mat[3,3] <- length(S2_starr_res[[3]])
S2_mat[3,4] <- length(S2_starr_res[[4]])

S2_mat[4,1] <- length(S2_bg[[1]])
S2_mat[4,2] <- length(S2_bg[[2]])
S2_mat[4,3] <- length(S2_bg[[3]])
S2_mat[4,4] <- length(S2_bg[[4]])

S2_bar <- t(S2_mat)
S2_bar_aves <- S2_bar

S2_bar_aves[,1] <- S2_bar[,1]/sum(S2_bar[,1])
S2_bar_aves[,2] <- S2_bar[,2]/sum(S2_bar[,2])
S2_bar_aves[,3] <- S2_bar[,3]/sum(S2_bar[,3])
S2_bar_aves[,4] <- S2_bar[,4]/sum(S2_bar[,4])

cs <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

pdf(paste0("~/Paper_replots/Plots/Figure3A.pdf"),
width = 22, height = 8, pointsize = 14)
par(mar=c(4,4,4,10.5), mfrow = c(1,2), cex = 1.2)
par(xpd = NA)

barplot(BG3_bar_aves, col=cs, xlab="Result Classification", ylab="Percentage of Results in Region", main="BG3 Annotation Comparisons", yaxt = "none")
axis(2, seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10), "%"), las=2)
legend(5.25,0.75,rownames(BG3_bar_aves), fill=cs, bty = "n")

barplot(S2_bar_aves, col=cs, xlab="Result Classification", ylab="Percentage of Results in Region", main="S2 Annotation Comparisons", yaxt = "none")
axis(2, seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10), "%"), las=2)
legend(5.25,0.75,rownames(S2_bar_aves), fill=cs, bty = "n")

dev.off()

BG3_hmap <- BG3_bar
BG3_hmap_aves <- BG3_bar_aves


BG3_hmap_comp <- BG3_hmap_aves[,1:3]
BG3_hmap_comp[1,] <- log2(BG3_hmap_comp[1,]/BG3_hmap_aves[1,4])
BG3_hmap_comp[2,] <- log2(BG3_hmap_comp[2,]/BG3_hmap_aves[2,4])
BG3_hmap_comp[3,] <- log2(BG3_hmap_comp[3,]/BG3_hmap_aves[3,4])
BG3_hmap_comp[4,] <- log2(BG3_hmap_comp[4,]/BG3_hmap_aves[4,4])

S2_hmap <- S2_bar
S2_hmap_aves <- S2_bar_aves


S2_hmap_comp <- S2_hmap_aves[,1:3]
S2_hmap_comp[1,] <- log2(S2_hmap_comp[1,]/S2_hmap_aves[1,4])
S2_hmap_comp[2,] <- log2(S2_hmap_comp[2,]/S2_hmap_aves[2,4])
S2_hmap_comp[3,] <- log2(S2_hmap_comp[3,]/S2_hmap_aves[3,4])
S2_hmap_comp[4,] <- log2(S2_hmap_comp[4,]/S2_hmap_aves[4,4])

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
custom_at <- seq(-5, 5, by = (10/30))
cols_contrast <- rev(c(colorRampPalette(c(cbbPalette[7], "white", cbbPalette[6]))(30)))
library(gridExtra)

v1 <- levelplot(BG3_hmap_comp,
  at = custom_at,
  main = list(label = "Predicted Region Annotation Comparisons BG3", cex = 2.4),
  xlab = list(label = "log2 Enrichment Difference", cex = 1.8),
  ylab = list(label = "Predicted By", cex = 1.8),
  scales=list(x=list(rot=90, cex=1.2), y=list(cex=1.2)),
  col.regions = cols_contrast)

  v2 <- levelplot(S2_hmap_comp,
    at = custom_at,
    main = list(label = "Predicted Region Annotation Comparisons S2", cex = 2.4),
    xlab = list(label = "log2 Enrichment Difference", cex = 1.8),
    ylab = list(label = "Predicted By", cex = 1.8),
    scales=list(x=list(rot=90, cex=1.2), y=list(cex=1.2)),
    col.regions = cols_contrast)

pdf(paste0("/home/jw18713/Paper_replots/Plots/Figure3A_HM.pdf"),
width = 22, height = 8, pointsize = 14)
par(mar=c(4,4,4,10.5), cex = 1.2)
grid.arrange(v1, v2, nrow = 1)
dev.off()
