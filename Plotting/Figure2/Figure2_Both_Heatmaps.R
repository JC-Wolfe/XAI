#	Heat map for marks BG3
#	Heat map for marks S2


library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)
library(gplots)
library(lattice)


BG3_results <- get(load("/home/jw18713/Project1/data/model_predictions/BG3_FuzzyGR.Rda"))
S2_results <- get(load("/home/jw18713/Project1/data/model_predictions/S2_FuzzyGR.Rda"))
marks <- mcols(BG3_results)[1:26]
S2_enhancers <- get(load("grown_predicted_enhancers_S2.Rda"))
S2_starr <- get(load("starr_S2.Rda"))
BG3_enhancers <- get(load("grown_predicted_enhancers_BG3.Rda"))
BG3_starr <- get(load("starr.Rda"))

BG3_starr_and_predicted <- reduce(c(BG3_enhancers[BG3_enhancers %over% BG3_starr],
                                BG3_starr[BG3_starr %over% BG3_enhancers]))
BG3_fuzzy <- BG3_enhancers[!BG3_enhancers %over% BG3_starr]
BG3_starr_specific <- BG3_starr[!BG3_starr %over% BG3_enhancers]



BG3_both <- BG3_results[BG3_results %over% starr_and_predicted]
BG3_fuzzy_only <- BG3_results[BG3_results %over% BG3_fuzzy]
BG3_STARR_only <- BG3_results[BG3_results %over% BG3_starr_specific]
BG3_neither <- BG3_results[!BG3_results %over% BG3_both & !BG3_results %over% BG3_fuzzy_only & !BG3_results %over% BG3_STARR_only]


BG3_map_matrix <- matrix(0, nrow = 4, ncol=length(marks))
colnames(BG3_map_matrix) <- names(marks)
rownames(BG3_map_matrix) <- c("Both", "XAI", "STARR-seq", "Neither")

BG3_map_matrix[1,] <- apply(mcols(BG3_both)[1:26], 2, mean)
BG3_map_matrix[2,] <- apply(mcols(BG3_fuzzy_only)[1:26], 2, mean)
BG3_map_matrix[3,] <- apply(mcols(BG3_STARR_only)[1:26], 2, mean)
BG3_map_matrix[4,] <- apply(mcols(BG3_neither)[1:26], 2, mean)

BG3_expected <- apply(mcols(BG3_results)[1:26], 2, mean)

BG3_map_matrix[1,] <- log2(BG3_map_matrix[1,]/BG3_expected)
BG3_map_matrix[2,] <- log2(BG3_map_matrix[2,]/BG3_expected)
BG3_map_matrix[3,] <- log2(BG3_map_matrix[3,]/BG3_expected)
BG3_map_matrix[4,] <- log2(BG3_map_matrix[4,]/BG3_expected)

# Everything below here is S2 stuff

S2_starr_and_predicted <- reduce(c(S2_enhancers[S2_enhancers %over% S2_starr],
                                S2_starr[S2_starr %over% S2_enhancers]))
S2_fuzzy <- S2_enhancers[!S2_enhancers %over% S2_starr]
S2_starr_specific <- S2_starr[!S2_starr %over% S2_enhancers]



S2_both <- S2_results[S2_results %over% starr_and_predicted]
S2_fuzzy_only <- S2_results[S2_results %over% S2_fuzzy]
S2_STARR_only <- S2_results[S2_results %over% S2_starr_specific]
S2_neither <- S2_results[!S2_results %over% S2_both & !S2_results %over% S2_fuzzy_only & !S2_results %over% S2_STARR_only]


S2_map_matrix <- matrix(0, nrow = 4, ncol=length(marks))
colnames(S2_map_matrix) <- names(marks)
rownames(S2_map_matrix) <- c("Both", "XAI", "STARR-seq", "Neither")

S2_map_matrix[1,] <- apply(mcols(S2_both)[1:26], 2, mean)
S2_map_matrix[2,] <- apply(mcols(S2_fuzzy_only)[1:26], 2, mean)
S2_map_matrix[3,] <- apply(mcols(S2_STARR_only)[1:26], 2, mean)
S2_map_matrix[4,] <- apply(mcols(S2_neither)[1:26], 2, mean)

S2_expected <- apply(mcols(S2_results)[1:26], 2, mean)

S2_map_matrix[1,] <- log2(S2_map_matrix[1,]/S2_expected)
S2_map_matrix[2,] <- log2(S2_map_matrix[2,]/S2_expected)
S2_map_matrix[3,] <- log2(S2_map_matrix[3,]/S2_expected)
S2_map_matrix[4,] <- log2(S2_map_matrix[4,]/S2_expected)




custom_at <- seq(-0.5, 0.5, by = (1/15))


cols_contrast <- rev(c(colorRampPalette(c(cbbPalette[7], "white", cbbPalette[6]))(15)))

pdf(file = "/home/jw18713/Project1/Paper_Plots/Figure2/Figure2A.pdf")
levelplot(t(BG3_map_matrix),
  at = custom_at,
  main = "BG3 Histone Modification Enrichment Change",
  xlab = "log2 Enrichment Difference",
  ylab = "Predicted By",
  scales=list(x=list(rot=90)),
  col.regions = cols_contrast)
dev.off()


pdf(file = "/home/jw18713/Project1/Paper_Plots/Figure2/Figure2B.pdf")
  levelplot(t(S2_map_matrix),
    at = custom_at,
    main = "S2 Histone Modification Enrichment Change",
    xlab = "log2 Enrichment Difference",
    ylab = "Predicted By",
    scales=list(x=list(rot=90)),
    col.regions = cols_contrast)
dev.off()

























pdf("heattest.pdf")
heatmap.2(map_matrix,
        Colv=FALSE, # We want to maintain column order
        Rowv=FALSE,
        key=TRUE,
        breaks = seq(0,1,0.02),
        symbreaks=FALSE,
        symkey = FALSE,
        col=rev(Colors),
        main = "Grouped Heatmap",
        dendrogram="none",
        scale="none", # because we already normalize
        trace="none",
        density.info="none",
        margins = c(6,8),
        #lmat = lmatm,
        #lwid = c(1.5, 4, 14),
        keysize=2)

rasterImage(image,
            xleft, ybottom, xright, ytop,
            angle = 0, interpolate = TRUE, …)
dev.off()
