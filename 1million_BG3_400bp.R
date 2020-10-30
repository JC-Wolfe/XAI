# Importing libraries
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)

all_data <- get(load("/home/jw18713/Project1/R_objects/BG3_tiles_400bp.Rda")) #Loading in full S2 dataset as all_data


# I need to figure out the ratios of chromosomes to be included in the training and testing data
chr_ratios <- width(seqnames(all_data))/sum(width(seqnames(all_data)))

# Now that I have the ratios I need to figure out how many samples I need from each chromosome to
# proportionally represent them in the 1 million data points I will be using. I'm rounding to be
# multiples of 100
nsamples <- round(chr_ratios * 1000000)
pos_ratio <- sum(all_data$STARR_seq_binary)/length(all_data)
npos <- round(nsamples * pos_ratio)
nneg <- nsamples - npos

# Now I need to build a set of samples to work with

# Splitting the chromosomes for ease of sampling
chr_split <- split(all_data, seqnames(all_data))
chr_split_resampled <- GRanges()


# Now I need to randomly sample the correct ratio of positive and negative samples and put them into a GRange object
for (i in seq_along(chr_split)){
  pos_neg_split <- split(chr_split[[i]], chr_split[[i]]$STARR_seq_binary)
  neg_STARR <- pos_neg_split[[1]]
  pos_STARR <- pos_neg_split[[2]]
  p1 <- sample(pos_STARR, npos[i])
  n1 <- sample(neg_STARR, nneg[i])
  chr_split_resampled <- c(chr_split_resampled, p1, n1)
}

# Resorting the new GRange object
chr_split_resampled <- sortSeqlevels(chr_split_resampled) #Seqlevels must be correctly sorted for sort to work!
chr_split_resampled <- sort(chr_split_resampled)

save(chr_split_resampled, file = "/home/jw18713/Project1/R_objects/1million_norm_BG3_400bp.Rda")
write.csv(as.data.frame(chr_split_resampled), file = "/home/jw18713/Project1/data/1million_norm_BG3_400bp.csv", row.names = F)
