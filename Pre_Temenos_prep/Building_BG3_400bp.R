# Importing libraries
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)

#Setting the right working directory (S2_histone)
setwd("/home/jw18713/machine_learning_project/datasets/BG3_histone")

# Creating a list of all of the folders that I need (excluding non .wig files)
folders <- dir()

# Using paste0 to create a list of file paths to folders containing each .wig file
paths <- paste0(getwd(),"/",folders,"/signal_data_files/")

# Creating a list of "0"s to be filled with the relevant files and file names
files <- rep("0", length(folders))
name <- rep("0", length(folders))

# A for loop to loop over every directory contained in paths that:
# Sets the working directory to the current element in paths
# Locates the smoothed file using grep and saves it to a temporary variable
# Changes the "i"th element of names to the correct dataset name
# Changes the "i"th element of files to the path to the relevant smoothed .wig file
for(i in seq_along(paths)){
setwd(paths[i])
smooth <- dir()[grep("smoothedM",dir())]
name[i] <- strsplit(smooth,":")[[1]][1]
files[i] <- paste0(paths[i],smooth)
}

# Getting the dm6 genome
genome <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6)

# Choosing a bin size
binsize <- 10

# Creation of the GRange object to contain ChIP metadata
chrs <- seqnames(Dmelanogaster)[1:5]
dm6_tiles <- GRanges()
seqlevels(dm6_tiles) <- chrs
for(i in 1:length(chrs)){
  print(chrs[i])
  start <- seq(1, length(Dmelanogaster[[i]]), by=binsize)
  buffer <- GRanges(seqnames=chrs[i], ranges=IRanges(start = start, end = start + 9), strand="*")
  dm6_tiles <- c(dm6_tiles, buffer)
}

# Chain for lift over
chain <- import.chain("/home/jw18713/machine_learning_project/dm3ToDm6.over.chain")

for(i in seq_along(files)){
  gro <- readGeneric(files[i], skip=1, sep=" ", meta.cols=list(score=4))
  seqlevelsStyle(gro) <- "UCSC"
  gro <- unlist(liftOver(gro, chain))
  # WE HAVE LIFT OVER!
  overlaps <- findOverlaps(dm6_tiles,gro)
  scores <- rep(0, length(dm6_tiles))
  scores[as.vector(queryHits(overlaps))] <- gro$score[as.vector(subjectHits(overlaps))]
  print(name[i])
  mcols(dm6_tiles)[i] <- scores
  names(mcols(dm6_tiles))[i] <- name[i]
}

STARRrep1 <- read.delim("/home/jw18713/machine_learning_project/datasets/BG3_starr_seq/BG3_peakSummits.txt", header = T, sep = "\t")
STARRgro <- GRanges(seqnames = STARRrep1$chr, ranges = IRanges(STARRrep1$summit-200, STARRrep1$summit+200)) #Creating a GRange object for STARR seq data (summit +- 100 bases)
STARRgro <- unlist(liftOver(STARRgro, chain))
overlaps <- findOverlaps(dm6_tiles, STARRgro)
scores <- rep(0, length(dm6_tiles))
scores[as.vector(queryHits(overlaps))] <- 1
dm6_tiles$STARR_seq_binary <- scores

sortdf <- mcols(dm6_tiles)[,sort(names(mcols(dm6_tiles)))]
mcols(dm6_tiles) <- sortdf

# Getting min and max values to de-normalise outputs later
output_max <- max(mcols(dm6_tiles)[length(mcols(dm6_tiles))][,1])
output_min <- min(mcols(dm6_tiles)[length(mcols(dm6_tiles))][,1])
save(output_max, file="/home/jw18713/Project1/R_objects/BG3_denorm_max_400bp.Rda")
save(output_min, file="/home/jw18713/Project1/R_objects/BG3_denorm_max_400bp.Rda")

# Scaling the dataset
for (c in seq_along(mcols(dm6_tiles))){
  max_value = max(mcols(dm6_tiles)[c][,1])
  min_value = min(mcols(dm6_tiles)[c][,1])
  mcols(dm6_tiles)[c][,1] <- (mcols(dm6_tiles)[c][,1] - min_value) / (max_value - min_value)
}


save(dm6_tiles, file="/home/jw18713/Project1/R_objects/BG3_tiles_400bp.Rda")

dfwrite <- as.data.frame(dm6_tiles)
write.csv(dfwrite, file="/home/jw18713/Project1/data/norm_BG3_400bp.csv", row.names=F)
