
# Importing libraries
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)
predicted <- get(load("/home/jw18713/Project1/data/model_predictions/S2_FuzzyGR.Rda"))

.sumReadsConfPerc1 <- function(methylationData){   #Function to compute sums within below function
  return(sum(methylationData$'Conf..Perc.1'))
}
.ConfPercAvg <- function(methylationData, regions){

  overlaps <- findOverlaps(methylationData, regions, ignore.strand = TRUE)
  methylationDataContextList <- IRanges::splitAsList(methylationData[queryHits(overlaps)],  subjectHits(overlaps))
  regionsIndexes <- as.integer(names(methylationDataContextList))

  regions$sumConfPerc1 <- rep(0, times=length(regions)) #Confidence percentage 1

  if(length(regionsIndexes) > 0){
    regions$sumConfPerc1[regionsIndexes] <- sapply(methylationDataContextList,.sumReadsConfPerc1) #Confidence percentage 1
  }

  regions$avgConfPerc1 <- regions$sumConfPerc1 / (width(regions)/10)

  return(regions)
}

over80 <- predicted[predicted$Conf..Perc.1 >= 0.8]
initial_regions <- .ConfPercAvg(predicted, reduce(over80))

combined <- reduce(over80, min.gapwidth = 100)
overlaps <- findOverlaps(initial_regions, combined)

glist <- split(initial_regions, subjectHits(overlaps))

save(glist, file="grangelist_S2.Rda")



loop_data <- rep(0, len = glist)
for (i in seq_along(glist)){
print(i)
loop_data[i] <- length(glist[[i]])
}
save(loop_data, file = "loop_data100_S2.Rda")
