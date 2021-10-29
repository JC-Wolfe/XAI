library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)

# Getting the dm6 genome
genome <- getSeq(BSgenome.Dmelanogaster.UCSC.dm6)

# Choosing a bin size
binsize <- 10

# Creation of the GRange object to contain ChIP metadata
chrs <- seqnames(Dmelanogaster)[1:5]
ann_build <- GRanges()
seqlevels(ann_build) <- chrs
for(i in 1:length(chrs)){
  print(chrs[i])
  start <- seq(1, length(Dmelanogaster[[i]]), by=binsize)
  buffer <- GRanges(seqnames=chrs[i], ranges=IRanges(start = start, end = start + 9), strand="*")
  ann_build <- c(ann_build, buffer)
}

# Importing the annotated GRange object
gff3 <- import("~/dm6_promoter_correction/dm6_annotations/Drosophila_melanogaster.BDGP6.22.96.chr.gff3.gz")
seqlevelsStyle(gff3) = "UCSC" #Setting sequence levels to UCSC
gff3 <- gff3[seqnames(gff3) %in% seqnames(Dmelanogaster)[1:5]]

# Creating and filling the annotation meta column
ann_build$annotations <- "intergenic"

genes <- gff3[gff3$type == "gene"]
exons <- gff3[gff3$type == "exon"]
fiveprime <- gff3[gff3$type == "five_prime_UTR"]
threeprime <- gff3[gff3$type == "three_prime_UTR"]
sense_promoter <- gff3[gff3$type == "gene" & strand(gff3) == "+"]
start(sense_promoter) <- start(sense_promoter) - 250
end(sense_promoter) <- start(sense_promoter) + 249
antisense_promoter <- gff3[gff3$type == "gene" & strand(gff3) == "-"]
end(antisense_promoter) <- end(antisense_promoter) + 250
start(antisense_promoter) <- end(antisense_promoter) - 249
promoters <- reduce(c(sense_promoter, antisense_promoter))

overlaps1 <- findOverlaps(ann_build, genes) #Finding where the mRNA annotations are
ann_build$annotations[queryHits(overlaps1)] <- "mRNA"
overlaps2 <- findOverlaps(ann_build, exons)
ann_build$annotations[queryHits(overlaps2)] <- "exon"
ann_build$annotations[ann_build$annotations == "mRNA"] <- "intron"
overlaps3 <- findOverlaps(ann_build, fiveprime)
ann_build$annotations[queryHits(overlaps3)] <- "five_prime_UTR"
overlaps4 <- findOverlaps(ann_build, threeprime)
ann_build$annotations[queryHits(overlaps4)] <- "three_prime_UTR"
overlaps5 <- findOverlaps(ann_build, promoters)
ann_build$annotations[queryHits(overlaps5)] <- "promoter"

# Attempting First Intron Stuff
sense_transcripts <- gff3[gff3$type == "mRNA" & strand(gff3) == "+"]
antisense_transcripts <- gff3[gff3$type == "mRNA" & strand(gff3) == "-"]

# All introns in the annoataion
all_introns <- reduce(ann_build[ann_build$annotations == "intron"])

# Overlaps of sense and antisense transcripts and introns
sense_over <- findOverlaps(sense_transcripts, all_introns)
antisense_over <- findOverlaps(antisense_transcripts, all_introns)

# Individual transcipts for finding introns
sense_targets <- unique(queryHits(sense_over))
antisense_targets <- unique(queryHits(antisense_over))

# Object to hold sense first introns
sense_fi <- GRanges()
antisense_fi <- GRanges()

# Finding first introns by unique subject hits
for (i in seq(1, length(sense_targets))){
  p <- all_introns[subjectHits(sense_over[queryHits(sense_over)==sense_targets[i]])]
  sense_fi <- c(sense_fi, p[1])
}

# And first introns for antisense genes (look like last introns)
for (i in seq(1, length(antisense_targets))){
  p <- all_introns[subjectHits(antisense_over[queryHits(antisense_over)==antisense_targets[i]])]
  antisense_fi <- c(antisense_fi, p[length(p)])
}

first_introns <- reduce(c(sense_fi, antisense_fi))

overlaps6 <- findOverlaps(ann_build, first_introns)
ann_build$annotations[queryHits(overlaps6)] <- "first_intron"

dm6_annotation <- ann_build

save(dm6_annotation, file = "~/dm6_promoter_correction/dm6_annotations/dm6_annotation.Rda")
