library(GenomicRanges)
library(rtracklayer)

load("/home/jw18713/Project1/data/BG3_enriched_contacts_GR.Rda") # gr_contacts
load("/home/jw18713/dm6_promoter_correction/dm6_annotations/dm6_annotation.Rda")

reformatting <- GRanges(seqnames = seqnames(gr_contacts),
                        ranges = IRanges(start(gr_contacts),
                        end(gr_contacts)))
reformatting$contact_chr <- gr_contacts$chr
reformatting$contact_start <- gr_contacts$start
reformatting$contact_end <- gr_contacts$end
reformatting$enrichment <- gr_contacts$enrichment

contact_check <- GRanges(seqnames = reformatting$contact_chr,
                          ranges = IRanges(reformatting$contact_start,
                          reformatting$contact_end))
seqlevelsStyle(contact_check) <- "UCSC"

promoters <- reduce(dm6_annotation[dm6_annotation$annotations == "promoter"])
prm <- rep(F, length(contact_check))

overlaps <- findOverlaps(promoters, contact_check)
prm[subjectHits(overlaps)] <- T

reformatting$promoter_contact <- prm

gr_contacts <- reformatting

save(gr_contacts, file = "~/dm6_promoter_correction/Contact_Maps/useful_HiC_dm6_contactMap_BG3.Rda")
