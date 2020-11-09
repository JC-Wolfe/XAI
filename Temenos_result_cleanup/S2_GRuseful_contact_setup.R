library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)
library(VennDiagram)

load("grown_predicted_enhancers_S2.Rda") # pred_enh
load("/home/jw18713/Project1/data/S2_enriched_contacts_GR.Rda") # gr_contacts
load(file = "/home/jw18713/annotated_GRange_FI_done.Rda") # ann_build

reformatting <- GRanges(seqnames = seqnames(gr_contacts),
                        ranges = IRanges(start(gr_contacts),
                        end(gr_contacts)))
reformatting$contact_chr <- gr_contacts$chr
reformatting$contact_start <- gr_contacts$start
reformatting$contact_end <- gr_contacts$end
reformatting$enrichment <- gr_contacts$enrichment



promoters <- reduce(ann_build[ann_build$annotations == "promoter"])
prm <- rep(F, length(contact_check))

contact_check <- GRanges(seqnames = reformatting$contact_chr,
                          ranges = IRanges(reformatting$contact_start,
                          reformatting$contact_end))
seqlevelsStyle(contact_check) <- "UCSC"

overlaps <- findOverlaps(promoters, contact_check)
prm[subjectHits(overlaps)] <- T

reformatting$promoter_contact <- prm

gr_contacts <- reformatting

save(gr_contacts, file = "useful_HiC_dm6_contactMap_S2.Rda")

load("/home/jw18713/Project1/data/model_predictions/S2_FuzzyGR.Rda") #grown_predicted_enhancers_S2
starr <- reduce(gro[gro$STARR_seq_binary == 1])
save(starr, file = "starr_S2.Rda")
