library(GenomicRanges)
library(rtracklayer)

gr_contacts <- get(load("~/dm6_promoter_correction/Contact_Maps/useful_HiC_dm6_contactMap_BG3.Rda"))

gr_contacts$contact_chr <- paste0("chr", as.character(gr_contacts$contact_chr))

contact_over5k <- gr_contacts[end(gr_contacts) + 5000 < gr_contacts$contact_start
                              |!as.character(seqnames(gr_contacts)) == gr_contacts$contact_chr]
p5k <- contact_over5k[contact_over5k$promoter_contact == TRUE] # Contact promoter over 5kb from enhancer

contact_under5k <- gr_contacts[end(gr_contacts) + 5000 > gr_contacts$contact_start
                              & as.character(seqnames(gr_contacts)) == gr_contacts$contact_chr]
p_under5k <- contact_under5k[contact_under5k$promoter_contact == TRUE] # Contact promoter under 5kb from enhancer

save(p5k, file = "~/dm6_promoter_correction/Objects/BG3_over5k_contacts.Rda")
save(p_under5k, file = "~/dm6_promoter_correction/Objects/BG3_under5k_contacts.Rda")
