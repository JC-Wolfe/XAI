
library(GenomicRanges)
library(rtracklayer)

.expression_extraction <- function(gr, expression_data){
  ovr <- findOverlaps(gr, expression_data)
  gr <- gr[queryHits(ovr)]
  gr$contact_chr <- expression_data[subjectHits(ovr)]$contact_chr
  gr$contact_start <- expression_data[subjectHits(ovr)]$contact_start
  gr$contact_end <- expression_data[subjectHits(ovr)]$contact_end
  gr$Gene_short_name <- expression_data[subjectHits(ovr)]$Gene_short_name
  gr$FPKM <- expression_data[subjectHits(ovr)]$FPKM
  return(gr)
}

distal_expression <- function(x, g_exp, distal_contacts, line){

  expression <- GRanges(seqnames = g_exp$Chromosome.arm, ranges=IRanges(g_exp$Start, g_exp$End))
  expression$Gene_short_name <- g_exp$Gene.short.name

  if (line == "BG3"){
    expression$FPKM <- g_exp$BG3_FPKM
  } else if (line == "S2"){
    expression$FPKM <- g_exp$S2DRSC_FPKM
  }

  # So lets see where we have promoter contacts then
  exp_test <- GRanges(seqnames = distal_contacts$contact_chr, ranges = IRanges(distal_contacts$contact_start, distal_contacts$contact_end))
  exp_over <- findOverlaps(exp_test, expression)
  contacted_expressed <- distal_contacts[queryHits(exp_over)]
  contacted_expressed$Gene_short_name <- expression$Gene_short_name[subjectHits(exp_over)]
  contacted_expressed$FPKM <- expression$FPKM[subjectHits(exp_over)]

  # This is a GRange list that will contain all of the target enhancers position by position
  list_distal_only <- split(x, as.factor(x))

  old <- Sys.time()
  glist <- lapply(list_distal_only, .expression_extraction, expression_data = contacted_expressed)
  new <- Sys.time()
  print(old - new)
  return(glist)
}

# Proximal functions are down here

proximal_expression <- function(x, g_exp, line){

  expression <- GRanges(seqnames = g_exp$Chromosome.arm, ranges=IRanges(g_exp$Start, g_exp$End))
  expression$Gene_short_name <- g_exp$Gene.short.name

  if (line == "BG3"){
    expression$FPKM <- g_exp$BG3_FPKM
  } else if (line == "S2"){
    expression$FPKM <- g_exp$S2DRSC_FPKM
  }

  # So lets see where we have promoter contacts then
  enh_reach <- GRanges(seqnames = seqnames(x), ranges = IRanges(start(x)-5000, end(x)+5000))
  exp_proms <- expression
  start(exp_proms) <- start(exp_proms)-250
  end(exp_proms) <- start(exp_proms)+250
  exp_over <- findOverlaps(enh_reach, exp_proms)
  contacted_expressed <- enh_reach[queryHits(exp_over)]
  contacted_expressed$Gene_short_name <- exp_proms$Gene_short_name[subjectHits(exp_over)]
  contacted_expressed$FPKM <- exp_proms$FPKM[subjectHits(exp_over)]

  # This is a GRange list that will contain all of the target enhancers position by position
  list_proximal <- split(x, as.factor(x))

  old <- Sys.time()
  glist <- lapply(list_proximal, .expression_extraction, expression_data = contacted_expressed)
  new <- Sys.time()
  print(old - new)
  return(glist)
}
