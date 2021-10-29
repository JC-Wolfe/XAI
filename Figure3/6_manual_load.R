library(GenomicRanges)
library(rtracklayer)

source("~/dm6_promoter_correction/Scripts/Functional/expression_analysis_functon.R")

#-----------------------/Loading Putative Files/-------------------------------#
#------------------------------------------------------------------------------#
# BG3
# Distal
load("~/dm6_promoter_correction/Full_glists/putative/BG3_put_dist.Rda")
BG3_putative_distal_only <- expression_analysis(BG3_put_dist)
rm(BG3_put_dist)
# Proximal
load("~/dm6_promoter_correction/Full_glists/putative/BG3_put_prox.Rda")
BG3_putative_proximal <- expression_analysis(BG3_put_prox)
rm(BG3_put_prox)

# S2
# Distal
load("~/dm6_promoter_correction/Full_glists/putative/S2_put_dist.Rda")
S2_putative_distal_only <- expression_analysis(S2_put_dist)
rm(S2_put_dist)
# Proximal
load("~/dm6_promoter_correction/Full_glists/putative/S2_put_prox.Rda")
S2_putative_proximal <- expression_analysis(S2_put_prox)
rm(S2_put_prox)
#------------------------/End of Putative Files/-------------------------------#
#------------------------------------------------------------------------------#
save(BG3_putative_distal_only,
  file = "~/dm6_promoter_correction/Expression_lists/BG3_putative_distal_only.Rda")

save(BG3_putative_proximal,
  file = "~/dm6_promoter_correction/Expression_lists/BG3_putative_proximal.Rda")

save(S2_putative_distal_only,
  file = "~/dm6_promoter_correction/Expression_lists/S2_putative_distal_only.Rda")

save(S2_putative_proximal,
  file = "~/dm6_promoter_correction/Expression_lists/S2_putative_proximal.Rda")

#------------------------/Loading Background Files/----------------------------#
#------------------------------------------------------------------------------#
# BG3 background
# Distal
load("~/dm6_promoter_correction/Full_glists/bg/BG3_bg_dist.Rda")
BG3_background_distal_only <- expression_analysis(BG3_bg_dist)
rm(BG3_bg_dist)
# Proximal
load("~/dm6_promoter_correction/Full_glists/bg/BG3_bg_prox.Rda")
BG3_background_proximal <- expression_analysis(BG3_bg_prox)
rm(BG3_bg_prox)

# S2 Background
# Distal
load("~/dm6_promoter_correction/Full_glists/bg/S2_bg_dist.Rda")
S2_background_distal_only <- expression_analysis(S2_bg_dist)
rm(S2_bg_dist)
# Proximal
load("~/dm6_promoter_correction/Full_glists/bg/S2_bg_prox.Rda")
S2_background_proximal <- expression_analysis(S2_bg_prox)
rm(S2_bg_prox)
#------------------------/End of Background Files/-----------------------------#
#------------------------------------------------------------------------------#
save(BG3_background_distal_only,
  file = "~/dm6_promoter_correction/Expression_lists/BG3_background_distal_only.Rda")

save(BG3_background_proximal,
  file = "~/dm6_promoter_correction/Expression_lists/BG3_background_proximal.Rda")

save(S2_background_distal_only,
  file = "~/dm6_promoter_correction/Expression_lists/S2_background_distal_only.Rda")

save(S2_background_proximal,
  file = "~/dm6_promoter_correction/Expression_lists/S2_background_proximal.Rda")

#--------------------------/Loading Common Files/------------------------------#
#------------------------------------------------------------------------------#
# BG3
# Distal
load("~/dm6_promoter_correction/Full_glists/common/BG3_common_dist.Rda")
BG3_common_distal_only <- expression_analysis(BG3_common_dist)
rm(BG3_common_dist)
# Proximal
load("~/dm6_promoter_correction/Full_glists/common/BG3_common_prox.Rda")
BG3_common_proximal <- expression_analysis(BG3_common_prox)
rm(BG3_common_prox)

# S2
# Distal
load("~/dm6_promoter_correction/Full_glists/common/S2_common_dist.Rda")
S2_common_distal_only <- expression_analysis(S2_common_dist)
rm(S2_common_dist)
# Proximal
load("~/dm6_promoter_correction/Full_glists/common/S2_common_prox.Rda")
S2_common_proximal <- expression_analysis(S2_common_prox)
rm(S2_common_prox)
#--------------------------/End of Common Files/-------------------------------#
#------------------------------------------------------------------------------#
save(BG3_common_distal_only,
  file = "~/dm6_promoter_correction/Expression_lists/BG3_common_distal_only.Rda")

save(BG3_common_proximal,
  file = "~/dm6_promoter_correction/Expression_lists/BG3_common_proximal.Rda")

save(S2_common_distal_only,
  file = "~/dm6_promoter_correction/Expression_lists/S2_common_distal_only.Rda")

save(S2_common_proximal,
  file = "~/dm6_promoter_correction/Expression_lists/S2_common_proximal.Rda")

#----------------------/Loading STARR-seq Only Files/--------------------------#
#------------------------------------------------------------------------------#
# BG3
# Distal
load("~/dm6_promoter_correction/Full_glists/starr/BG3_starr_dist.Rda")
BG3_starr_distal_only <- expression_analysis(BG3_starr_dist)
rm(BG3_starr_dist)
# Proximal
load("~/dm6_promoter_correction/Full_glists/starr/BG3_starr_prox.Rda")
BG3_starr_proximal <- expression_analysis(BG3_starr_prox)
rm(BG3_starr_prox)

# S2
# Distal
load("~/dm6_promoter_correction/Full_glists/starr/S2_starr_dist.Rda")
S2_starr_distal_only <- expression_analysis(S2_starr_dist)
rm(S2_starr_dist)
# Proximal
load("~/dm6_promoter_correction/Full_glists/starr/S2_starr_prox.Rda")
S2_starr_proximal <- expression_analysis(S2_starr_prox)
rm(S2_starr_prox)
#----------------------/End of STARR-seq Only Files/---------------------------#
#------------------------------------------------------------------------------#
save(BG3_starr_distal_only,
  file = "~/dm6_promoter_correction/Expression_lists/BG3_starr_distal_only.Rda")

save(BG3_starr_proximal,
  file = "~/dm6_promoter_correction/Expression_lists/BG3_starr_proximal.Rda")

save(S2_starr_distal_only,
  file = "~/dm6_promoter_correction/Expression_lists/S2_starr_distal_only.Rda")

save(S2_starr_proximal,
  file = "~/dm6_promoter_correction/Expression_lists/S2_starr_proximal.Rda")
