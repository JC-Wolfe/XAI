# Get the functions
source("~/dm6_promoter_correction/Scripts/Functional/glist_functions.R")

# Load the BG3 objects
# Proximal and Distal is 1
# Distal Only is 2
# Proximal Only is 3
# Neither is 4
# So distal_only will be list[[2]]
# proximal_all will be sort(c(list[[1]],list[[3]]))


# Background
load("~/dm6_promoter_correction/Class_split_glists/BG3_bg_list.Rda")
load("~/dm6_promoter_correction/Class_split_glists/S2_bg_list.Rda")

# Loading common resources

# loading p5k (BG3 promoter contacts over 5kb away for distal function)
BG3_distal_contacts <- get(load("~/dm6_promoter_correction/BG3_over5k_contacts.Rda"))
# loading p5k for S2
S2_distal_contacts <- get(load("~/dm6_promoter_correction/S2_over5k_contacts.Rda"))

# expression data for both cell lines
expression_data <- read.csv("~/dm6_promoter_correction/gene_expression_data.csv")

# And finally, running the things themselves
# Putative Enhancers

# For BG3 bg enhancers
# Distal
BG3_bg_dist <- distal_expression(BG3_bg[[2]], g_exp = expression_data,
  distal_contacts = BG3_distal_contacts, line = "BG3")
save(BG3_bg_dist, file = "~/dm6_promoter_correction/Full_glists/bg/BG3_bg_dist.Rda")
# Proximal
BG3_bg_prox <- proximal_expression(sort(c(BG3_bg[[1]],BG3_bg[[3]])),
  g_exp = expression_data, line = "BG3")
save(BG3_bg_prox, file = "~/dm6_promoter_correction/Full_glists/bg/BG3_bg_prox.Rda")

# For S2 bg enhancers
# Distal
S2_bg_dist <- distal_expression(S2_bg[[2]], g_exp = expression_data,
  distal_contacts = S2_distal_contacts, line = "S2")
save(S2_bg_dist, file = "~/dm6_promoter_correction/Full_glists/bg/S2_bg_dist.Rda")
# Proximal
S2_bg_prox <- proximal_expression(sort(c(S2_bg[[1]],S2_bg[[3]])),
  g_exp = expression_data, line = "S2")
save(S2_bg_prox, file = "~/dm6_promoter_correction/Full_glists/bg/S2_bg_prox.Rda")
