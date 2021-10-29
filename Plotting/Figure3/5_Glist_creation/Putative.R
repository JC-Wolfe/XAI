# Get the functions
source("~/dm6_promoter_correction/Scripts/Functional/glist_functions.R")

# Load the BG3 objects
# Proximal and Distal is 1
# Distal Only is 2
# Proximal Only is 3
# Neither is 4
# So distal_only will be list[[2]]
# proximal_all will be sort(c(list[[1]],list[[3]]))

# Putative
load("~/dm6_promoter_correction/Class_split_glists/BG3_putative_list.Rda")
load("~/dm6_promoter_correction/Class_split_glists/S2_putative_list.Rda")

# Loading common resources

# loading p5k (BG3 promoter contacts over 5kb away for distal function)
BG3_distal_contacts <- get(load("~/dm6_promoter_correction/BG3_over5k_contacts.Rda"))
# loading p5k for S2
S2_distal_contacts <- get(load("~/dm6_promoter_correction/S2_over5k_contacts.Rda"))

# expression data for both cell lines
expression_data <- read.csv("~/dm6_promoter_correction/gene_expression_data.csv")

# And finally, running the things themselves
# Putative Enhancers

# For BG3 putative enhancers
# Distal
BG3_put_dist <- distal_expression(BG3_putative[[2]], g_exp = expression_data,
  distal_contacts = BG3_distal_contacts, line = "BG3")
save(BG3_put_dist, file = "~/dm6_promoter_correction/Full_glists/putative/BG3_put_dist.Rda")
# Proximal
BG3_put_prox <- proximal_expression(sort(c(BG3_putative[[1]],BG3_putative[[3]])),
  g_exp = expression_data, line = "BG3")
save(BG3_put_prox, file = "~/dm6_promoter_correction/Full_glists/putative/BG3_put_prox.Rda")

# For S2 putative enhancers
# Distal
S2_put_dist <- distal_expression(S2_putative[[2]], g_exp = expression_data,
  distal_contacts = S2_distal_contacts, line = "S2")
save(S2_put_dist, file = "~/dm6_promoter_correction/Full_glists/putative/S2_put_dist.Rda")
# Proximal
S2_put_prox <- proximal_expression(sort(c(S2_putative[[1]],S2_putative[[3]])),
  g_exp = expression_data, line = "S2")
save(S2_put_prox, file = "~/dm6_promoter_correction/Full_glists/putative/S2_put_prox.Rda")
