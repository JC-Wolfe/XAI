
# Get the functions
source("~/Project1/Paper_Source_Scripts/glist_functions.R")

# Load the BG3 objects
# Proximal and Distal is 1
# Distal Only is 2
# Proximal Only is 3
# Neither is 4
# So distal_only will be list[[2]]
# proximal_all will be sort(c(list[[1]],list[[3]]))

# STARR-seq
load("BG3_starr_res_list.Rda")
load("S2_starr_res_list.Rda")

# Loading common resources

# loading p5k (BG3 promoter contacts over 5kb away for distal function)
BG3_distal_contacts <- get(load("/home/jw18713/Project1/Paper_Rda_Stuff/BG3/over5k_contacts.Rda"))
# loading p5k for S2
S2_distal_contacts <- get(load("/home/jw18713/Project1/Paper_Rda_Stuff/S2/over5k_contacts.Rda"))

# expression data for both cell lines
expression_data <- read.csv("/home/jw18713/Project1/data/gene_expression_data.csv")

# And finally, running the things themselves
# Putative Enhancers

# For BG3 starr enhancers
# Distal
BG3_starr_dist <- distal_expression(BG3_starr_res[[2]], g_exp = expression_data,
  distal_contacts = BG3_distal_contacts, line = "BG3")
save(BG3_starr_dist, file = "~/rerun_glists/starr/BG3_starr_dist.Rda")
# Proximal
BG3_starr_prox <- proximal_expression(sort(c(BG3_starr_res[[1]],BG3_starr_res[[3]])),
  g_exp = expression_data, line = "BG3")
save(BG3_starr_prox, file = "~/rerun_glists/starr/BG3_starr_prox.Rda")

# For S2 starr enhancers
# Distal
S2_starr_dist <- distal_expression(S2_starr_res[[2]], g_exp = expression_data,
  distal_contacts = S2_distal_contacts, line = "S2")
save(S2_starr_dist, file = "~/rerun_glists/starr/S2_starr_dist.Rda")
# Proximal
S2_starr_prox <- proximal_expression(sort(c(S2_starr_res[[1]],S2_starr_res[[3]])),
  g_exp = expression_data, line = "S2")
save(S2_starr_prox, file = "~/rerun_glists/starr/S2_starr_prox.Rda")
