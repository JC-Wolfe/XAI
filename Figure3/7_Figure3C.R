load("~/dm6_promoter_correction/Expression_lists/S2_starr_proximal.Rda")
load("~/dm6_promoter_correction/Expression_lists/S2_starr_distal_only.Rda")
load("~/dm6_promoter_correction/Expression_lists/BG3_starr_proximal.Rda")
load("~/dm6_promoter_correction/Expression_lists/BG3_starr_distal_only.Rda")
load("~/dm6_promoter_correction/Expression_lists/S2_common_proximal.Rda")
load("~/dm6_promoter_correction/Expression_lists/S2_common_distal_only.Rda")
load("~/dm6_promoter_correction/Expression_lists/BG3_common_proximal.Rda")
load("~/dm6_promoter_correction/Expression_lists/BG3_common_distal_only.Rda")
load("~/dm6_promoter_correction/Expression_lists/S2_background_proximal.Rda")
load("~/dm6_promoter_correction/Expression_lists/S2_background_distal_only.Rda")
load("~/dm6_promoter_correction/Expression_lists/BG3_background_proximal.Rda")
load("~/dm6_promoter_correction/Expression_lists/BG3_background_distal_only.Rda")
load("~/dm6_promoter_correction/Expression_lists/S2_putative_proximal.Rda")
load("~/dm6_promoter_correction/Expression_lists/S2_putative_distal_only.Rda")
load("~/dm6_promoter_correction/Expression_lists/BG3_putative_proximal.Rda")
load("~/dm6_promoter_correction/Expression_lists/BG3_putative_distal_only.Rda")
load("~/dm6_promoter_correction/Neither_matrix.Rda")

all_calc <- function(x, y, neither){
  # Calculating all with contacts over all total (with and without)
  all <- sum(x$count_with_contacts, y$count_with_contacts)/
    sum(x$count_without_contacts, y$count_without_contacts,
    x$count_with_contacts,y$count_with_contacts, neither)
  return(all)
}

BG3_matrix <- matrix(0, 4, 3, dimnames = list(c("Putative", "Common",
  "STARR-seq only", "Whole Genome"), c("Proximal", "Distal Only", "All")))

BG3_matrix[1,1] <- BG3_putative_proximal$contact_proportion
BG3_matrix[1,2] <- BG3_putative_distal_only$contact_proportion
BG3_matrix[1,3] <- all_calc(BG3_putative_proximal, BG3_putative_distal_only,
  neither = Neither_matrix[1,1])

BG3_matrix[2,1] <- BG3_common_proximal$contact_proportion
BG3_matrix[2,2] <- BG3_common_distal_only$contact_proportion
BG3_matrix[2,3] <- all_calc(BG3_common_proximal, BG3_common_distal_only,
  neither = Neither_matrix[1,2])

BG3_matrix[3,1] <- BG3_starr_proximal$contact_proportion
BG3_matrix[3,2] <- BG3_starr_distal_only$contact_proportion
BG3_matrix[3,3] <- all_calc(BG3_starr_proximal, BG3_starr_distal_only,
  neither = Neither_matrix[1,3])

BG3_matrix[4,1] <- BG3_background_proximal$contact_proportion
BG3_matrix[4,2] <- BG3_background_distal_only$contact_proportion
BG3_matrix[4,3] <- all_calc(BG3_background_proximal, BG3_background_distal_only,
  neither = Neither_matrix[1,4])


#---------------------------S2 contact counts----------------------------------#
#------------------------------------------------------------------------------#

S2_matrix <- matrix(0, 4, 3, dimnames = list(c("Putative", "Common",
  "STARR-seq only", "Whole Genome"), c("Proximal", "Distal Only", "All")))
S2_matrix[1,1] <- S2_putative_proximal$contact_proportion
S2_matrix[1,2] <- S2_putative_distal_only$contact_proportion
S2_matrix[1,3] <- all_calc(S2_putative_proximal, S2_putative_distal_only,
  neither = Neither_matrix[2,1])

S2_matrix[2,1] <- S2_common_proximal$contact_proportion
S2_matrix[2,2] <- S2_common_distal_only$contact_proportion
S2_matrix[2,3] <- all_calc(S2_common_proximal, S2_common_distal_only,
  neither = Neither_matrix[2,2])

S2_matrix[3,1] <- S2_starr_proximal$contact_proportion
S2_matrix[3,2] <- S2_starr_distal_only$contact_proportion
S2_matrix[3,3] <- all_calc(S2_starr_proximal, S2_starr_distal_only,
  neither = Neither_matrix[2,3])

S2_matrix[4,1] <- S2_background_proximal$contact_proportion
S2_matrix[4,2] <- S2_background_distal_only$contact_proportion
S2_matrix[4,3] <- all_calc(S2_background_proximal, S2_background_distal_only,
  neither = Neither_matrix[2,4])


#---------------------------------BG3 Stats (All)------------------------------#
#------------------------------------------------------------------------------#
all_counts <- function(x, y, neither){
  without <- sum(x$count_without_contacts, y$count_without_contacts, neither)
  with <- sum(x$count_with_contacts, y$count_with_contacts)
  return(c(without = without, with = with))
}

# Putative
BG3_putative_all <- all_counts(BG3_putative_proximal, BG3_putative_distal_only,
  neither = Neither_matrix[1,1])
# Common
BG3_common_all <- all_counts(BG3_common_proximal, BG3_common_distal_only,
  neither = Neither_matrix[1,2])
# STARR-seq Only
BG3_starr_all <- all_counts(BG3_starr_proximal, BG3_starr_distal_only,
  neither = Neither_matrix[1,3])
# Background
BG3_background_all <- all_counts(BG3_background_proximal, BG3_background_distal_only,
  neither = Neither_matrix[1,4])

BG3_all_stats <- matrix("", 4, 4, dimnames = list(c(rownames(BG3_matrix)),
rownames(BG3_matrix)))

same <- matrix(c(1,1,2,2,3,3,4,4), 4, 2 , byrow = T)

# Filling putative/putative etc with -
BG3_all_stats[same] <- "-"

# Putative/Common
BG3_all_stats[2,1] <- fisher.test(matrix(c(BG3_putative_all, BG3_common_all),
  ncol = 2, byrow = T))$p.value

# Putative/STARR-seq only
BG3_all_stats[3,1] <- fisher.test(matrix(c(BG3_putative_all, BG3_starr_all),
  ncol = 2, byrow = T))$p.value

# Putative/Whole Genome
BG3_all_stats[4,1] <- fisher.test(matrix(c(BG3_putative_all, BG3_background_all),
  ncol = 2, byrow = T))$p.value

# Common/STARR-seq only
BG3_all_stats[3,2] <- fisher.test(matrix(c(BG3_common_all, BG3_starr_all),
  ncol = 2, byrow = T))$p.value

# Common/Whole Genome
BG3_all_stats[4,2] <- fisher.test(matrix(c(BG3_common_all, BG3_background_all),
  ncol = 2, byrow = T))$p.value

# STARR-seq/Whole Genome
BG3_all_stats[4,3] <- fisher.test(matrix(c(BG3_starr_all, BG3_background_all),
  ncol = 2, byrow = T))$p.value

# Writing stats matrix to a csv file
write.csv(BG3_all_stats,
  file = "~/dm6_promoter_correction/Stats/3C_stats/BG3_All.csv")


#---------------------------BG3 Stats (Distal Only)----------------------------#
#------------------------------------------------------------------------------#

# Putative
BG3_putative_d <- c(BG3_putative_distal_only$count_without_contacts,
  BG3_putative_distal_only$count_with_contacts)
# Common
BG3_common_d <- c(BG3_common_distal_only$count_without_contacts,
  BG3_common_distal_only$count_with_contacts)
# STARR-seq Only
BG3_starr_d <- c(BG3_starr_distal_only$count_without_contacts,
  BG3_starr_distal_only$count_with_contacts)
# Background
BG3_background_d <- c(BG3_background_distal_only$count_without_contacts,
  BG3_background_distal_only$count_with_contacts)

BG3_d_stats <- matrix("", 4, 4, dimnames = list(c(rownames(BG3_matrix)),
rownames(BG3_matrix)))

same <- matrix(c(1,1,2,2,3,3,4,4), 4, 2 , byrow = T)

# Filling putative/putative etc with -
BG3_d_stats[same] <- "-"

# Putative/Common
BG3_d_stats[2,1] <- fisher.test(matrix(c(BG3_putative_d, BG3_common_d),
  ncol = 2, byrow = T))$p.value

# Putative/STARR-seq only
BG3_d_stats[3,1] <- fisher.test(matrix(c(BG3_putative_d, BG3_starr_d),
  ncol = 2, byrow = T))$p.value

# Putative/Whole Genome
BG3_d_stats[4,1] <- fisher.test(matrix(c(BG3_putative_d, BG3_background_d),
  ncol = 2, byrow = T))$p.value

# Common/STARR-seq only
BG3_d_stats[3,2] <- fisher.test(matrix(c(BG3_common_d, BG3_starr_d),
  ncol = 2, byrow = T))$p.value

# Common/Whole Genome
BG3_d_stats[4,2] <- fisher.test(matrix(c(BG3_common_d, BG3_background_d),
  ncol = 2, byrow = T))$p.value

# STARR-seq/Whole Genome
BG3_d_stats[4,3] <- fisher.test(matrix(c(BG3_starr_d, BG3_background_d),
  ncol = 2, byrow = T))$p.value

# Writing stats matrix to a csv file
write.csv(BG3_d_stats,
  file = "~/dm6_promoter_correction/Stats/3C_stats/BG3_Distal_Only.csv")


#-----------------------------BG3 Stats (Proximal)-----------------------------#
#------------------------------------------------------------------------------#

# Putative
BG3_putative_pr <- c(BG3_putative_proximal$count_without_contacts,
  BG3_putative_proximal$count_with_contacts)
# Common
BG3_common_pr <- c(BG3_common_proximal$count_without_contacts,
  BG3_common_proximal$count_with_contacts)
# STARR-seq Only
BG3_starr_pr <- c(BG3_starr_proximal$count_without_contacts,
  BG3_starr_proximal$count_with_contacts)
# Background
BG3_background_pr <- c(BG3_background_proximal$count_without_contacts,
  BG3_background_proximal$count_with_contacts)

BG3_pr_stats <- matrix("", 4, 4, dimnames = list(c(rownames(BG3_matrix)),
rownames(BG3_matrix)))

same <- matrix(c(1,1,2,2,3,3,4,4), 4, 2 , byrow = T)

# Filling putative/putative etc with -
BG3_pr_stats[same] <- "-"

# Putative/Common
BG3_pr_stats[2,1] <- fisher.test(matrix(c(BG3_putative_pr, BG3_common_pr),
  ncol = 2, byrow = T))$p.value

# Putative/STARR-seq only
BG3_pr_stats[3,1] <- fisher.test(matrix(c(BG3_putative_pr, BG3_starr_pr),
  ncol = 2, byrow = T))$p.value

# Putative/Whole Genome
BG3_pr_stats[4,1] <- fisher.test(matrix(c(BG3_putative_pr, BG3_background_pr),
  ncol = 2, byrow = T))$p.value

# Common/STARR-seq only
BG3_pr_stats[3,2] <- fisher.test(matrix(c(BG3_common_pr, BG3_starr_pr),
  ncol = 2, byrow = T))$p.value

# Common/Whole Genome
BG3_pr_stats[4,2] <- fisher.test(matrix(c(BG3_common_pr, BG3_background_pr),
  ncol = 2, byrow = T))$p.value

# STARR-seq/Whole Genome
BG3_pr_stats[4,3] <- fisher.test(matrix(c(BG3_starr_pr, BG3_background_pr),
  ncol = 2, byrow = T))$p.value

# Writing stats matrix to a csv file
write.csv(BG3_pr_stats,
  file = "~/dm6_promoter_correction/Stats/3C_stats/BG3_Proximal.csv")

#
#---------------------------------S2 Stats (All)-------------------------------#
#------------------------------------------------------------------------------#
all_counts <- function(x, y, neither){
  without <- sum(x$count_without_contacts, y$count_without_contacts, neither)
  with <- sum(x$count_with_contacts, y$count_with_contacts)
  return(c(without = without, with = with))
}

# Putative
S2_putative_all <- all_counts(S2_putative_proximal, S2_putative_distal_only,
  neither = Neither_matrix[1,1])
# Common
S2_common_all <- all_counts(S2_common_proximal, S2_common_distal_only,
  neither = Neither_matrix[1,2])
# STARR-seq Only
S2_starr_all <- all_counts(S2_starr_proximal, S2_starr_distal_only,
  neither = Neither_matrix[1,3])
# Background
S2_background_all <- all_counts(S2_background_proximal, S2_background_distal_only,
  neither = Neither_matrix[1,4])

S2_all_stats <- matrix("", 4, 4, dimnames = list(c(rownames(S2_matrix)),
rownames(S2_matrix)))

same <- matrix(c(1,1,2,2,3,3,4,4), 4, 2 , byrow = T)

# Filling putative/putative etc with -
S2_all_stats[same] <- "-"

# Putative/Common
S2_all_stats[2,1] <- fisher.test(matrix(c(S2_putative_all, S2_common_all),
  ncol = 2, byrow = T))$p.value

# Putative/STARR-seq only
S2_all_stats[3,1] <- fisher.test(matrix(c(S2_putative_all, S2_starr_all),
  ncol = 2, byrow = T))$p.value

# Putative/Whole Genome
S2_all_stats[4,1] <- fisher.test(matrix(c(S2_putative_all, S2_background_all),
  ncol = 2, byrow = T))$p.value

# Common/STARR-seq only
S2_all_stats[3,2] <- fisher.test(matrix(c(S2_common_all, S2_starr_all),
  ncol = 2, byrow = T))$p.value

# Common/Whole Genome
S2_all_stats[4,2] <- fisher.test(matrix(c(S2_common_all, S2_background_all),
  ncol = 2, byrow = T))$p.value

# STARR-seq/Whole Genome
S2_all_stats[4,3] <- fisher.test(matrix(c(S2_starr_all, S2_background_all),
  ncol = 2, byrow = T))$p.value

# Writing stats matrix to a csv file
write.csv(S2_all_stats,
  file = "~/dm6_promoter_correction/Stats/3C_stats/S2_All.csv")


#---------------------------S2 Stats (Distal Only)-----------------------------#
#------------------------------------------------------------------------------#

# Putative
S2_putative_d <- c(S2_putative_distal_only$count_without_contacts,
  S2_putative_distal_only$count_with_contacts)
# Common
S2_common_d <- c(S2_common_distal_only$count_without_contacts,
  S2_common_distal_only$count_with_contacts)
# STARR-seq Only
S2_starr_d <- c(S2_starr_distal_only$count_without_contacts,
  S2_starr_distal_only$count_with_contacts)
# Background
S2_background_d <- c(S2_background_distal_only$count_without_contacts,
  S2_background_distal_only$count_with_contacts)

S2_d_stats <- matrix("", 4, 4, dimnames = list(c(rownames(S2_matrix)),
rownames(S2_matrix)))

same <- matrix(c(1,1,2,2,3,3,4,4), 4, 2 , byrow = T)

# Filling putative/putative etc with -
S2_d_stats[same] <- "-"

# Putative/Common
S2_d_stats[2,1] <- fisher.test(matrix(c(S2_putative_d, S2_common_d),
  ncol = 2, byrow = T))$p.value

# Putative/STARR-seq only
S2_d_stats[3,1] <- fisher.test(matrix(c(S2_putative_d, S2_starr_d),
  ncol = 2, byrow = T))$p.value

# Putative/Whole Genome
S2_d_stats[4,1] <- fisher.test(matrix(c(S2_putative_d, S2_background_d),
  ncol = 2, byrow = T))$p.value

# Common/STARR-seq only
S2_d_stats[3,2] <- fisher.test(matrix(c(S2_common_d, S2_starr_d),
  ncol = 2, byrow = T))$p.value

# Common/Whole Genome
S2_d_stats[4,2] <- fisher.test(matrix(c(S2_common_d, S2_background_d),
  ncol = 2, byrow = T))$p.value

# STARR-seq/Whole Genome
S2_d_stats[4,3] <- fisher.test(matrix(c(S2_starr_d, S2_background_d),
  ncol = 2, byrow = T))$p.value

# Writing stats matrix to a csv file
write.csv(S2_d_stats,
  file = "~/dm6_promoter_correction/Stats/3C_stats/S2_Distal_Only.csv")

#-----------------------------S2 Stats (Proximal)------------------------------#
#------------------------------------------------------------------------------#

S2_putative_pr <- c(S2_putative_proximal$count_without_contacts,
  S2_putative_proximal$count_with_contacts)
# Common
S2_common_pr <- c(S2_common_proximal$count_without_contacts,
  S2_common_proximal$count_with_contacts)
# STARR-seq Only
S2_starr_pr <- c(S2_starr_proximal$count_without_contacts,
  S2_starr_proximal$count_with_contacts)
# Background
S2_background_pr <- c(S2_background_proximal$count_without_contacts,
  S2_background_proximal$count_with_contacts)

S2_pr_stats <- matrix("", 4, 4, dimnames = list(c(rownames(S2_matrix)),
rownames(S2_matrix)))

same <- matrix(c(1,1,2,2,3,3,4,4), 4, 2 , byrow = T)

# Filling putative/putative etc with -
S2_pr_stats[same] <- "-"

# Putative/Common
S2_pr_stats[2,1] <- fisher.test(matrix(c(S2_putative_pr, S2_common_pr),
  ncol = 2, byrow = T))$p.value

# Putative/STARR-seq only
S2_pr_stats[3,1] <- fisher.test(matrix(c(S2_putative_pr, S2_starr_pr),
  ncol = 2, byrow = T))$p.value

# Putative/Whole Genome
S2_pr_stats[4,1] <- fisher.test(matrix(c(S2_putative_pr, S2_background_pr),
  ncol = 2, byrow = T))$p.value

# Common/STARR-seq only
S2_pr_stats[3,2] <- fisher.test(matrix(c(S2_common_pr, S2_starr_pr),
  ncol = 2, byrow = T))$p.value

# Common/Whole Genome
S2_pr_stats[4,2] <- fisher.test(matrix(c(S2_common_pr, S2_background_pr),
  ncol = 2, byrow = T))$p.value

# STARR-seq/Whole Genome
S2_pr_stats[4,3] <- fisher.test(matrix(c(S2_starr_pr, S2_background_pr),
  ncol = 2, byrow = T))$p.value

# Writing stats matrix to a csv file
write.csv(S2_pr_stats,
  file = "~/dm6_promoter_correction/Stats/3C_stats/S2_Proximal.csv")

#-------------------------------------Plotting---------------------------------#
#------------------------------------------------------------------------------#

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
  "#D55E00", "#CC79A7")
cols <- c("white", cbPalette[c(5,2,7)])


pdf(file = "~/dm6_promoter_correction/Redone_plots/Figure3C.pdf", width = 21,
  height = 6, pointsize = 14)
par(mar=c(4,3.5,4,13), mfrow=c(1,2), cex = 1.2)
par(xpd = NA)

  barplot(BG3_matrix[4:1,], beside=T, horiz = T, col = cols,
      xlim = c(0,1), xlab = "% of enhancers that contact expressed genes", xaxt = "n",
      main = "Enhancers in BG3 contact expressed genes")
    axis(1, seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10), "%"), las=2)
    legend(1.1,12.5,rownames(BG3_matrix), fill=rev(cols), bty='n')


  barplot(S2_matrix[4:1,], beside=T, horiz = T, col = cols,
          xlim = c(0,1), xlab = "% of enhancers that contact expressed genes", xaxt = "n",
          main = "Enhancers in S2 contact expressed genes")
  axis(1, seq(0, 1, 0.1), labels = paste0(seq(0, 100, 10), "%"), las=2)
  legend(1.1,12.5,rownames(BG3_matrix), fill=rev(cols), bty='n')


dev.off()
