# Importing libraries
library(BSgenome.Dmelanogaster.UCSC.dm6)
library(genomation)
library(GenomicRanges)
library(rtracklayer)
predicted <- get(load("/home/jw18713/Project1/data/model_predictions/BG3_FuzzyGR.Rda"))

load("loop_data100.Rda")
load("grangelist_BG3.Rda")

.sumReadsConfPerc1 <- function(methylationData){   #Function to compute sums within below function
  return(sum(methylationData$'Conf..Perc.1'))
}
.ConfPercAvg <- function(methylationData, regions){

  overlaps <- findOverlaps(methylationData, regions, ignore.strand = TRUE)
  methylationDataContextList <- IRanges::splitAsList(methylationData[queryHits(overlaps)],  subjectHits(overlaps))
  regionsIndexes <- as.integer(names(methylationDataContextList))

  regions$sumConfPerc1 <- rep(0, times=length(regions)) #Confidence percentage 1

  if(length(regionsIndexes) > 0){
    regions$sumConfPerc1[regionsIndexes] <- sapply(methylationDataContextList,.sumReadsConfPerc1) #Confidence percentage 1
  }

  regions$avgConfPerc1 <- regions$sumConfPerc1 / (width(regions)/10)

  return(regions)
}

one <- glist[loop_data == 1]
two <- glist[loop_data == 2]
three <- glist[loop_data == 3]
more <- glist[loop_data > 3]

# Merging and combining regions with 2 bins to test

two_check <- two
for (i in seq_along(two_check)){
  two_check[[i]] <- reduce(two_check[[i]], min.gapwidth = 101)
}

two_check <- unlist(two_check)
two_merged <- .ConfPercAvg(predicted, two_check)
merge_over80 <- two_merged[two_merged$avgConfPerc1 >= 0.8]

to_keep <- reduce(c(unlist(one), unlist(two), merge_over80))

# Merging and combining all 3 bins where 3 can be tested

three_all <- three
for (i in seq_along(three_all)){
  three_all[[i]] <- reduce(three_all[[i]], min.gapwidth = 101)
}

three_all <- unlist(three_all)
three_all_conf <- .ConfPercAvg(predicted, three_all)
three_all_over80 <- three_all_conf[three_all_conf$avgConfPerc1 >= 0.8]
three_all_under80 <- three_all_conf[three_all_conf$avgConfPerc1 < 0.8]

to_keep <- c(to_keep, three_all_over80)
three_check_left <- three[three %over% three_all_under80]

# Now doing the left side of regions with 3 bins that were below threshold

three_left <- three_check_left
for (i in seq_along(three_left)){
  three_left[[i]] <- reduce(three_left[[i]][1:2], min.gapwidth = 101)
}

three_left <- unlist(three_left)
three_left_conf <- .ConfPercAvg(predicted, three_left)
three_left_over80 <- three_left_conf[three_left_conf$avgConfPerc1 >= 0.8]
three_left_under80 <- three_left_conf[three_left_conf$avgConfPerc1 < 0.8]

to_keep <- c(to_keep, three_left_over80)
three_check_right <- three[three %over% three_left_under80]

# Now doing the right side of regions with 3 bins that were below threshold

three_right <- three_check_right
for (i in seq_along(three_right)){
  three_right[[i]] <- reduce(three_right[[i]][2:3], min.gapwidth = 101)
}

three_right <- unlist(three_right)
three_right_conf <- .ConfPercAvg(predicted, three_right)
three_right_over80 <- three_right_conf[three_right_conf$avgConfPerc1 >= 0.8]
three_right_under80 <- three_right_conf[three_right_conf$avgConfPerc1 < 0.8]

to_keep <- c(to_keep, three_right_over80, unlist(three))
to_keep <- reduce(to_keep)

# Ok, this is everything up to and including 3. From here it gets harder.

func_test <- function(glist, predicted){
  results <- glist
  for (i in seq_along(glist)){
    print (i)
    reduced <- .ConfPercAvg(predicted, reduce(glist[[i]], min.gapwidth = 101))
    if (reduced$avgConfPerc1 >= 0.8){
      results[[i]] <- reduced
    }
  }
  return(results)
}

more_tested <- func_test(more, predicted)
length_vector <- rep(0, length(more_tested))

for (i in seq_along(length_vector)){
  length_vector[i] <- length(more_tested[[i]])
}

to_keep <- reduce(c(to_keep, unlist(more_tested[length_vector == 1])))
awkward <- more_tested[!length_vector == 1]


# Ok, that got everything that combines easily. Now for the awkward ones.

awkward_count <- rep(0, length(awkward))

for (i in seq_along(awkward_count)){
  awkward_count[i] <- length(awkward[[i]])
}

four <- awkward[awkward_count == 4] # 21 of these
five <- awkward[awkward_count == 5] # 3 of these
six <- awkward[awkward_count == 6] # 4 of these

# Alright, it's time to do four first
# We know these don't combine in total, so lets do 1:3 (left), and 2:4 (right)

four_left <- four

for (i in seq_along(four_left)){
  four_left[[i]] <- four_left[[i]][1:3]
}

four_right <- four

for (i in seq_along(four_right)){
  four_right[[i]] <- four_right[[i]][2:4]
}

# Now that I have left and right split it's time to make them converge
four_left <- reduce(unlist(four_left), min.gapwidth = 101)
four_left_conf <- .ConfPercAvg(predicted, four_left)
four_l1_over80 <- four_left_conf[four_left_conf$avgConfPerc1 >= 0.8] # This gives me 7
four_l1_under80 <- four_left_conf[four_left_conf$avgConfPerc1 < 0.8]

four_right <- reduce(unlist(four_right), min.gapwidth = 101)
four_right_conf <- .ConfPercAvg(predicted, four_right)
four_r1_over80 <- four_right_conf[four_right_conf$avgConfPerc1 >= 0.8] # This gives me 3
four_r1_under80 <- four_right_conf[four_right_conf$avgConfPerc1 < 0.8]

to_keep <- c(to_keep, four_l1_over80, four_r1_over80) # Here I've added those 10
# to the to keep pile. This is the longest either of these will be.

# Now for some smaller fragments
four_left2 <- four[four %over% four_l1_under80]
four_right2 <- four[four %over% four_r1_under80]

for (i in seq_along(four_left2)){
  four_left2[[i]] <- four_left2[[i]][1:2]
}

for (i in seq_along(four_right2)){
  four_right2[[i]] <- four_right2[[i]][3:4]
}

four_left2 <- reduce(unlist(four_left2), min.gapwidth = 101)
four_left_conf2 <- .ConfPercAvg(predicted, four_left2)
four_l2_over80 <- four_left_conf2[four_left_conf2$avgConfPerc1 >= 0.8] # This gives me 8

four_right2 <- reduce(unlist(four_right2), min.gapwidth = 101)
four_right_conf2 <- .ConfPercAvg(predicted, four_right2)
four_r2_over80 <- four_right_conf2[four_right_conf2$avgConfPerc1 >= 0.8] # This gives me 7

to_keep <- reduce(c(to_keep, four_l2_over80, four_r2_over80, unlist(four)))

# What combinations can we have?
# 1:4 c1
# 2:5 c2 1 result here
# 1:3 l1 1 result here
# 3:5 r1
# 1:2 s1
# 2:3 s2
# 3:4 s3
# 4:5 s4

c1 <- five
c2 <- five

for (i in seq_along(c1)){
  c1[[i]] <- c1[[i]][1:4]
}

for (i in seq_along(c2)){
  c2[[i]] <- c2[[i]][2:5]
}

# Now that I have center 1 and center 2 split it's time to make them converge
c1 <- reduce(unlist(c1), min.gapwidth = 101)
c1_conf <- .ConfPercAvg(predicted, c1)
c1_over80 <- c1_conf[c1_conf$avgConfPerc1 >= 0.8]
c1_under80 <- c1_conf[c1_conf$avgConfPerc1 < 0.8]

c2 <- reduce(unlist(c2), min.gapwidth = 101)
c2_conf <- .ConfPercAvg(predicted, c2)
c2_over80 <- c2_conf[c2_conf$avgConfPerc1 >= 0.8]
c2_under80 <- c2_conf[c2_conf$avgConfPerc1 < 0.8]

to_keep <- c(to_keep, c1_over80, c2_over80)

# Ok, now we have to do the sets of 3 bins, there are currently two left

l1 <- five[five %over% c1_under80]
r1 <- five[five %over% c2_under80]

for (i in seq_along(l1)){
  l1[[i]] <- l1[[i]][1:3]
}

for (i in seq_along(r1)){
  r1[[i]] <- r1[[i]][3:5]
}

# Now that I have left 1 and right 1 split it's time to make them converge
l1 <- reduce(unlist(l1), min.gapwidth = 101)
l1_conf <- .ConfPercAvg(predicted, l1)
l1_over80 <- l1_conf[l1_conf$avgConfPerc1 >= 0.8]
l1_under80 <- l1_conf[l1_conf$avgConfPerc1 < 0.8]

r1 <- reduce(unlist(r1), min.gapwidth = 101)
r1_conf <- .ConfPercAvg(predicted, r1)
r1_over80 <- r1_conf[r1_conf$avgConfPerc1 >= 0.8]
r1_under80 <- r1_conf[r1_conf$avgConfPerc1 < 0.8]

to_keep <- c(to_keep, l1_over80, r1_over80)

# Ok, now I have pairs of two, and I'm looking for one more thing

s1 <- five
s2 <- five
s3 <- five
s4 <- five

for (i in seq_along(s1)){
  s1[[i]] <- s1[[i]][1:2]
}

for (i in seq_along(s2)){
  s2[[i]] <- s2[[i]][2:3]
}

for (i in seq_along(s3)){
  s3[[i]] <- s3[[i]][3:4]
}

for (i in seq_along(s4)){
  s4[[i]] <- s4[[i]][4:5]
}

# And on to looking for options

s1 <- reduce(unlist(s1), min.gapwidth = 101)
s1_conf <- .ConfPercAvg(predicted, s1)
s1_over80 <- s1_conf[s1_conf$avgConfPerc1 >= 0.8]
s1_under80 <- s1_conf[s1_conf$avgConfPerc1 < 0.8]

s2 <- reduce(unlist(s2), min.gapwidth = 101)
s2_conf <- .ConfPercAvg(predicted, s2)
s2_over80 <- s2_conf[s2_conf$avgConfPerc1 >= 0.8]
s2_under80 <- s2_conf[s2_conf$avgConfPerc1 < 0.8]

s3 <- reduce(unlist(s3), min.gapwidth = 101)
s3_conf <- .ConfPercAvg(predicted, s3)
s3_over80 <- s3_conf[s3_conf$avgConfPerc1 >= 0.8]
s3_under80 <- s3_conf[s3_conf$avgConfPerc1 < 0.8]

s4 <- reduce(unlist(s4), min.gapwidth = 101)
s4_conf <- .ConfPercAvg(predicted, s4)
s4_over80 <- s4_conf[s4_conf$avgConfPerc1 >= 0.8]
s4_under80 <- s4_conf[s4_conf$avgConfPerc1 < 0.8]

# Finally, getting some more to keep

to_keep <- reduce(c(to_keep, s1, s2, s3, s4, unlist(five)))

# Ok, on to the sixes

# So what divisions do we have?
# 1:5 sixes1
# 2:6 sixes2
# 1:4 sixes3
# 3:6 sixes4
# 1:3 sixes5
# 4:6 sixes6
# 1:2 sixes7
# 2:3 sixes8
# 3:4 sixes9
# 4:5 sixes10
# 5:6 sixes11

sixes1 <- six
sixes2 <- six
sixes3 <- six
sixes4 <- six
sixes5 <- six
sixes6 <- six
sixes7 <- six
sixes8 <- six
sixes9 <- six
sixes10 <- six
sixes11 <- six


for (i in seq_along(sixes1)){
  sixes1[[i]] <- sixes1[[i]][1:2]
}


for (i in seq_along(sixes2)){
  sixes2[[i]] <- sixes2[[i]][1:2]
}


for (i in seq_along(sixes3)){
  sixes3[[i]] <- sixes3[[i]][1:2]
}


for (i in seq_along(sixes4)){
  sixes4[[i]] <- sixes4[[i]][1:2]
}


for (i in seq_along(sixes5)){
  sixes5[[i]] <- sixes5[[i]][1:2]
}


for (i in seq_along(sixes6)){
  sixes6[[i]] <- sixes6[[i]][1:2]
}


for (i in seq_along(sixes7)){
  sixes7[[i]] <- sixes7[[i]][1:2]
}


for (i in seq_along(sixes8)){
  sixes8[[i]] <- sixes8[[i]][1:2]
}


for (i in seq_along(sixes9)){
  sixes9[[i]] <- sixes9[[i]][1:2]
}


for (i in seq_along(sixes10)){
  sixes10[[i]] <- sixes10[[i]][1:2]
}


for (i in seq_along(sixes11)){
  sixes11[[i]] <- sixes11[[i]][1:2]
}

# And now part 2

sixes1 <- reduce(unlist(sixes1), min.gapwidth = 101)
sixes1_conf <- .ConfPercAvg(predicted, sixes1)
sixes1_over80 <- sixes1_conf[sixes1_conf$avgConfPerc1 >= 0.8]
sixes1_under80 <- sixes1_conf[sixes1_conf$avgConfPerc1 < 0.8]

sixes2 <- reduce(unlist(sixes2), min.gapwidth = 101)
sixes2_conf <- .ConfPercAvg(predicted, sixes2)
sixes2_over80 <- sixes2_conf[sixes2_conf$avgConfPerc1 >= 0.8]
sixes2_under80 <- sixes2_conf[sixes2_conf$avgConfPerc1 < 0.8]

sixes3 <- reduce(unlist(sixes3), min.gapwidth = 101)
sixes3_conf <- .ConfPercAvg(predicted, sixes3)
sixes3_over80 <- sixes3_conf[sixes3_conf$avgConfPerc1 >= 0.8]
sixes3_under80 <- sixes3_conf[sixes3_conf$avgConfPerc1 < 0.8]

sixes4 <- reduce(unlist(sixes4), min.gapwidth = 101)
sixes4_conf <- .ConfPercAvg(predicted, sixes4)
sixes4_over80 <- sixes4_conf[sixes4_conf$avgConfPerc1 >= 0.8]
sixes4_under80 <- sixes4_conf[sixes4_conf$avgConfPerc1 < 0.8]

sixes5 <- reduce(unlist(sixes5), min.gapwidth = 101)
sixes5_conf <- .ConfPercAvg(predicted, sixes5)
sixes5_over80 <- sixes5_conf[sixes5_conf$avgConfPerc1 >= 0.8]
sixes5_under80 <- sixes5_conf[sixes5_conf$avgConfPerc1 < 0.8]

sixes6 <- reduce(unlist(sixes6), min.gapwidth = 101)
sixes6_conf <- .ConfPercAvg(predicted, sixes6)
sixes6_over80 <- sixes6_conf[sixes6_conf$avgConfPerc1 >= 0.8]
sixes6_under80 <- sixes6_conf[sixes6_conf$avgConfPerc1 < 0.8]

sixes7 <- reduce(unlist(sixes7), min.gapwidth = 101)
sixes7_conf <- .ConfPercAvg(predicted, sixes7)
sixes7_over80 <- sixes7_conf[sixes7_conf$avgConfPerc1 >= 0.8]
sixes7_under80 <- sixes7_conf[sixes7_conf$avgConfPerc1 < 0.8]

sixes8 <- reduce(unlist(sixes8), min.gapwidth = 101)
sixes8_conf <- .ConfPercAvg(predicted, sixes8)
sixes8_over80 <- sixes8_conf[sixes8_conf$avgConfPerc1 >= 0.8]
sixes8_under80 <- sixes8_conf[sixes8_conf$avgConfPerc1 < 0.8]

sixes9 <- reduce(unlist(sixes9), min.gapwidth = 101)
sixes9_conf <- .ConfPercAvg(predicted, sixes9)
sixes9_over80 <- sixes9_conf[sixes9_conf$avgConfPerc1 >= 0.8]
sixes9_under80 <- sixes9_conf[sixes9_conf$avgConfPerc1 < 0.8]

sixes10 <- reduce(unlist(sixes10), min.gapwidth = 101)
sixes10_conf <- .ConfPercAvg(predicted, sixes10)
sixes10_over80 <- sixes10_conf[sixes10_conf$avgConfPerc1 >= 0.8]
sixes10_under80 <- sixes10_conf[sixes10_conf$avgConfPerc1 < 0.8]

sixes11 <- reduce(unlist(sixes11), min.gapwidth = 101)
sixes11_conf <- .ConfPercAvg(predicted, sixes11)
sixes11_over80 <- sixes11_conf[sixes11_conf$avgConfPerc1 >= 0.8]
sixes11_under80 <- sixes11_conf[sixes11_conf$avgConfPerc1 < 0.8]

to_keep <- reduce(c(to_keep,
  sixes1_over80,
  sixes2_over80,
  sixes3_over80,
  sixes4_over80,
  sixes5_over80,
  sixes6_over80,
  sixes7_over80,
  sixes8_over80,
  sixes9_over80,
  sixes10_over80,
  sixes11_over80,
  unlist(six)))

pred_enh <- to_keep
save(pred_enh, file = "grown_predicted_enhancers_BG3.Rda")
