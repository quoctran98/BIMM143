#' ---
#' title: "Class 5: R Graphics"
#' author: "Quoc Tran"
#' date: "April 16th, 2019"
#' ---

#Class 5 R graphics

#2A. Line Plot
weight <- read.table("bimm143_05_rstats/weight_chart.txt", header = TRUE)
plot(weight$Age, weight$Weight, type = "b", title = "Weight over Time", xlab = "Age (months)", ylab = "Weight (kg)")

#2B. Bar Plot
feature <- read.table("bimm143_05_rstats/feature_counts.txt", header = TRUE, sep = "\t")
par(mar = c(5,12,4,2))
barplot(feature$Count, names.arg = feature$Feature, horiz = TRUE, las = 1, title = "Feature Counts", xlab = "Count")

#3A. Providing Vector Colors
male_female <- read.table("bimm143_05_rstats/male_female_counts.txt", header = TRUE, sep = "\t")
par(mar = c(6,4,4,2))
barplot(male_female$Count, names.arg = male_female$Sample, las = 2 , col = c("blue2", "red2"), ylab = "Count")
