# =========================================================
# Script: 05_data_split.R
# Description:
# Performs a stratified 70/30 train-test split based on event status
#
# Input:
# - data_n19254.csv
#
# Output:
# - train_data_seer.csv
# - test_data_seer.csv
# =========================================================
### Data Split ###
# Load required libraries
library(survival)
library(caret)
library(dplyr)

data <- read.csv("data_n19254.csv")

# Create event indicator (1 = dead, 0 = alive)
data$event_status <- ifelse(data$Vital.status.recode == "Dead", 1, 0)

# Stratified 70-30 split by event status
set.seed(123)  # for reproducibility
train_index <- createDataPartition(data$event_status, p = 0.7, list = FALSE)

# Create training and testing datasets
train_data <- data[train_index, ]
test_data  <- data[-train_index, ]

# Display sample sizes and event proportions
cat("Train set size:", nrow(train_data), "\n")
cat("Test set size:", nrow(test_data), "\n")
cat("Event proportion in train:", round(mean(train_data$event_status), 3), "\n")
cat("Event proportion in test :", round(mean(test_data$event_status), 3), "\n")

# Save train and test sets as CSV
write.csv(train_data, "train_data_seer.csv", row.names = FALSE)
write.csv(test_data, "test_data_seer.csv", row.names = FALSE)

cat("✅ Train and test datasets saved successfully.\n")



