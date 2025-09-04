#!/usr/bin/env Rscript

## This script is used to:
## Select the optimized mutational event threshold.

library(data.table)
library(dplyr)
library(ggpubr)
library(pROC)

all_set <- fread("/work/users/l/i/lingbo1/Mtb/WHO_denovo/ms_table_figures/roc_used_mutations.tsv", header = T)


all_set[, label_bin := ifelse(true_label == "positive", 1, 0)]

cutoffs <- seq(min(all_set$count), max(all_set$count))

roc_metrics <- lapply(cutoffs, function(cut) {
  all_set[, pred := ifelse(count >= cut, 1, 0)]
  TP <- sum(all_set$pred == 1 & all_set$label_bin == 1)
  FP <- sum(all_set$pred == 1 & all_set$label_bin == 0)
  FN <- sum(all_set$pred == 0 & all_set$label_bin == 1)
  TN <- sum(all_set$pred == 0 & all_set$label_bin == 0)
  
  TPR <- TP / (TP + FN)
  FPR <- FP / (FP + TN)
  specificity <- TN / (TN + FP)
  youden_J <- TPR + specificity - 1
  
  return(data.frame(
    cutoff = cut,
    TPR = TPR,
    FPR = FPR,
    TP = TP,
    FP = FP,
    TN = TN,
    FN = FN,
    youden_J = youden_J
  ))
})

roc_df <- do.call(rbind, roc_metrics)

best_row <- roc_df[which.max(roc_df$youden_J), ]

roc_df$text_label <- ifelse(roc_df$cutoff <=20, roc_df$cutoff, NA)

roc <- ggplot(roc_df, aes(x = FPR, y = TPR)) +
  # geom_path(color = "#f38400") +
  # geom_point(color = "#f38400", size=2.5) +
  geom_path(color = "#D62728") +
  geom_point(color = "#D62728", size=2.5) +
  # geom_text(aes(label = text_label), vjust = -0.5, size = 4, color = "black") +
  geom_abline(slope = 1, intercept = 0, color="#00468bff", linetype=2) + 
  geom_point(data = best_row, aes(x = FPR, y = TPR), color = "#00468bff", size = 3.5, shape = 24, fill  = "#00468bff") +
  # geom_text(data = best_row, aes(label = paste0("Best: ", cutoff)), 
  #           vjust = 1.5, color = "darkgreen", size = 4) +
  labs(title = NULL,
       x = "False Positive Rate",
       y = "True Positive Rate") +
  coord_fixed() +  # ensures 1:1 aspect ratio
  theme_classic() + 
  theme(axis.text.x = element_text(size = 17, hjust = 1),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_text(size = 18, face = "bold"),
        axis.title.y = element_text(size = 18, face = "bold"))

