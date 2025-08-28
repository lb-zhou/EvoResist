#!/usr/bin/env Rscript

## This script is used to:
## Test whether prediction specificity differs significantly among Lineages 1, 2, 3, and 4.


library(data.table)
library(dplyr)

get_drug_variant_list <- function(drug_name) {
  folder_path <- paste0("/work/users/l/i/lingbo1/Mtb/WHO_denovo/results/", drug_name, "/lineage_effect")
  file_list <- list.files(folder_path, full.names = TRUE)
  #file_list <- file_list[!grepl("find_variant_ids_with_metadata.log|find_variant_ids_with_cooccurring_mutations.log", file_list)]
  file_list <- file_list[!grepl(".log$|withG1.tsv$", file_list)]
  return(file_list)
}

get_output_2 <- function(drug_name) {
  file_list <- get_drug_variant_list(drug_name)
  print(paste0("DRUG: ", drug_name))
  
  results_list <- list()
  
  for (variant in file_list) {
    input <- fread(variant, header = TRUE)
    
    # Skip empty files
    if (nrow(input) == 0) {
      next
    }
    
    input$Lineage <- ifelse(input$Lineage == "LINEAGEBOV", NA, input$Lineage)
    input <- input[input$Lineage %in% c("LINEAGE1", "LINEAGE2", "LINEAGE3", "LINEAGE4"), ]
    tab <- table(input$Phenotype, input$Lineage, useNA = "no")
    
    if (nrow(tab) >= 2 && ncol(tab) >= 2) {
      if (all(tab >= 5)) {
        ct <- chisq.test(tab)
        pval <- ct$p.value
        test_type <- "Chi-squared test"
      } else {
        ft <- fisher.test(tab)
        pval <- ft$p.value
        test_type <- "Fisher's exact test"
      }
    } else {
      next
    }
    
    if (!is.na(pval) && pval < 0.05) {
      cat(paste0("\n", variant, ":\n"))
      cat(paste0("Using ", test_type, ", p = ", round(pval, 4), "\n"))
      print(table(input$Phenotype, input$Lineage, useNA = "ifany"))
      print(table(input$Phenotype, input$Sublineage, useNA = "ifany"))
      
      spec_vec <- NA
      if (all(c("R", "S") %in% rownames(tab))) {
        R_counts <- tab["R", ]
        S_counts <- tab["S", ]
        specificity <- round(R_counts / (R_counts + S_counts), 3)
        cat("Specificity (R / total) for each lineage:\n")
        print(specificity)
        spec_vec <- paste(names(specificity), specificity, sep = ":", collapse = "; ")
      } else {
        cat("Cannot compute specificity: R or S phenotype missing in data.\n")
      }
      
      results_list[[variant]] <- data.table(
        variant = variant,
        p_value = pval,
        test_type = test_type,
        specificity = spec_vec
      )
    }
  }
  
  if (length(results_list) > 0) {
    return(rbindlist(results_list, use.names = TRUE, fill = TRUE))
  } else {
    return(data.table())  # return empty table if no significant results
  }
}

raw_results <- data.frame()
for (drug_name in c("RIF", "INH", "EMB", "PZA", "LFX", "MFX", "BDQ", "AMK", "STM", "ETO", "KAN", "CAP", "LZD")) {
  drug_df <- get_output_2(drug_name)
  raw_results <- rbind(raw_results, drug_df)
}