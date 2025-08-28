#!/bin/env Rscript

## This script is used to:
## Test if the prevance of variants within primary drug-resistant genes differ statistically significant among Lineage 1, Lineage 2, Lineage 3, and Lineage 4.

args <- commandArgs(trailingOnly = TRUE)

gene_name <- args[1]
only_nonsy <- args[2]

only_nonsy <- as.logical(only_nonsy)


library(data.table)
library(dplyr)
library(ggplot2)
library(RColorBrewer)


input_dir <- paste0("/work/users/l/i/lingbo1/Mtb/WHO_denovo/results/mutations_lineage/", gene_name, "_mutations.tsv")
input <- fread(input_dir, header = T)

## mutation info:
mutation_anno <- fread("/work/users/l/i/lingbo1/Mtb/WHO_denovo/source_data/all_ann_convergent_flt_v3.txt", header = T, fill=TRUE)
colnames(mutation_anno) <- c("gene", "codon", "anno", "position", "ref", "alt", "event_number")
mutation_anno <- mutation_anno %>% filter(!grepl("-", position))
mutation_anno$position <- as.integer(mutation_anno$position)

input_anno <- merge(input, mutation_anno, by.x=c("Position", "Ref", "Alt"), by.y = c("position", "ref", "alt"), all.x = T)


if ( (gene_name %in% c("rpoB", "embB", "katG", "pncA", "gyrA", "gyrB", "Rv0678", "pepQ", "atpE", "rpsL", "gid", "ethA", "tlyA", "rplC")) & only_nonsy) {
  
  input_anno <- input_anno %>% filter(grepl("Nonsynonymous", anno))

  dt <- as.data.table(input_anno)
  
}

if (gene_name == "inhA" & !only_nonsy) {
  input_anno <- input_anno %>%
    mutate(codon = ifelse(Position < 1674202,
                          paste0("-", 1674202 - Position),
                          codon))
  
  input_anno_used <- input_anno %>%
    filter(Position < 1674202 | (Position >= 1674202 & grepl("Nonsynonymous", anno)))
  
  
  dt <- as.data.table(input_anno_used)
}

if (gene_name == "eis" & !only_nonsy) {
  input_anno <- input_anno %>%
    mutate(codon = ifelse(Position > 2715332,
                          paste0("-", Position - 2715332),
                          codon))
  
  input_anno_used <- input_anno %>%
    filter(Position > 2715332 | (Position <= 2715332 & grepl("Nonsynonymous", anno)))
  
  
  dt <- as.data.table(input_anno_used)
}

if (gene_name == "rrs" & !only_nonsy) {
  input_anno <- input_anno %>%
    mutate(codon = Position)
  
  dt <- as.data.table(input_anno)
}

if (gene_name == "rrl" & !only_nonsy) {
  input_anno <- input_anno %>%
    mutate(codon = Position)
  
  dt <- as.data.table(input_anno)
}


top_pooled <- dt[, .N, by = codon][order(-N)][1:30]
top_codon_set <- top_pooled$codon


top_all <- dt[codon %in% top_codon_set, .N, by = codon][
  , total := sum(N)][
    , proportion := N / total][
      , Lineage := "Pooled"][
        , .(Lineage, codon, proportion)]


top_by_lineage <- dt[codon %in% top_codon_set, .N, by = .(Lineage, codon)][
  , total := sum(N), by = Lineage][
    , proportion := N / total][
      , .(Lineage, codon, proportion)]

top_combined <- rbind(top_all, top_by_lineage)


heat_data <- dcast(top_combined, codon ~ Lineage, value.var = "proportion", fill = 0)
heat_long <- melt(heat_data, id.vars = "codon", variable.name = "Lineage", value.name = "proportion")


codon_order <- top_pooled$codon

heat_long[, Lineage := fifelse(Lineage == "Pooled", "Pooled", gsub("LINEAGE", "L", Lineage))]


dt_binary <- dt[codon %in% codon_order, .(ID, Lineage, codon)]
n_total_lineage <- unique(dt_binary[, .(ID, Lineage)])[, .N, by = Lineage]

test_results <- list()
for (codon_i in codon_order) {
  has_mut <- unique(dt_binary[codon == codon_i, .(ID, Lineage)])
  n_present <- has_mut[, .N, by = Lineage]
  
  codon_tab <- merge(n_total_lineage, n_present, by = "Lineage", all.x = TRUE)
  setnames(codon_tab, c("Lineage", "n_total", "n_present"))
  codon_tab[is.na(n_present), n_present := 0]
  codon_tab[, n_absent := n_total - n_present]
  
  mat <- as.matrix(codon_tab[, .(n_present, n_absent)])
  rownames(mat) <- codon_tab$Lineage
  mat <- t(mat)
  
  test_result <- tryCatch({
    temp_test <- suppressWarnings(chisq.test(mat))  # suppress warning for now
    if (all(temp_test$expected >= 5)) {
      list(test = temp_test, type = "Chi-sq")
    } else {
      NULL
    }
  }, error = function(e) {
    message(sprintf("Skipping test for %s due to error: %s", codon_i, e$message))
    NULL
  })
  
  
  if (is.null(test_result)) {

    test_results[[codon_i]] <- data.table(
      codon = codon_i,
      p_value = NA_real_,
      test_type = NA_character_
    )
  } else {
    test_results[[codon_i]] <- data.table(
      codon = codon_i,
      p_value = test_result$test$p.value,
      test_type = test_result$type
    )
  }

}

codon_test_results <- rbindlist(test_results, use.names = TRUE, fill = TRUE)


codon_test_results[, p_adj := NA_real_]
codon_test_results[!is.na(p_value),
                   p_adj := p.adjust(p_value[!is.na(p_value)], method = "BH")]

codon_test_results[, star := ifelse(is.na(p_value), "-",
                                    ifelse(p_value < 1e-4, "***",
                                           ifelse(p_value < 1e-3, "**",
                                                  ifelse(p_value < 5e-2, "*", ""))))]


heat_long <- merge(heat_long, codon_test_results[, .(codon, star)], by = "codon", all.x = TRUE)
heat_long[, codon := factor(codon, levels = codon_order)]
heat_long[, Lineage := as.factor(Lineage)]

p <- ggplot(heat_long, aes(x = Lineage, y = codon, fill = proportion)) +
  geom_tile(color = "white", linewidth = 0.2) +
  # heatmap cells
  scale_fill_distiller(palette = "Reds", direction = 1, limits = c(0,1)) +
  labs(x = NULL, y = NULL, fill = "Proportion") +
  theme_classic() +
  geom_text(
    data = heat_long[Lineage == "Pooled" & star != ""],
    aes(x = length(levels(heat_long$Lineage)) + 0.5,
        y = codon, label = star),
    inherit.aes = FALSE, hjust = 0, vjust= 0.8, size = 8
  ) + 
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 17, angle = -90, hjust = 0.5),
    axis.text.y = element_text(size = 18),
    axis.title.x = element_text(size = 18),
    axis.title.y = element_text(size = 18),
    plot.margin = margin(5, 30, 5, 5)
  ) + 
  coord_cartesian(clip = "off") 