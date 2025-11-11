#!/usr/bin/env Rscript
# -------------------------------------------------------------------------
# Script: murine-mirna-target-dotplot.R
# Author: Bicna Song | Nov 2025
# Description:
#   Generate dot plots of miRNA–immune target interactions for a specified
#   cluster table provided as a command-line argument.
#   Example:
#     Rscript murine-mirna-target-dotplot.R \
#       results/cluster1-mirna-immune-target-interactions.tsv
# -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(rprojroot)
  library(tidyverse)
  library(ggplot2)
  library(forcats)
})

#  Parse argument and derive cluster ID
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Usage: Rscript murine-mirna-target-dotplot.R <path_to_interaction_table.tsv>")
}
infile <- args[1]
cluster_id <- str_extract(basename(infile), "cluster\\d+")

#  Paths and theme
root_dir <- find_root(has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "murine-mirna-target-interactions")
plot_dir <- file.path(analysis_dir, "plots")
source(file.path(root_dir, "figures", "theme.R"))

#  Load and classify data
df <- read_tsv(infile, show_col_types = FALSE)

df <- df %>%
  mutate(category = case_when(
    str_detect(Term, regex("T cell|CD8|alpha-beta", ignore_case = TRUE)) ~ "T cell-related pathways",
    str_detect(Term, regex("myeloid|macrophage", ignore_case = TRUE))     ~ "Myeloid-related pathways",
    str_detect(Term, regex("antigen receptor", ignore_case = TRUE))       ~ "Antigen Presentation",
    str_detect(Term, regex("interferon-beta", ignore_case = TRUE))        ~ "Interferon-beta Production",
    TRUE                                                                  ~ "Other immune pathways"
  ))

#  Helper to build and save plots
make_plot <- function(data, title_suffix, outfile, color_values, height = 18) {
  p <- ggplot(data, aes(x = miRNA, y = Genes)) +
    geom_point(aes(size = abs(log2FoldChange_b7h3_vs_untr_Day21),
                   color = expr_pattern),
               alpha = 0.85) +
    scale_color_manual(values = color_values) +
    scale_size_continuous(name = "|log2FC|", range = c(2, 7)) +
    labs(
      x = "miRNA",
      y = "Immune-related target gene",
      color = "Expression pattern",
      title = paste("miRNA–Immune Target Interactions", title_suffix)
    ) +
    facet_grid(category ~ ., scales = "free_y", space = "free_y", switch = "y") +
    theme_Publication() +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
  
  ggsave(file.path(plot_dir, outfile), p, width = 12, height = height)
}

#  Prepare reordered data
plot_df_full <- df %>%
  mutate(
    Genes = fct_reorder(Genes, log2FoldChange_b7h3_vs_untr_Day21, .fun = median),
    miRNA = fct_reorder(miRNA, log2FoldChange_b7h3_vs_untr_Day21, .fun = median)
  )

# Full plot (all expression patterns)
make_plot(
  plot_df_full,
  sprintf("(Full – %s)", cluster_id),
  sprintf("%s-immune-mirna-target-dotplot-full.pdf", cluster_id),
  color_values = c(
    "B7H3 down, Day 14" = "#F8766D",
    "B7H3 down, Day 21" = "#B22222",
    "B7H3 up, Day 14"   = "#7FD67F",
    "B7H3 up, Day 21"   = "#228B22",
    "B7H3 up, multiple" = "#006400"
  ),
  height = 18
)

# Filtered plot (downregulated targets only)
plot_df_filtered <- df %>%
  filter(expr_pattern %in% c("B7H3 down, Day 14", "B7H3 down, Day 21")) %>%
  mutate(
    Genes = fct_reorder(Genes, log2FoldChange_b7h3_vs_untr_Day21, .fun = median),
    miRNA = fct_reorder(miRNA, log2FoldChange_b7h3_vs_untr_Day21, .fun = median)
  )

make_plot(
  plot_df_filtered,
  sprintf("(Downregulated – %s)", cluster_id),
  sprintf("%s-immune-mirna-target-dotplot-filtered.pdf", cluster_id),
  color_values = c(
    "B7H3 down, Day 14" = "#F8766D",
    "B7H3 down, Day 21" = "#B22222"
  ),
  height = 10
)

message(sprintf("[%s] Plots for %s completed successfully.", Sys.time(), cluster_id))
