#!/usr/bin/env Rscript
# -------------------------------------------------------------------------
# Script: 01-murine-mirna-target-interactions.R
# Author: Bicna Song
# Date: 2025-11
# Description:
#   Assess coordinated differential expression between murine miRNAs and
#   their predicted immune-related target genes across clusters.
#   This script integrates miRNA target predictions (miRanda), immune-related
#   GO enrichment results, and RNA-seq differential expression data to identify
#   significant miRNA–target interactions for downstream visualization.
# -------------------------------------------------------------------------

### Setup
suppressPackageStartupMessages({
  library(rprojroot)
  library(tidyverse)
  library(ggplot2)
  library(gtools)
})

### Paths
root_dir <- find_root(has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "murine-mirna-target-interactions")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

source(file.path(root_dir, "figures", "theme.R"))

### Input files
miranda_output_file <- file.path(
  root_dir, "analyses", "murine-novel-mirna-target-prediction",
  "results", "mouse_miranda_output_parsed_anno.csv"
)

de_rna_file <- file.path(
  root_dir, "analyses", "murine-rna-expression",
  "results", "rna-differential-expression-deseq2-b7h3-stop-vs-untreated-by-timepoint.tsv"
)

# Load shared data

### Load miRNA–target predictions (miRanda)
miranda_df <- read_csv(miranda_output_file, show_col_types = FALSE) %>%
  filter(Max_Score >= 145, Max_Energy <= -20) %>%
  distinct(Seq1, GeneSymbol, .keep_all = TRUE)

### Load differential expression results
de_rna_df <- read_tsv(de_rna_file, show_col_types = FALSE) %>%
  filter(!expr_pattern %in% c("Not DE", "Other DE"))

# Function to process each cluster
process_cluster <- function(cluster_id) {
  message("Processing cluster ", cluster_id, " ...")
  
  # Load immune-related GO enrichment results for this cluster
  cluster_file <- file.path(
    root_dir, "analyses", "murine-mirna-go-enrichment",
    "results", sprintf("cluster%s-mirna-target-go-enr-immune-terms-reduced.tsv", cluster_id)
  )
  
  if (!file.exists(cluster_file)) {
    warning("File not found for cluster ", cluster_id)
    return(NULL)
  }
  
  cluster_go_df <- read_tsv(cluster_file, show_col_types = FALSE) %>%
    filter(!is.na(Term)) %>%
    select(miRNA, Term, Genes) %>%
    separate_rows(Genes, sep = ",\\s*") %>%
    group_by(miRNA, Genes) %>%
    summarize(Term = str_c(Term, collapse = ";"), .groups = "drop")
  
  # Join GO results with miRanda predictions
  cluster_go_df <- cluster_go_df %>%
    left_join(miranda_df, by = c("miRNA" = "Seq1", "Genes" = "GeneSymbol")) %>%
    filter(!is.na(GeneID))
  
  # Join with RNA DE results
  cluster_go_df <- cluster_go_df %>%
    left_join(de_rna_df, by = c("Genes" = "gene_symbol")) %>%
    filter(!is.na(expr_pattern))
  
  # Save outputs
  write_tsv(cluster_go_df,
            file.path(results_dir, sprintf("cluster%s-mirna-immune-target-interactions.tsv", cluster_id)))

  message("✓ Saved results for cluster ", cluster_id)
  return(cluster_go_df)
}

# Run for all clusters
target_clusters <- c(1, 2, 3, 4, 5)
cluster_results <- map(target_clusters, process_cluster)

# Session info
sessionInfo()

