# Assess DMG cluster 6 miRNA - neuronal target interactions
# Ryan Corbett | Oct 2025
# Updated by Bicna Song | Jan 2026 

# This script assesses correlations between DMG cluster 6 miRNA and neuronal pathway target expression, and assesses coordinated DE of targets to identify significant associations

# Load packages
library(tidyverse)
library(ggplot2)
library(gtools)

# set dir paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "human-neuronal-expression")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

source(file.path(root_dir, "figures", "theme.R"))

# set file paths
dmg_cluster6_go_file <- file.path(results_dir,
                                  "DIPG_or_DMG-cluster6-mirna-target-go-enr-neuronal-terms-reduced.tsv")

mirtarbase_file <- file.path(root_dir, "analyses",
                             "human-mirna-go-enrichment",
                             "results",
                             "human-mirna-target-predictions-miRTarBase.tsv.gz")

hist_file <- file.path(root_dir, "analyses",
                       "histology-preprocessing",
                       "results",
                       "histologies.tsv")

mirna_tpm_file <- file.path(root_dir,
                            "analyses",
                            "immune-deconvolution",
                            "results",
                             "mirna-tpm.rds")

target_tpm_file <- file.path(data_dir,
                      "gene-expression-rsem-tpm-collapsed.all.rds")

de_mirna_healthy_file <- file.path(root_dir, 
                                   "analyses",
                                   "human-rna-expression",
                                   "results",
                                   "DESeq2_DIPG or DMG_vs_healthyNormal.csv")

de_mirna_adjacent_file <- file.path(root_dir, 
                                   "analyses",
                                   "human-rna-expression",
                                   "results",
                                   "DESeq2_DIPG or DMG_paired_full.csv")

# Wrangle data

dmg_cluster6_go_df <- read_tsv(dmg_cluster6_go_file) %>%
  dplyr::filter(!is.na(Term)) %>%
  dplyr::select(miRNA, Term, Genes) %>%
  separate_rows(Genes, sep = ",\\s*") %>%
  group_by(miRNA, Genes) %>%
  summarise(Term = str_c(Term, collapse = ";"))

mirtarbase_df <- read_tsv(mirtarbase_file) %>%
  dplyr::filter(`Species (miRNA)` == "hsa") %>% 
  distinct(miRNA, `Target Gene`, .keep_all = TRUE) %>%
  dplyr::select(miRNA, `Target Gene`,
                Experiments, `Support Type`)

dmg_cluster6_go_df <- dmg_cluster6_go_df %>%
  left_join(mirtarbase_df,
            by = c("miRNA",
                   "Genes" = "Target Gene")) %>%
  dplyr::filter(!is.na(`Support Type`))


# Load hist and extract DMG miRNA and RNA IDs

hist <- read_tsv(hist_file)

dmg_mirna_ids <- hist %>% 
  dplyr::filter(experimental_strategy == "miRNA-Seq",
                sample_type == "Tumor",
                histology == "DIPG or DMG",
                !external_sample_id %in% c("1-1855-CC1", 
                                           "4-1958-CC1")) %>%
  arrange(external_sample_id) %>%
  pull(external_sample_id)

dmg_rna_ids <- hist %>% 
  dplyr::filter(experimental_strategy == "RNA-Seq",
                sample_type == "Tumor",
                histology == "DIPG or DMG") %>%
  arrange(external_sample_id) %>%
  pull(Bioassay_ID)

# load TPM files 

mirna_tpm <- readRDS(mirna_tpm_file) %>%
  dplyr::select(any_of(c(dmg_mirna_ids)))

target_tpm <- readRDS(target_tpm_file) %>%
  dplyr::select(any_of(c(dmg_rna_ids)))


# add empty pearson correlation coefficient and p-value columns to cluster 6 df

dmg_cluster6_go_df <- dmg_cluster6_go_df %>%
  dplyr::mutate(mirna_target_pearson_r = NA_integer_,
                mirna_target_pearson_p = NA_integer_)

for (i in 1:nrow(dmg_cluster6_go_df)){
  
  miRNA <- dmg_cluster6_go_df$miRNA[i]
  target <- dmg_cluster6_go_df$Genes[i]
  
  if (!target %in% rownames(target_tpm)) next
  
  miRNA_tpms <- unlist(mirna_tpm[miRNA, dmg_mirna_ids])
  target_tpms <- unlist(target_tpm[target, dmg_rna_ids])
  
  dmg_cluster6_go_df$mirna_target_pearson_r[i] <- cor.test(log2(miRNA_tpms + 0.01),
                                                           log2(target_tpms + 0.01),
                                                           method = "pearson")$estimate
  
  dmg_cluster6_go_df$mirna_target_pearson_p[i] <- cor.test(log2(miRNA_tpms + 0.01),
                                                           log2(target_tpms + 0.01),
                                                           method = "pearson")$p.value
  
}

# Load DE results

de_mirna_healthy_df <- read_csv(de_mirna_healthy_file) %>%
  dplyr::select(gene_symbol, log2FoldChange,
                padj) %>%
  dplyr::rename(dmg_vs_healthy_log2FC = log2FoldChange,
                dmg_vs_healthy_padj = padj)

de_mirna_adj_df <- read_csv(de_mirna_adjacent_file) %>%
  dplyr::select(gene_symbol, log2FoldChange,
                padj) %>%
  dplyr::rename(dmg_vs_adj_log2FC = log2FoldChange,
                dmg_vs_adj_padj = padj)

dmg_cluster6_go_df <- dmg_cluster6_go_df %>%
  left_join(de_mirna_healthy_df, 
            by = c("Genes" = "gene_symbol")) %>%
  left_join(de_mirna_adj_df,
            by = c("Genes" = "gene_symbol"))

# filter for sig interactions (neg mirna-target correlations or target sig downregulated)
sig_interactions_df <- dmg_cluster6_go_df %>%
  dplyr::filter((mirna_target_pearson_r < 0 & mirna_target_pearson_p < 0.05) |
                  (dmg_vs_healthy_log2FC < 0 & dmg_vs_healthy_padj < 0.05) |
                  (dmg_vs_adj_log2FC < 0 & dmg_vs_adj_padj < 0.05))

# save full and filtered files
write_tsv(dmg_cluster6_go_df,
          file.path(results_dir,
                    "DIPG_or_DMG-cluster6-mirna-neuronal-target-interactions.tsv"))

write_tsv(sig_interactions_df,
          file.path(results_dir,
                    "DIPG_or_DMG-cluster6-mirna-neuronal-target-sig-interactions.tsv"))

# print session info
sessionInfo()
