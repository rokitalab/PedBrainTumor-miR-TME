# Assess MB cluster 1 miRNA - neuronal target interactions
# Ryan Corbett | Oct 2025
# Updated by Bicna Song | Jan 2026 

# This script assesses correlations between MB cluster 1 miRNA and immune pathway target expression, and assesses coordinated DE of targets to identify significant associations

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

# figure theme
source(file.path(root_dir, "figures", "theme.R"))

# set file paths
mb_cluster1_go_file <- file.path(results_dir,
                                  "MB-cluster1-mirna-target-go-enr-neuronal-terms-reduced.tsv")

# mirTarbase targets
mirtarbase_file <- file.path(root_dir, "analyses",
                             "human-mirna-go-enrichment",
                             "results",
                             "human-mirna-target-predictions-miRTarBase.tsv.gz")

hist_file <- file.path(root_dir, "analyses",
                       "histology-preprocessing",
                       "results",
                       "histologies.tsv")

# TPM files
mirna_tpm_file <- file.path(root_dir,
                            "analyses",
                            "immune-deconvolution",
                            "results",
                            "mirna-tpm.rds")

target_tpm_file <- file.path(data_dir,
                             "gene-expression-rsem-tpm-collapsed.all.rds")

# DE bulk RNA-seq files
de_mirna_healthy_file <- file.path(root_dir, 
                                   "analyses",
                                   "human-rna-expression",
                                   "results",
                                   "DESeq2_MB_vs_healthyNormal.csv")

de_mirna_adjacent_file <- file.path(root_dir, 
                                    "analyses",
                                    "human-rna-expression",
                                    "results",
                                    "DESeq2_MB_paired_full.csv")

# miranda-predicted targets for novel miRNAs
miranda_file <- file.path(root_dir, "analyses",
                          "human-novel-mirna-target-prediction",
                          "results",
                          "human_miranda_output_parsed_anno.csv")

# Wrangle data

mb_cluster1_go_df <- read_tsv(mb_cluster1_go_file) %>%
  dplyr::filter(!is.na(Term)) %>%
  dplyr::select(miRNA, Term, Genes) %>%
  # expand to one target gene per row
  separate_rows(Genes, sep = ",\\s*") %>%
  # collapse enriched terms per target
  group_by(miRNA, Genes) %>%
  summarise(Term = str_c(Term, collapse = ";"))

mirtarbase_df <- read_tsv(mirtarbase_file) %>%
  dplyr::filter(`Species (miRNA)` == "hsa") %>% 
  distinct(miRNA, `Target Gene`, .keep_all = TRUE) %>%
  dplyr::select(miRNA, `Target Gene`,
                Experiments, `Support Type`)

miranda_df <- read_csv(miranda_file) %>% 
  distinct(Seq1, GeneSymbol, .keep_all = TRUE) %>%
  dplyr::rename(miRNA = Seq1,
                `Target Gene` = GeneSymbol)

# merge annotated and novel miRNA target predictions
target_df <- mirtarbase_df %>%
  bind_rows(miranda_df) %>%
  dplyr::mutate(`Support Type` = case_when(
    is.na(`Support Type`) ~ "miRanda",
    TRUE ~ `Support Type`
  ))

# Join target with GO data
mb_cluster1_go_df <- mb_cluster1_go_df %>%
  left_join(target_df,
            by = c("miRNA",
                   "Genes" = "Target Gene")) %>%
  # only retain genes that are predicted as target
  dplyr::filter(!is.na(`Support Type`))


# Load hist and extract MB miRNA and RNA IDs

hist <- read_tsv(hist_file)

mb_mirna_ids <- hist %>% 
  dplyr::filter(experimental_strategy == "miRNA-Seq",
                sample_type == "Tumor",
                histology == "MB") %>%
  arrange(external_sample_id) %>%
  pull(external_sample_id)

mb_rna_ids <- hist %>% 
  dplyr::filter(experimental_strategy == "RNA-Seq",
                sample_type == "Tumor",
                histology == "MB") %>%
  arrange(external_sample_id) %>%
  pull(Bioassay_ID)

# load TPM files 

mirna_tpm <- readRDS(mirna_tpm_file) %>%
  dplyr::select(any_of(c(mb_mirna_ids)))

target_tpm <- readRDS(target_tpm_file) %>%
  dplyr::select(any_of(c(mb_rna_ids)))

# add empty pearson correlation coefficient and p-value columns to cluster 6 df

mb_cluster1_go_df <- mb_cluster1_go_df %>%
  dplyr::mutate(mirna_target_pearson_r = NA_integer_,
                mirna_target_pearson_p = NA_integer_)

# loop through miRNA-target pairs
for (i in 1:nrow(mb_cluster1_go_df)){
  
  miRNA <- mb_cluster1_go_df$miRNA[i]
  target <- mb_cluster1_go_df$Genes[i]
  
  # check if target is in TPM mat
  if (!target %in% rownames(target_tpm)) next
  
  # get TPMs
  miRNA_tpms <- unlist(mirna_tpm[miRNA, mb_mirna_ids])
  target_tpms <- unlist(target_tpm[target, mb_rna_ids])
  
  # calculate correlation stats
  mb_cluster1_go_df$mirna_target_pearson_r[i] <- cor.test(log2(miRNA_tpms + 0.01),
                                                           log2(target_tpms + 0.01),
                                                           method = "pearson")$estimate
  
  mb_cluster1_go_df$mirna_target_pearson_p[i] <- cor.test(log2(miRNA_tpms + 0.01),
                                                           log2(target_tpms + 0.01),
                                                           method = "pearson")$p.value
  
}

# Load DE results

de_mirna_healthy_df <- read_csv(de_mirna_healthy_file) %>%
  dplyr::select(gene_symbol, log2FoldChange,
                padj) %>%
  dplyr::rename(mb_vs_healthy_log2FC = log2FoldChange,
                mb_vs_healthy_padj = padj)

de_mirna_adj_df <- read_csv(de_mirna_adjacent_file) %>%
  dplyr::select(gene_symbol, log2FoldChange,
                padj) %>%
  dplyr::rename(mb_vs_adj_log2FC = log2FoldChange,
                mb_vs_adj_padj = padj)

# join DE results to GO target df 
mb_cluster1_go_df <- mb_cluster1_go_df %>%
  left_join(de_mirna_healthy_df, 
            by = c("Genes" = "gene_symbol")) %>%
  left_join(de_mirna_adj_df,
            by = c("Genes" = "gene_symbol"))

# filter for sig interactions (neg mirna-target correlations or target sig downregulated)
sig_interactions_df <- mb_cluster1_go_df %>%
  dplyr::filter((mirna_target_pearson_r < 0 & mirna_target_pearson_p < 0.05) |
                  (mb_vs_healthy_log2FC < 0 & mb_vs_healthy_padj < 0.05) |
                  (mb_vs_adj_log2FC < 0 & mb_vs_adj_padj < 0.05))

# save full and filtered files
write_tsv(mb_cluster1_go_df,
          file.path(results_dir,
                    "MB-cluster1-mirna-immune-target-interactions.tsv"))

write_tsv(sig_interactions_df,
          file.path(results_dir,
                    "MB-cluster1-mirna-immune-target-sig-interactions.tsv"))

# print session info
sessionInfo()
