# Plot human DMG cluster 6 mirna orthologs in mouse samples
#
# Ryan Corbett
#
# Oct 2025

# This script plots human DMG cluster 6 miRNA orthologs that are also DE in B7H3 treated mice

### Setup
suppressPackageStartupMessages({
  library(rprojroot)
  library(tidyverse)
  library(ggplot2)
  library(gtools)
  library(ggsci)
})

### Paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "human-mirna-expr-in-mouse")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

source(file.path(root_dir, "figures", "theme.R"))


## set file paths
dmg_cluster_file <- file.path(root_dir, "analyses",
                              "mirna-clustering",
                              "results",
                              "DIPG or DMG-de-mirna-cluster-membership-immune-scores.tsv")

mmu_meta_file <- file.path(root_dir,
                           "analyses", "murine-mirna-differential-expression",
                           "results", "sample-metadata.tsv")

mmu_mir_tpm_file <- file.path(root_dir, "analyses",
                              "murine-mirna-clustering",
                              "results",
                              "murine-mirna-tpm.rds")

mmu_de_file <- file.path(root_dir, "analyses",
                           "murine-mirna-differential-expression",
                           "results",
                           "mirna-differential-expression-deseq2-b7h3-stop-vs-untreated-by-timepoint.tsv")

# Load hsa DMG miRNA cluster df and filter for cluster 6 miRNAs
dmg_hsa_mir_df <- read_tsv(dmg_cluster_file)

dmg_clust6_mirs <- dmg_hsa_mir_df %>%
  dplyr::filter(row_cluster == 6,
                grepl("hsa", miRNA)) %>%
  # remove hsa prefix from miRNA ID
  dplyr::mutate(miRNA = str_remove(miRNA, "hsa-")) %>%
  pull(miRNA)

# Load merged mouse differential expression results
mmu_de_df <- read_tsv(mmu_de_file) %>%
  # rm species prefix from miRNA IDs
  dplyr::mutate(mirna_id = str_remove(mirna_id, "mmu-"))

# Get DE cluster6 orthologs 
de_clust6_orthos <- mmu_de_df %>%
  dplyr::filter(mirna_id %in% dmg_clust6_mirs) %>%
  dplyr::filter(expr_pattern != "Not DE") %>%
  pull(mirna_id)

# Load mouse sample metadata
mmu_meta_df <- read_tsv(mmu_meta_file)

# Load TPM matrix and filter for DE cluster6 orthologs
mmu_mirna_tpm_df <- readRDS(mmu_mir_tpm_file) %>%
  rownames_to_column("mirna_id") %>%
  dplyr::mutate(mirna_id = str_remove(mirna_id, "mmu-")) %>%
  column_to_rownames("mirna_id") %>%
  t() %>%
  as.data.frame() %>%
  dplyr::select(any_of(de_clust6_orthos)) %>%
  rownames_to_column("sample_id") %>%
  pivot_longer(-sample_id,
               names_to = "mirna_id",
               values_to = "tpm") %>%
  dplyr::mutate(log2_tpm = log2(tpm + 1)) %>%
  left_join(mmu_meta_df, 
            by = c("sample_id" = "id")) %>%
  # define high and low CAR T responders
  dplyr::mutate(car_t_response = case_when(
    treatment == "B7H3" & sample_id %in% c("Day14-Tube3", "Day21-Tube9") ~ "Low",
    treatment == "B7H3" ~ "High",
    TRUE ~ "N/A"
  ))

# generate TPM plots
expr_plot <- mmu_mirna_tpm_df %>%
  ggplot(aes(x = time,
             y = log2_tpm,
             fill = treatment,
             color = treatment,
             group = treatment)) +
  geom_point(aes(shape = car_t_response),
             size = 3,
             alpha = 0.6,
             position = position_dodge(width = 0.3)) +
  # error bars per treatment
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               width = 0.2,
               position = position_dodge(width = 0.3)) +
  scale_shape_manual(
    values = c("High" = 24,
               "Low"  = 25,
               "N/A"  = 21)) +
  scale_color_npg() +
  scale_fill_npg() +
  facet_wrap(~mirna_id, ncol = 3, scales = "free_y") +
  labs(x = NULL,
       y = "log2-TPM",
       color = "Treatment",
       shape = "CAR-T response") +
  guides(fill = "none") +
  theme_Publication()

# save plot
ggsave(file.path(plot_dir,
                 "hsa-dmg-clust6-mouse-ortholog-de-tpm.pdf"),
       width = 8, height = 5)

# print session info
sessionInfo()
