# Regenerate miRNA clustering heatmap with neuronal gobp correlation coefficient as annotation
# The script has been adpated from /mirna-clustering/01-de-mirna-clustering.R
# Author: Bicna Song | Jan 2026

### Setup
library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(rprojroot)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(cluster)
library(purrr)
library(factoextra)

### Paths and Files
root_dir <- find_root(has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "human-neuronal-expression")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

miRNA_tpm <- file.path(root_dir,
                       "analyses",
                       "immune-deconvolution",
                       "results",
                       "mirna-tpm.rds")
histology_file <- file.path(root_dir,
                            "analyses",
                            "histology-preprocessing",
                            "results",
                            "histologies.tsv")

### Wrangle data

### miRNA TPM
mirna_tpm <- readRDS(miRNA_tpm)
mirna_tpm <- log2(as.matrix(mirna_tpm) + 1)

### Histology
histology_df <- read_delim(histology_file, show_col_types = FALSE) %>%
  dplyr::filter(external_sample_id != "1-1855-CC1")

### Define identified neuronal-related GO terms
neuronal_gobp <- c(
  "GOBP_NEUROTRANSMITTER_TRANSPORT",
  "GOBP_POSTSYNAPTIC_MEMBRANE_ORGANIZATION",
  "GOBP_NEURON_CELL_CELL_ADHESION",
  "GOBP_SYNAPTIC_VESICLE_TRANSPORT",
  "GOBP_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC",
  "GOBP_REGULATION_OF_SYNAPTIC_TRANSMISSION_GABAERGIC",
  "GOBP_SYNAPTIC_TRANSMISSION_CHOLINERGIC"
)

### Load the correlation analysis result between DE miRNAs and GSVA scores of neuronal-related GO terms
neuronal_rho <- read_tsv(file.path(results_dir, "mirna_neuronal_gobp_spearman.tsv")) %>%
  dplyr::filter(Pathway %in% neuronal_gobp) %>%
  pivot_wider(id_cols = miRNA,
              names_from = Pathway,
              values_from = rho)

### Subset histologies file and pull miRNA sample IDS
group_hist <- histology_df %>%
  filter(experimental_strategy == "miRNA-Seq",
         histology == "DIPG or DMG")

group_ids <- group_hist %>%
  pull(external_sample_id)

### Subset miRNA tpm matrix
miRNA_expr <- mirna_tpm  [, group_ids, drop = FALSE]

### Normalize miRNA TPMs
miRNA_zscores <- t(apply(miRNA_expr, 1, function(x) (x - mean(x)) / sd(x)))
miRNA_zscores <- miRNA_zscores[rowSums(is.nan(miRNA_zscores)) == 0, ]

### Load DE miRNAs
sig_DE_miRNA_list <- read.csv(
  file.path(
    root_dir,
    "analyses",
    "human-mirna-expression",
    "results",
    "DIPG or DMG_sig_DE_miRNA_list.csv"
  )
)

de_mirnas <- sig_DE_miRNA_list$miRNA

### Subset miRNA expr matrix
miRNA_expr_sub <- miRNA_zscores[rownames(miRNA_zscores) %in% de_mirnas, , drop = FALSE]

### Create miRNA annotation df
mirna_anno <- data.frame(miRNA = rownames(miRNA_expr_sub)) %>%
  left_join(sig_DE_miRNA_list, by = "miRNA") %>%
  # add neuronal gobp correlation coefficients
  left_join(neuronal_rho) %>%
  # add column indicating if miRNA is annotated or novel
  dplyr::mutate(Annotated = case_when(grepl("hsa", miRNA) ~ "Yes", TRUE ~ "No")) %>%
  dplyr::mutate(Annotated = factor(Annotated, levels = c("Yes", "No"))) %>%
  distinct(miRNA, .keep_all = TRUE)

mirna_anno$direction <- factor(mirna_anno$direction, levels = c("up", "down"))

anno_cols <- c("Yes" = "#1b9e77", "No" = "#d95f02")
dir_cols  <- c(up = "#E41A1C", down = "#377EB8")

### Create sample annotation df
sample_anno <- data.frame(external_sample_id = colnames(miRNA_expr_sub)) %>%
  left_join(group_hist %>% dplyr::select(sample_type, external_sample_id))

### Left annotation object for discrete miRNA annotations
ra_left <- rowAnnotation(
  Annotated = mirna_anno$Annotated,
  Direction  = mirna_anno$direction,
  col = list(Annotated = anno_cols, Direction = dir_cols),
  annotation_name_gp = gpar(fontsize = 10),
  annotation_legend_param = list(
    Direction  = list(title = "Direction", at = levels(mirna_anno$direction)),
    Annotated  = list(title = "Annotated", at = levels(mirna_anno$Annotated))
  )
)

### Right annotation object for correlation annotations
ra_right <- rowAnnotation(
  "Neurotransmitter\nreceptors r" = anno_barplot(
    mirna_anno$GOBP_NEUROTRANSMITTER_TRANSPORT,
    # numeric vector, one value per row
    gp = gpar(fill = "steelblue"),
    # bar fill color
    width = unit(2, "cm")
  ),
  "Postsynaptic\nscaffolding r" = anno_barplot(
    mirna_anno$GOBP_POSTSYNAPTIC_MEMBRANE_ORGANIZATION,
    # numeric vector, one value per row
    gp = gpar(fill = "steelblue4"),
    # bar fill color
    width = unit(2, "cm")
  ),
  "Cell\nadhesion r" = anno_barplot(
    mirna_anno$GOBP_NEURON_CELL_CELL_ADHESION,
    # numeric vector, one value per row
    gp = gpar(fill = "steelblue4"),
    # bar fill color
    width = unit(2, "cm")
  ),
  "Synaptic\nvesicle\ntrafficking r" = anno_barplot(
    mirna_anno$GOBP_SYNAPTIC_VESICLE_TRANSPORT,
    # numeric vector, one value per row
    gp = gpar(fill = "steelblue4"),
    # bar fill color
    width = unit(2, "cm")
  ),
  "Glutamatergic\nr" = anno_barplot(
    mirna_anno$GOBP_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC,
    # numeric vector, one value per row
    gp = gpar(fill = "steelblue4"),
    # bar fill color
    width = unit(2, "cm")
  ),
  "GABAergic\nr" = anno_barplot(
    mirna_anno$GOBP_REGULATION_OF_SYNAPTIC_TRANSMISSION_GABAERGIC,
    # numeric vector, one value per row
    gp = gpar(fill = "steelblue4"),
    # bar fill color
    width = unit(2, "cm")
  ),
  "Cholinergic\nr" = anno_barplot(
    mirna_anno$GOBP_SYNAPTIC_TRANSMISSION_CHOLINERGIC,
    # numeric vector, one value per row
    gp = gpar(fill = "steelblue4"),
    # bar fill color
    width = unit(2, "cm")
  ),
  annotation_name_gp = gpar(fontsize = 8)
)

### Sample annotation
ha_all <- HeatmapAnnotation(
  Condition = sample_anno$sample_type,
  col = list(
    Condition = c(
      "Tumor" = "black",
      "Normal-Tumor Adjacent" = "red3",
      "Normal" = "green4"
    )
  ),
  annotation_name_gp = gpar(fontsize = 10)
)

col_fun <- colorRamp2(c(-4, 0, 4), c("navyblue", "white", "orangered"))

set.seed(123)

### Define cluster number
clusters <- 8

### Generate heatmap
mirna_ht <- Heatmap(
  miRNA_expr_sub,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  row_split = clusters,
  column_split = 2,
  col = col_fun,
  na_col = "gray",
  name = "TPM z-score",
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = ha_all,
  left_annotation = ra_left,
  right_annotation = ra_right,
  use_raster = TRUE,
  heatmap_legend_param = list(
    legend_gp = gpar(fontsize = 10),
    labels_gp = gpar(fontsize = 10)
  )
)

### Save heatmap
pdf(
  file.path(plot_dir, "DIPG_or_DMG_mirna_neuronal_gobp_heatmap.pdf"),
  width = 15,
  height = 8
)

mirna_ht <- draw(mirna_ht)

dev.off()

### create miRNA df that includes cluster assignment
row_idx <- unlist(row_order(mirna_ht))
row_idx_list <- row_order(mirna_ht)
row_cluster_by_order <- rep(seq_along(row_idx_list), lengths(row_idx_list))

mirna_clusters <- tibble::tibble(miRNA = rownames(miRNA_expr_sub)[row_idx], row_cluster = row_cluster_by_order)

# add other annotation columns
mirna_clusters <- mirna_clusters %>%
  left_join(mirna_anno)

# write to output
write_tsv(
  mirna_clusters,
  file.path(
    results_dir,
    "DIPG_or_DMG_mirna_neuronal_cluster_memebership.tsv"
  )
)

### Session Info
sessionInfo()
