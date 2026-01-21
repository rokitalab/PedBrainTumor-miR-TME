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

### Tumor types to run
tumor_types <- c("DIPG or DMG", "MB")

for (tumor in tumor_types) {
  
  message("Running miRNA neuronal heatmap for: ", tumor)
  
  tumor_tag <- str_replace_all(tumor, " ", "_")
  
  ### Load correlation results
  neuronal_rho <- read_tsv(
    file.path(
      results_dir,
      paste0(
        str_replace_all(tumor, " ", "_"),
        "_mirna_neuronal_gobp_spearman.tsv"
      )
    ),
    show_col_types = FALSE
  ) %>%
    filter(Pathway %in% neuronal_gobp) %>%
    pivot_wider(
      id_cols = miRNA,
      names_from = Pathway,
      values_from = rho
    )
  
  ### Subset histology
  group_hist <- histology_df %>%
    filter(
      experimental_strategy == "miRNA-Seq",
      histology == tumor
    )
  
  group_ids <- intersect(colnames(mirna_tpm), group_hist$external_sample_id)
  
  miRNA_expr <- mirna_tpm[, group_ids, drop = FALSE]
  
  ### Normalize TPMs
  miRNA_zscores <- t(apply(miRNA_expr, 1, function(x) (x - mean(x)) / sd(x)))
  miRNA_zscores <- miRNA_zscores[rowSums(is.na(miRNA_zscores)) == 0, ]
  
  ### Load DE miRNAs (tumor-specific)
  sig_DE_miRNA_list <- read.csv(
    file.path(
      root_dir,
      "analyses",
      "human-mirna-expression",
      "results",
      paste0(tumor, "_sig_DE_miRNA_list.csv")
    ),
    stringsAsFactors = FALSE
  )
  
  de_mirnas <- sig_DE_miRNA_list$miRNA
  
  miRNA_expr_sub <- miRNA_zscores[
    rownames(miRNA_zscores) %in% de_mirnas,
    ,
    drop = FALSE
  ]
  
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
  
  ### Sample annotation
  sample_anno <- data.frame(external_sample_id = colnames(miRNA_expr_sub)) %>%
    left_join(group_hist %>% select(sample_type, external_sample_id),
              by = "external_sample_id")
  
  ### Annotations (unchanged)
  ra_left <- rowAnnotation(
    Annotated = mirna_anno$Annotated,
    Direction = mirna_anno$direction,
    col = list(
      Annotated = c("Yes" = "#1b9e77", "No" = "#d95f02"),
      Direction = c(up = "#E41A1C", down = "#377EB8")
    )
  )
  
  ### Right annotation object for correlation annotations
  ra_right <- rowAnnotation(
    "neurotransmitter\ntransport r" = anno_barplot(
      mirna_anno$GOBP_NEUROTRANSMITTER_TRANSPORT,
      # numeric vector, one value per row
      gp = gpar(fill = "steelblue"),
      # bar fill color
      width = unit(2, "cm")
    ),
    "postsynaptic membrane\norganization r" = anno_barplot(
      mirna_anno$GOBP_POSTSYNAPTIC_MEMBRANE_ORGANIZATION,
      # numeric vector, one value per row
      gp = gpar(fill = "steelblue4"),
      # bar fill color
      width = unit(2, "cm")
    ),
    "neuron cell \ncell adhesion r" = anno_barplot(
      mirna_anno$GOBP_NEURON_CELL_CELL_ADHESION,
      # numeric vector, one value per row
      gp = gpar(fill = "steelblue4"),
      # bar fill color
      width = unit(2, "cm")
    ),
    "synaptic vesicle\ntransport r" = anno_barplot(
      mirna_anno$GOBP_SYNAPTIC_VESICLE_TRANSPORT,
      # numeric vector, one value per row
      gp = gpar(fill = "steelblue4"),
      # bar fill color
      width = unit(2, "cm")
    ),
    " synaptic transmission\nglutamatergic r" = anno_barplot(
      mirna_anno$GOBP_SYNAPTIC_TRANSMISSION_GLUTAMATERGIC,
      # numeric vector, one value per row
      gp = gpar(fill = "steelblue4"),
      # bar fill color
      width = unit(2, "cm")
    ),
    " regulation of synaptic\ntransmission gabaergic r" = anno_barplot(
      mirna_anno$GOBP_REGULATION_OF_SYNAPTIC_TRANSMISSION_GABAERGIC,
      # numeric vector, one value per row
      gp = gpar(fill = "steelblue4"),
      # bar fill color
      width = unit(2, "cm")
    ),
    "  synaptic transmission\ncholinergic r" = anno_barplot(
      mirna_anno$GOBP_SYNAPTIC_TRANSMISSION_CHOLINERGIC,
      # numeric vector, one value per row
      gp = gpar(fill = "steelblue4"),
      # bar fill color
      width = unit(2, "cm")
    ),
    annotation_name_gp = gpar(fontsize = 10),
    annotation_name_rot = 45
  )
  
  ha_all <- HeatmapAnnotation(
    Condition = sample_anno$sample_type,
    col = list(
      Condition = c(
        "Tumor" = "black",
        "Normal-Tumor Adjacent" = "red3",
        "Normal" = "green4"
      )
    )
  )
  
  col_fun <- colorRamp2(c(-4, 0, 4), c("navyblue", "white", "orangered"))
  
  set.seed(123)
  
  # define cluster number
  clusters <- ifelse(tumor == "DIPG or DMG",
                     8, 6)
  
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
    file.path(plot_dir, paste0(tumor_tag, "_mirna_neuronal_gobp_heatmap.pdf")),
    width = 15,
    height = 8
  )
  draw(mirna_ht)
  dev.off()
  
  ### Cluster membership
  mirna_ht <- draw(mirna_ht)
  row_idx_list <- row_order(mirna_ht)
  row_idx <- unlist(row_idx_list)
  row_cluster <- rep(seq_along(row_idx_list), lengths(row_idx_list))
  
  mirna_clusters <- tibble(
    miRNA = rownames(miRNA_expr_sub)[row_idx],
    row_cluster = row_cluster
  ) %>%
    left_join(mirna_anno, by = "miRNA")
  
  write_tsv(
    mirna_clusters,
    file.path(
      results_dir,
      paste0(tumor_tag, "_mirna_neuronal_cluster_membership.tsv")
    )
  )
}

### Session Info
sessionInfo()


