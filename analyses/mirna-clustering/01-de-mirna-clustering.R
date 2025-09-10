# Cluster DIPG/DMG and MB DE miRNAs by expression
#
# Ryan Corbett
#
# Sep 2025

### Overview

# Cluster DE miRNAs by normalized expression patterns 
# Append immune cell fraction and score correlation coefficients to determine which clusters are associated with immune cell fractions
  
  
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
analysis_dir <- file.path(root_dir, "analyses", "miRNA-clustering")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")
input_dir <- file.path(analysis_dir, "input")

miRNA_tpm <- file.path(root_dir, 
                       "analyses",
                       "immune-deconvolution",
                       "results",
                       "mirna-tpm.rds")
histology_file <- file.path(root_dir, "analyses", 
                            "histology-preprocessing", 
                            "results", "histologies.tsv")

### Load Data

# miRNA TPM
mirna_tpm <- readRDS(miRNA_tpm)
mirna_tpm <- log2(as.matrix(mirna_tpm) + 1)

# Histology
histology_df <- read_delim(histology_file, show_col_types = FALSE)


groups <- c("DIPG or DMG", "MB")
normal_ids <- c("B1","B2","B3","B4","B5")

for (group in groups){
  
  xcell_scores <- read_tsv(file.path(root_dir, 
                                     "analyses",
                                     "immune-deconvolution",
                                     "results",
                                     glue::glue("{group}-de-mirna-xcell-fraction-correlations.tsv"))) %>%
    dplyr::filter(cell_type %in% c("immune score", "microenvironment score",
                                   "Macrophage M1", "Macrophage M2")) %>% 
    dplyr::mutate(cell_type = glue::glue("{cell_type} (x)")) %>%
    pivot_wider(id_cols = miRNA,
                names_from = cell_type,
                values_from = pearson_r)
  
  quantiseq_scores <- read_tsv(file.path(root_dir, 
                                         "analyses",
                                         "immune-deconvolution",
                                         "results",
                                         glue::glue("{group}-de-mirna-quantiseq-fraction-correlations.tsv"))) %>%
    dplyr::filter(cell_type %in% c("NK cell", "T cell regulatory (Tregs)",
                                   "B cell", "T cell CD4+ (non-regulatory)")) %>% 
    dplyr::mutate(cell_type = glue::glue("{cell_type} (q)")) %>%
    pivot_wider(id_cols = miRNA,
                names_from = cell_type,
                values_from = pearson_r)
  
  
  ### Select DIPG/DMG tumor samples and subset matrices
  
  group_hist <- histology_df %>%
    filter(experimental_strategy == "RNA-Seq",
           histology == group)
  
  group_ids <- group_hist %>%
    pull(external_sample_id)
  
  # Subset GSVA and miRNA matrices
  miRNA_expr <- mirna_tpm  [, group_ids, drop=FALSE]
  
  miRNA_zscores <- t(apply(miRNA_expr, 1, function(x) (x - mean(x)) / sd(x)))
  miRNA_zscores <- miRNA_zscores[rowSums(is.nan(miRNA_zscores)) == 0,]
  
  ### Correlate miRNA Expression with Pathway Scores
  
  # Select DE miRNAs of interest (170 genes)
  sig_DE_miRNA_list <- read.csv(file.path(root_dir, "analyses", 
                                          "human-mirna-expression",
                                          "results", 
                                         glue::glue("{group}_sig_DE_miRNA_list.csv")))
  
  de_mirnas <- sig_DE_miRNA_list$miRNA
  miRNA_expr_sub <- miRNA_zscores[rownames(miRNA_zscores) %in% de_mirnas, , drop=FALSE]
  
  ### Heatmap: all pairs
  
  # Annotation mapping for "all pairs" view 
  mirna_anno <- data.frame(miRNA = rownames(miRNA_expr_sub)) %>%
    left_join(sig_DE_miRNA_list, by = "miRNA") %>%
    left_join(xcell_scores) %>%
    left_join(quantiseq_scores) %>%
    dplyr::mutate(Annotated = case_when(
      grepl("hsa", miRNA) ~ "Yes",
      TRUE ~ "No"
    )) %>%
    dplyr::mutate(Annotated = factor(Annotated,
                                    levels = c("Yes", "No"))) %>%
    distinct(miRNA, .keep_all = TRUE)
  
  mirna_anno$direction <- factor(mirna_anno$direction,
                                 levels = c("up", "down"))
  
  anno_cols <- c("Yes" = "#1b9e77",
            "No" = "#d95f02")
  dir_cols  <- c(up = "#E41A1C", down = "#377EB8")
  
  sample_anno <- data.frame(external_sample_id = colnames(miRNA_expr_sub)) %>%
    left_join(group_hist %>% dplyr::select(sample_type,
                                          external_sample_id))
  
  # Color mapping for correlation
  ra_left <- rowAnnotation(
    Annotated = mirna_anno$Annotated,
    Direction  = mirna_anno$direction,
    col = list(Annotated = anno_cols, Direction = dir_cols),
    annotation_name_gp = gpar(fontsize = 10),
    annotation_legend_param = list(
      Direction  = list(title = "Direction",  at = levels(mirna_anno$direction)),
      Annotated  = list(title = "Annotated",  at = levels(mirna_anno$Annotated))
    )
  )
  
  ra_right <- rowAnnotation(
    "Immune\nscore (x)" = anno_barplot(
      mirna_anno$`immune score (x)`,     # numeric vector, one value per row
      gp = gpar(fill = "steelblue"), # bar fill color
      width = unit(2, "cm")),
    "TME\nscore (x)" = anno_barplot(
      mirna_anno$`microenvironment score (x)`,     # numeric vector, one value per row
      gp = gpar(fill = "steelblue4"), # bar fill color
      width = unit(2, "cm")),
    "Macrophage\nM1 (x)" = anno_barplot(
      mirna_anno$`Macrophage M1 (x)`,     # numeric vector, one value per row
      gp = gpar(fill = "steelblue4"), # bar fill color
      width = unit(2, "cm")),
    "Macrophage\nM2 (x)" = anno_barplot(
      mirna_anno$`Macrophage M2 (x)`,     # numeric vector, one value per row
      gp = gpar(fill = "steelblue4"), # bar fill color
      width = unit(2, "cm")),
    "B cell\n(q)" = anno_barplot(
      mirna_anno$`B cell (q)`,     # numeric vector, one value per row
      gp = gpar(fill = "steelblue4"), # bar fill color
      width = unit(2, "cm")),
    "NK cell\n(q)" = anno_barplot(
      mirna_anno$`NK cell (q)`,     # numeric vector, one value per row
      gp = gpar(fill = "steelblue4"), # bar fill color
      width = unit(2, "cm")),
    "Treg\n(q)" = anno_barplot(
      mirna_anno$`T cell regulatory (Tregs) (q)`,     # numeric vector, one value per row
      gp = gpar(fill = "steelblue4"), # bar fill color
      width = unit(2, "cm")),
    "CD4+ Tcell\n(q)" = anno_barplot(
      mirna_anno$`T cell CD4+ (non-regulatory) (q)`,     # numeric vector, one value per row
      gp = gpar(fill = "steelblue4"), # bar fill color
      width = unit(2, "cm")),
    annotation_name_gp = gpar(fontsize = 10)
  )
  
  ha_all <- HeatmapAnnotation(
    Condition = sample_anno$sample_type,
    col = list(Condition = c("Tumor" = "black",
                             "Normal-Tumor Adjacent" = "red3",
                             "Normal" = "green4")),
    annotation_name_gp = gpar(fontsize = 10)
  )
  
  col_fun <- colorRamp2(c(-4, 0, 4), c("navyblue", "white", "orangered"))
  
  ### Determining Optimal Clusters
  
  set.seed(123)
  
  clusters <- ifelse(group == "DIPG or DMG",
                     9, 6)
  
  ht2 <- Heatmap(
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
    heatmap_legend_param = list(legend_gp = gpar(fontsize = 10),
                                labels_gp = gpar(fontsize = 10))
  )
  
  ht <- ifelse(group == "DIPG or DMG",
                    8, 5)
  
  pdf(file.path(plot_dir,
                glue::glue("{group}-de-mirna-heatmap.pdf")),
      width = 15, height = ht)
  
  ht2 <- draw(ht2)
  
  dev.off()
  
  row_idx <- unlist(row_order(ht2))
  col_idx <- unlist(column_order(ht2))
  
  row_idx_list <- row_order(ht2)
  col_idx_list <- column_order(ht2)
  
  row_cluster_by_order <- rep(seq_along(row_idx_list), lengths(row_idx_list))
  col_cluster_by_order <- rep(seq_along(col_idx_list), lengths(col_idx_list))
  
  ## Pathway (row) clusters
  mirna_clusters <- tibble::tibble(
    miRNA = rownames(miRNA_expr_sub)[row_idx],
    row_cluster = row_cluster_by_order
  )
  
  mirna_clusters <- mirna_clusters %>%
    left_join(mirna_anno) 
  
  write_tsv(mirna_clusters,
            file.path(results_dir,
                      glue::glue("{group}-de-mirna-cluster-membership-immune-scores.tsv")))
  
}

### Session Info

sessionInfo()
