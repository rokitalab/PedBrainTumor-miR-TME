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

# Set file paths

miRNA_tpm <- file.path(root_dir, 
                       "analyses",
                       "immune-deconvolution",
                       "results",
                       "mirna-tpm.rds")
histology_file <- file.path(root_dir, "analyses", 
                            "histology-preprocessing", 
                            "results", "histologies.tsv")

# Wrangle data

# miRNA TPM
mirna_tpm <- readRDS(miRNA_tpm)
mirna_tpm <- log2(as.matrix(mirna_tpm) + 1)

# Histology
histology_df <- read_delim(histology_file, show_col_types = FALSE) %>%
  dplyr::filter(external_sample_id != "1-1855-CC1")

# define histology groups to cluster & plot
groups <- c("DIPG or DMG", "MB")

# loop through groups
for (group in groups){
  
  # load histology-specific miRNA-immune cell fraction correlation results
  
  # XCell results
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
  
  # quanti-seq results
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
  
  
  # T cell marker gene corelation results
  Tcell_gene_scores <- read_tsv(file.path(root_dir, 
                                         "analyses",
                                         "immune-deconvolution",
                                         "results",
                                         glue::glue("{group}-de-mirna-tcell-marker-gene-correlations.tsv"))) %>%
    dplyr::filter(Gene_symbol %in% c("CD3D", "CD3G", "CD4", "CD8A", "IFNG",
                                   "GZMA", "GZMB", "PRF1")) %>% 
    pivot_wider(id_cols = miRNA,
                names_from = Gene_symbol,
                values_from = pearson_r)
  
  
  # Subset histologies file and pull miRNA sample IDS
  group_hist <- histology_df %>%
    filter(experimental_strategy == "miRNA-Seq",
           histology == group)
  
  group_ids <- group_hist %>%
    pull(external_sample_id)
  
  # Subset miRNA tpm matrix
  miRNA_expr <- mirna_tpm  [, group_ids, drop=FALSE]
  
  # normalize miRNA TPMs
  miRNA_zscores <- t(apply(miRNA_expr, 1, function(x) (x - mean(x)) / sd(x)))
  miRNA_zscores <- miRNA_zscores[rowSums(is.nan(miRNA_zscores)) == 0,]
  
  # Load DE miRNAs
  sig_DE_miRNA_list <- read.csv(file.path(root_dir, "analyses", 
                                          "human-mirna-expression",
                                          "results", 
                                         glue::glue("{group}_sig_DE_miRNA_list.csv")))
  
  de_mirnas <- sig_DE_miRNA_list$miRNA
  
  # subset miRNA expr matrix
  miRNA_expr_sub <- miRNA_zscores[rownames(miRNA_zscores) %in% de_mirnas, , drop=FALSE]
  
  # Create miRNA annotation df
  mirna_anno <- data.frame(miRNA = rownames(miRNA_expr_sub)) %>%
    left_join(sig_DE_miRNA_list, by = "miRNA") %>%
    # add immune cell fraction correlation coefficients
    left_join(xcell_scores) %>%
    left_join(quantiseq_scores) %>%
    left_join(Tcell_gene_scores) %>%
    # add column indicating if miRNA is annotated or novel
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
  
  # Create sample annotation df
  sample_anno <- data.frame(external_sample_id = colnames(miRNA_expr_sub)) %>%
    left_join(group_hist %>% dplyr::select(sample_type,
                                          external_sample_id))
  
  # Left annotation object for discrete miRNA annotations 
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
  
  # Right annotation object for correlation annotations 
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
    "CD8A\nTPM r" = anno_barplot(
      mirna_anno$CD8A,     # numeric vector, one value per row
      gp = gpar(fill = "steelblue4"), # bar fill color
      width = unit(2, "cm")),
    "CD4\nTPM r" = anno_barplot(
      mirna_anno$CD4,     # numeric vector, one value per row
      gp = gpar(fill = "steelblue4"), # bar fill color
      width = unit(2, "cm")),
    annotation_name_gp = gpar(fontsize = 10)
  )
  
  # sample annotation 
  ha_all <- HeatmapAnnotation(
    Condition = sample_anno$sample_type,
    col = list(Condition = c("Tumor" = "black",
                             "Normal-Tumor Adjacent" = "red3",
                             "Normal" = "green4")),
    annotation_name_gp = gpar(fontsize = 10)
  )
  
  col_fun <- colorRamp2(c(-4, 0, 4), c("navyblue", "white", "orangered"))
  
  set.seed(123)
  
  # define cluster number
  clusters <- ifelse(group == "DIPG or DMG",
                     8, 6)
  
  # Generate heatmap
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
    heatmap_legend_param = list(legend_gp = gpar(fontsize = 10),
                                labels_gp = gpar(fontsize = 10))
  )
  
  # Save heatmap
  ht <- ifelse(group == "DIPG or DMG",
                    8, 5)
  
  pdf(file.path(plot_dir,
                glue::glue("{group}-de-mirna-heatmap.pdf")),
      width = 15, height = ht)
  
  mirna_ht <- draw(mirna_ht)
  
  dev.off()
  
  # create miRNA df that includes cluster assignment
  
  row_idx <- unlist(row_order(mirna_ht))
  row_idx_list <- row_order(mirna_ht)
  row_cluster_by_order <- rep(seq_along(row_idx_list), lengths(row_idx_list))

  mirna_clusters <- tibble::tibble(
    miRNA = rownames(miRNA_expr_sub)[row_idx],
    row_cluster = row_cluster_by_order
  )
  
  # add other annotation columns
  mirna_clusters <- mirna_clusters %>%
    left_join(mirna_anno) 
  
  # write to output
  write_tsv(mirna_clusters,
            file.path(results_dir,
                      glue::glue("{group}-de-mirna-cluster-membership-immune-scores.tsv")))
  
}

### Session Info

sessionInfo()
