# Generate xCell2-scaled enrichment heatmaps for each reference dataset

# Author: Bicna Song

# Load libraries
suppressPackageStartupMessages({
  library(rprojroot)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
})

# Set path to module and results directories
root_dir <- find_root(has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "immune-deconvolution")
input_dir <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

# File paths
quantiseq_file <- file.path(results_dir, "quantiseq_output.rds")
xcell_file <- file.path(results_dir, "xcell_output.rds")
metadata_file <- file.path(root_dir, "analyses", 
                           "histology-preprocessing", 
                           "results", 
                           "histologies.tsv")

# wrangle data
quantiseq_results <- readRDS(quantiseq_file) %>%
  dplyr::filter(!histology %in% c("ATRT", "HGG")) %>%
  dplyr::mutate(histology = case_when(
    sample_type == "Tumor" ~ histology,
    TRUE ~ "Normal"
  ))

xcell_results <- readRDS(xcell_file) %>%
  dplyr::filter(!histology %in% c("ATRT", "HGG")) %>%
  dplyr::mutate(histology = case_when(
    sample_type == "Tumor" ~ histology,
    TRUE ~ "Normal"
  ))

# Filter to keep only selected immune-related cell types
remove <- c("immune score", "microenvironment score", 
            "stroma score", "uncharacterized cell",
            "Endothelial cell", "Cancer associated fibroblast",
            "Common lymphoid progenitor", "Common myeloid progenitor",
            "Granulocyte-monocyte progenitor", "Hematopoietic stem cell")

deconv_list <- list("Quantiseq" = quantiseq_results,
                    "XCell" = xcell_results)

for (method in names(deconv_list)){
  
  res <- deconv_list[[method]]

  # Scale scores within each reference
  scaled_results <- res %>%
    dplyr::filter(!cell_type %in% remove) %>% 
    group_by(cell_type) %>%
    mutate(Scaled_Score = if (sd(fraction, na.rm = TRUE) == 0) 0 else scale(fraction)[, 1]) %>%
    ungroup()
  
  # Save scaled results
  scaled_results_wide <- scaled_results %>%
    select(cell_type, Bioassay_ID, Scaled_Score) %>%
    pivot_wider(
      names_from = Bioassay_ID,
      values_from = Scaled_Score
    )
  
  # Save scaled results
  write.table(scaled_results_wide, file.path(results_dir, "xCell_scaled_fraction_mat.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  message("✅ Filtered and scaled enrichment scores saved.")
  
  # Generate heatmaps per reference
  col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

  # Custom metadata colors
  hist_cols <- c(
    "DIPG or DMG" = "#ff40d9",  # warm taupe
    "MB" = "#a340ff",   # muted orange
    "EPN" = "#2200ff",
    "LGG" = "#8f8fbf",
    "Normal" = "grey90"
  )
  
  type_cols <- c(
    "Tumor" = "#E64B35FF",
    "Normal" = "#4DBBD5FF",
    "Normal-Tumor Adjacent" = "#00A087FF"
  )
  
  # create matrix
  mat <- scaled_results %>%
    select(cell_type, Bioassay_ID, Scaled_Score) %>%
    pivot_wider(names_from = Bioassay_ID, values_from = Scaled_Score) %>%
    column_to_rownames("cell_type") %>%
    as.matrix()
  
  # Remove rows with zero variance
  non_constant <- apply(mat, 1, function(x) sd(x, na.rm = TRUE) > 0)
  mat <- mat[non_constant, , drop = FALSE]
  
  if (nrow(mat) < 2 || ncol(mat) < 2) {
    message(paste0("⚠️ Skipping ", ref, ": insufficient variable data for clustering."))
    next
  }
  
  # Match metadata for available samples
  meta_ref <- res %>%
    distinct(Bioassay_ID, .keep_all = TRUE) %>%
    filter(Bioassay_ID %in% colnames(mat)) %>%
    column_to_rownames("Bioassay_ID")
  
  # Ensure column order matches matrix columns
  meta_ref <- meta_ref[colnames(mat), , drop = FALSE]
  
  # Define annotations
  ha_col <- HeatmapAnnotation(
    `Sample Type` = meta_ref$sample_type,
    Histology = meta_ref$histology,
    col = list(
      `Sample Type` = type_cols,
      Histology = hist_cols
    ),
    annotation_name_side = "left",
    annotation_legend_param = list(
      `Sample Type` = list(title = "Sample Type", direction = "horizontal"),
      Histology = list(title = "Histology", direction = "horizontal")
    )
  )
  
  # Define clustering
  row_hclust <- hclust(dist(mat, method = "euclidean"), method = "ward.D")
  col_hclust <- hclust(dist(t(mat), method = "euclidean"), method = "ward.D")
  
  ht_clustered <- Heatmap(
    mat,
    name = "Z-score",
    col = col_fun,
    cluster_rows = row_hclust,
    cluster_columns = col_hclust,
    top_annotation = ha_col,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 7),
    column_title = glue::glue("{method} Immune Deconvolution"),
    heatmap_legend_param = list(title = "Z-score", legend_direction = "horizontal")
  )
  
  pdf_clustered <- file.path(plot_dir, paste0(method, "_scaled_heatmap_clustered.pdf"))
  pdf(pdf_clustered, width = 8, height = 6)
  draw(ht_clustered, merge_legend = TRUE,
       heatmap_legend_side = "right", annotation_legend_side = "bottom")
  dev.off()
  message(paste0("✅ Saved clustered heatmap: ", pdf_clustered))
  
}



