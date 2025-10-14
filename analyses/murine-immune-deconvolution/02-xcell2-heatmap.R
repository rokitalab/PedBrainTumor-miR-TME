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
analysis_dir <- file.path(root_dir, "analyses", "murine-immune-deconvolution")
input_dir <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

# File path
merged_file <- file.path(results_dir, "xCell2_all_references_merged.tsv")
metadata_file <- file.path(root_dir, "analyses", "murine-rna-expression", "input", "Haydar_Mouse_RNA_miRNA_manifest_IDs_assigned.tsv")

# Read merged results
merged_results <- read_tsv(merged_file)
metadata_df <- read_delim(metadata_file) %>%
  filter(experimental_strategy == "RNA-Seq") %>%
  select(external_sample_id, Day, Treatment, Bioassay_ID) %>%
  distinct(external_sample_id, .keep_all = TRUE) %>%
  rename(time = Day) %>%
  rename(treatment = Treatment) %>%
  mutate(
    time = paste0("Day", time),
    treatment = str_remove(treatment, "\\s*CAR\\b"),
    treatment = str_squish(treatment),
    treatment = if_else(treatment == "B7H3 STOP", "STOP", treatment)
  ) 

# Pivot to long format for scaling
scaled_results <- merged_results %>%
  pivot_longer(
    cols = -c(CellType, Reference),
    names_to = "Bioassay_ID",
    values_to = "Score"
  ) %>%
  group_by(Reference, CellType) %>%
  mutate(
    Scaled_Score = if (sd(Score, na.rm = TRUE) == 0) 0 
    else scale(Score)[, 1]
  ) %>%
  ungroup()

# Pivot back to wide format if needed (for plotting heatmap or saving)
scaled_results_wide <- scaled_results %>%
  select(CellType, Bioassay_ID, Scaled_Score, Reference) %>%
  pivot_wider(
    names_from = Bioassay_ID,
    values_from = Scaled_Score
  )

# Save scaled results
write.table(scaled_results_wide, file.path(results_dir, "xCell2_scaled_within_reference.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("✅ Scaled enrichment scores (z-scores) saved to: xCell2_scaled_within_reference.tsv")

# Generate heatmaps per reference
unique_refs <- unique(scaled_results$Reference)
col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))

# Custom metadata colors
time_cols <- c(
  "Day14" = "#A67C52",  # warm taupe
  "Day21" = "#FDB462"   # muted orange
)

treatment_cols <- c(
  "B7H3"      = "#E64B35FF",
  "STOP"      = "#4DBBD5FF",
  "Untreated" = "#00A087FF"
)

for (ref in unique_refs) {
  message(paste0("Generating ComplexHeatmap for reference: ", ref))
  
  mat <- scaled_results %>%
    filter(Reference == ref) %>%
    select(CellType, Bioassay_ID, Scaled_Score) %>%
    pivot_wider(names_from = Bioassay_ID, values_from = Scaled_Score) %>%
    column_to_rownames("CellType") %>%
    as.matrix()
  
  # Remove rows with zero variance
  non_constant <- apply(mat, 1, function(x) sd(x, na.rm = TRUE) > 0)
  mat <- mat[non_constant, , drop = FALSE]
  
  if (nrow(mat) < 2 || ncol(mat) < 2) {
    message(paste0("⚠️ Skipping ", ref, ": insufficient variable data for clustering."))
    next
  }
  
  # Match metadata for available samples
  meta_ref <- metadata_df %>%
    filter(Bioassay_ID %in% colnames(mat)) %>%
    column_to_rownames("Bioassay_ID")
  
  # Ensure column order matches matrix columns
  meta_ref <- meta_ref[colnames(mat), , drop = FALSE]
  
  # Define annotations
  ha_col <- HeatmapAnnotation(
    Time = meta_ref$time,
    Treatment = meta_ref$treatment,
    col = list(
      Time = time_cols,
      Treatment = treatment_cols
    ),
    annotation_name_side = "left",
    annotation_legend_param = list(
      Time = list(title = "Time", direction = "horizontal"),
      Treatment = list(title = "Treatment", direction = "horizontal")
    )
  )
  
  # Define clustering
  row_hclust <- hclust(dist(mat, method = "euclidean"), method = "ward.D2")
  col_hclust <- hclust(dist(t(mat), method = "euclidean"), method = "ward.D2")
  
  # Create heatmap
  ht <- Heatmap(
    mat,
    name = "Z-score",
    col = col_fun,
    cluster_rows = row_hclust,
    cluster_columns = col_hclust,
    top_annotation = ha_col,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 7),
    column_title = paste(ref, "xCell2 (Scaled)"),
    heatmap_legend_param = list(
      title = "Z-score",
      legend_direction = "horizontal"
    )
  )
  
  # Save to PDF
  pdf_file <- file.path(plot_dir, paste0("xCell2_scaled_heatmap_", ref, ".pdf"))
  pdf(pdf_file, width = 11, height = 8)
  draw(ht, merge_legend = TRUE, heatmap_legend_side = "right", annotation_legend_side = "bottom")
  dev.off()
  
  message(paste0("✅ Saved annotated ComplexHeatmap to: ", pdf_file))
}



