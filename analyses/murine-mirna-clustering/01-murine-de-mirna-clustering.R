# Cluster DE miRNAs from CAR-T versus untreated mice by expression
# Author: Bicna Song | Oct 2025 
  
### Setup
suppressPackageStartupMessages({
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
  library(glue)
})

### Paths
root_dir <- find_root(has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "murine-mirna-clustering")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

mirna_cts_file <- file.path(data_dir, "30-963755216-miRNA_expression.csv")
histology_file <- file.path(root_dir, "analyses", 
                            "murine-mirna-differential-expression", 
                            "results", "sample-metadata.tsv")

### Load & preprocess
mirna_cts <- read_csv(mirna_cts_file) %>%
  column_to_rownames("...1") %>% 
  dplyr::select(-`Day14-Tube1`) %>%
  select(-matches("Day28"))
  
names(mirna_cts) <- str_replace(names(mirna_cts), "RESUB-", "")
mirna_cts <- mirna_cts[rowSums(mirna_cts != 0) >= 5, ]

### TPM calculation

# Assume fixed miRNA length of 22 bases (0.022 kb)
mirna_length_kb <- 22 / 1000  

# Function to compute TPMs
calculate_tpm <- function(count_matrix, gene_length_kb) {
  rpk <- count_matrix / gene_length_kb
  tpm <- sweep(rpk, 2, colSums(rpk), FUN = "/") * 1e6
  return(tpm)
}

# calculate tpm
mirna_tpm <- calculate_tpm(mirna_cts, gene_length_kb = mirna_length_kb)

# save TPM file to output
saveRDS(mirna_tpm,
        file.path(results_dir,
                  "murine-mirna-tpm.rds"))

mirna_tpm <- log2(as.matrix(mirna_tpm) + 1)

### Metadata
histology_df <- read_delim(histology_file, show_col_types = FALSE) %>%
  filter(time != "Day28")

### Normalize (Z-score)
miRNA_zscores <- t(scale(t(mirna_tpm)))
miRNA_zscores <- miRNA_zscores[complete.cases(miRNA_zscores), ]

### Load DE miRNAs
DE_miRNA_list <- read_tsv(file.path(root_dir, "analyses", "murine-mirna-differential-expression", "results",
                                    "mirna-differential-expression-deseq2-b7h3-stop-vs-untreated-by-timepoint.tsv"),
                          show_col_types = FALSE)

sig_DE_miRNA_list <- DE_miRNA_list %>%
  mutate(expr_pattern = str_replace_all(expr_pattern, 
                                        c("Day 14" = "Day14",
                                          "Day 21" = "Day21",
                                          "Day 28" = "Day28"))) %>%
  filter(expr_pattern %in% c("B7H3 down, Day14","B7H3 down, Day21", "B7H3 down, multiple",
                             "B7H3 up, Day14", "B7H3 up, Day21", "B7H3 up, multiple")) %>%
  filter(mirna_id != "NovelmiRNA-1146")

write_csv(sig_DE_miRNA_list, file.path(results_dir, "mouse_sig_DE_miRNA_list.csv"))

de_mirnas <- sig_DE_miRNA_list$mirna_id
miRNA_expr_sub <- miRNA_zscores[rownames(miRNA_zscores) %in% de_mirnas, , drop=FALSE] 
  
### Annotations

# Expression pattern
expr_cols <- c(
  # Downregulated (blue shades)
  "B7H3 down, Day14"   = "#08306B",  # dark navy blue
  "B7H3 down, Day21"   = "#2171B5",  # medium blue
  "B7H3 down, multiple" = "#6BAED6",  # light blue
  
  # Upregulated (red shades)
  "B7H3 up, Day14"     = "#67000D",  # dark crimson
  "B7H3 up, Day21"     = "#CB181D",  # strong red
  "B7H3 up, multiple"   = "#FB6A4A"   # salmon red
)

# Time
time_cols <- c(
  "Day14" = "#A67C52",  # warm taupe
  "Day21" = "#FDB462"   # muted orange
)

# Treatment
treatment_cols <- c(
  "B7H3"      = "#E64B35FF",
  "STOP"      = "#4DBBD5FF",
  "untreated" = "#00A087FF"
)

# Annotation
anno_cols <- list(
  Annotated = c("Yes" = "#1B9E77", "No" = "#D95F02"),  # green / orange
  `Expression pattern` = expr_cols,
  Treatment = treatment_cols,
  Time = time_cols
)

mirna_anno <- data.frame(mirna_id = rownames(miRNA_expr_sub)) %>%
  left_join(sig_DE_miRNA_list, by = "mirna_id") %>%
  mutate(
    Annotated = if_else(grepl("mmu", mirna_id), "Yes", "No"),
    Annotated = factor(Annotated, levels = c("Yes", "No"))
  ) %>%
  distinct(mirna_id, .keep_all = TRUE)

ra_left <- rowAnnotation(
  Annotated = mirna_anno$Annotated,
  `Expression pattern` = mirna_anno$expr_pattern,
  col = anno_cols[c("Annotated", "Expression pattern")],
  annotation_name_gp = gpar(fontsize = 10)
)

ha_all <- HeatmapAnnotation(
  Time = histology_df$time,
  Treatment = histology_df$treatment,
  col = anno_cols[c("Time", "Treatment")],
  annotation_name_gp = gpar(fontsize = 10)
)

col_fun <- colorRamp2(c(-4, 0, 4), c("navyblue", "white", "orangered"))

### Heatmap
set.seed(123)
  
mirna_ht <- Heatmap(
  miRNA_expr_sub,
  name = "TPM z-score",
  col = col_fun,
  na_col = "gray",
  clustering_distance_rows = "spearman",
  clustering_distance_columns = "spearman",
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  top_annotation = ha_all,
  left_annotation = ra_left,
  use_raster = TRUE,
  heatmap_legend_param = list(legend_gp = gpar(fontsize = 10))
)
 
pdf(file.path(plot_dir, glue::glue("de-mirna-heatmap.pdf")), width = 15, height = 12)
mirna_ht <- draw(mirna_ht)
dev.off()
 
### Per-time point heatmaps
# Generate one heatmap per time point (e.g., Day 14, Day 21) using only the miRNAs significantly associated with that condition (up/down)
# Iterate for each time point
for (day in c("Day14", "Day21")) {
  
  message("Generating heatmap for ", day)
  
  # Filter miRNAs specific to this time point
  de_mirna_list_time <- sig_DE_miRNA_list %>%
    filter(expr_pattern %in% c(
      paste("B7H3 down,", day),
      "B7H3 down, multiple",
      paste("B7H3 up,", day),
      "B7H3 up, multiple"
    ))
  
  de_mirnas <- de_mirna_list_time$mirna_id
  miRNA_expr_sub <- miRNA_zscores[rownames(miRNA_zscores) %in% de_mirnas, , drop = FALSE]
  
  # Subset columns by time point
  miRNA_expr_sub <- miRNA_expr_sub[, grepl(day, colnames(miRNA_expr_sub)), drop = FALSE]
  
  # Row annotation (miRNA-level)
  mirna_anno <- data.frame(mirna_id = rownames(miRNA_expr_sub)) %>%
    left_join(de_mirna_list_time, by = "mirna_id") %>%
    mutate(
      Annotated = if_else(grepl("mmu", mirna_id), "Yes", "No"),
      Annotated = factor(Annotated, levels = c("Yes", "No"))
    ) %>%
    distinct(mirna_id, .keep_all = TRUE)
  
  ra_left <- rowAnnotation(
    Annotated = mirna_anno$Annotated,
    `Expression pattern` = mirna_anno$expr_pattern,
    col = anno_cols[c("Annotated", "Expression pattern")],
    annotation_name_gp = gpar(fontsize = 10)
  )
  
  # Column annotation (sample-level)
  histology_df_sub <- histology_df %>%
    filter(time == day)
  
  ha_all <- HeatmapAnnotation(
    Time = histology_df_sub$time,
    Treatment = histology_df_sub$treatment,
    col = anno_cols[c("Time", "Treatment")],
    annotation_name_gp = gpar(fontsize = 10)
  )
  
  # Color scale
  col_fun <- colorRamp2(c(-4, 0, 4), c("navyblue", "white", "orangered"))
  
  # Generate heatmap
  set.seed(123)
  mirna_ht <- Heatmap(
    miRNA_expr_sub,
    name = "TPM z-score",
    col = col_fun,
    na_col = "gray",
    clustering_distance_rows = "spearman",
    clustering_distance_columns = "spearman",
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    top_annotation = ha_all,
    left_annotation = ra_left,
    use_raster = TRUE,
    heatmap_legend_param = list(legend_gp = gpar(fontsize = 10))
  )
  
  # Save to PDF
  pdf(file.path(plot_dir, glue::glue("de-mirna-heatmap-{day}.pdf")),
      width = 15, height = 12)
  draw(mirna_ht)
  dev.off()
}

### Session Info

sessionInfo()
