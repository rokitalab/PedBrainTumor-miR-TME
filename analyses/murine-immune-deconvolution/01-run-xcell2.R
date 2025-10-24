# Perform immune deconvolution of mouse RNA-seq data

# Author: Bicna Song

# Load libraries
suppressPackageStartupMessages({
  library(rprojroot)
  library(xCell2)
})

# Set path to module and results directories
root_dir <- find_root(has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "murine-immune-deconvolution")
results_dir <- file.path(analysis_dir, "results")

# Load expression matrix
mouse_rna <- readRDS(file.path(data_dir, "mouse-gene-expression-rsem-tpm-collapsed.all.rds"))
genes_mix <- rownames(mouse_rna)

# Define function for xCell2 analysis
run_xcell2 <- function(ref_name, ref_obj, mix_data, min_shared, output_file) {
  genes_ref <- getGenesUsed(ref_obj)
  overlap <- round(length(intersect(genes_mix, genes_ref)) / length(genes_ref) * 100, 2)
  message(paste0("[", ref_name, "] Overlap: ", overlap, "%"))
  
  res <- xCell2Analysis(mix = mix_data, xcell2object = ref_obj, minSharedGenes = min_shared)
  res_df <- as.data.frame(res)
  res_df <- tibble::rownames_to_column(res_df, var = "CellType")
  res_df$Reference <- ref_name
  
  write.table(res_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  return(res_df)
}

# Run analyses for all reference datasets
data("MouseRNAseqData.xCell2Ref", package = "xCell2")
res_mouse <- run_xcell2("MouseRNAseqData", MouseRNAseqData.xCell2Ref, mouse_rna, 0.7, file.path(results_dir, "MouseRNAseqData_xCell2_results.tsv"))

data("TabulaMurisBlood.xCell2Ref", package = "xCell2")
res_tm <- run_xcell2("TabulaMurisBlood", TabulaMurisBlood.xCell2Ref, mouse_rna, 0.8, file.path(results_dir, "TabulaMurisBlood_xCell2_results.tsv"))

data("ImmGenData.xCell2Ref", package = "xCell2")
res_immgen <- run_xcell2("ImmGenData", ImmGenData.xCell2Ref, mouse_rna, 0.8, file.path(results_dir, "ImmGenData_xCell2_results.tsv"))

# Merge all results into one summary table
merged_results <- dplyr::bind_rows(res_mouse, res_tm, res_immgen)
write.table(merged_results, file = file.path(results_dir, "xCell2_all_references_merged.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("✅ All analyses complete. Results saved as individual .tsv files and combined summary table: xCell2_all_references_merged.tsv")