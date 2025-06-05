# Create a base histologies file for Haydar miRNA project

# Author: Bicna Song

# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(rprojroot)
})

# Set path to module and results directories
root_dir <- "haydar-mirna/"
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "histology-preprocessing")
input_dir <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")

# Read and adjust sample_metadata
sample_metadata <- read_delim(
  file.path(data_dir, "miRNA-sample-metadata.txt"),
  delim = "\t",
  col_types = cols()
)

sample_metadata$sample_id[51] <- "4-1958-CC1"
sample_metadata$sample_id[52] <- "2-1578-CC1"
sample_metadata$sample_id[53] <- "1-1855-CC1"
sample_metadata$sample_id[54] <- "3-1890-CC1"
sample_metadata$sample_id[55] <- "5-1234-left-pons"

# Subset only the columns needed for merging
metadata_sub <- sample_metadata %>%
  select(sample_id, sample_type, histology, primary_site)

# Read RNA-seq and miRNA-seq manifests
rna_manifest <- read_tsv(
  file.path(input_dir, "Haydar-RNAseq-manifest_IDs_assigned.tsv"),
  col_types = cols()
)
mirna_manifest <- read_tsv(
  file.path(input_dir, "Haydar-miRNAseq-manifest_IDs_assigned.tsv"),
  col_types = cols()
)

# Combine manifests and deduplicate by Bioassay_ID
merged_manifest <- bind_rows(rna_manifest, mirna_manifest) %>%
  distinct(Bioassay_ID, .keep_all = TRUE)

# Keep only the key columns needed for downstream merging
merged_manifest <- merged_manifest %>%
  select(
    Bioassay_ID,
    experimental_strategy,
    sample_id,
    external_patient_id,
    external_sample_id
  )

# Identify and swap Bioassay_IDs for 2058-N / 2058-T within each strategy

# RNA‐seq pair:
rna_2058N_id <- merged_manifest %>%
  filter(experimental_strategy == "RNA-Seq", external_sample_id == "2058-N") %>%
  pull(Bioassay_ID)

rna_2058T_id <- merged_manifest %>%
  filter(experimental_strategy == "RNA-Seq", external_sample_id == "2058-T") %>%
  pull(Bioassay_ID)

# miRNA‐seq pair:
mirna_2058N_id <- merged_manifest %>%
  filter(experimental_strategy == "miRNA-Seq", external_sample_id == "2058-N") %>%
  pull(Bioassay_ID)

mirna_2058T_id <- merged_manifest %>%
  filter(experimental_strategy == "miRNA-Seq", external_sample_id == "2058-T") %>%
  pull(Bioassay_ID)

# Perform swapping within each cohort
merged_manifest <- merged_manifest %>%
  mutate(
    Bioassay_ID = case_when(
      experimental_strategy == "RNA-Seq"   & external_sample_id == "2058-N" ~ rna_2058T_id,
      experimental_strategy == "RNA-Seq"   & external_sample_id == "2058-T" ~ rna_2058N_id,
      experimental_strategy == "miRNA-Seq" & external_sample_id == "2058-N" ~ mirna_2058T_id,
      experimental_strategy == "miRNA-Seq" & external_sample_id == "2058-T" ~ mirna_2058N_id,
      TRUE ~ Bioassay_ID
    )
  )

# Join histology-related metadata onto merged_manifest
merged_hist <- merged_manifest %>%
  left_join(
    metadata_sub,
    by = c("external_sample_id" = "sample_id")
  ) %>%
  mutate(
    histology = na_if(histology, "N/A")
  ) %>%
  select(
    Bioassay_ID,
    experimental_strategy,
    sample_id,
    external_patient_id,
    external_sample_id,
    sample_type,
    histology,
    primary_site,
    everything()
  )

# Write out the final CSV
write_csv(
  merged_hist,
  file.path(results_dir, "histologies.csv")
)




