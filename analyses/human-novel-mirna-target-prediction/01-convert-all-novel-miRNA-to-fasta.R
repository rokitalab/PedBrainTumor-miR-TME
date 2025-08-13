# Read a list of rna table files and convert to fasta
library(readr)
library(dplyr)
library(Biostrings)


root_dir     <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir     <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "human-novel-mirna-target-prediction")
results_dir  <- file.path(analysis_dir, "results")
output_fasta <- file.path(results_dir, "combined_all_novel_mirna.fa")

batch1_file  <- file.path(root_dir, "analyses", "human-mirna-expression", "results", "30-931106737-all_novel_miRNA_merged_id.tsv")
batch2_file  <- file.path(root_dir, "analyses", "human-mirna-expression", "results", "30-1075661268-all_novel_miRNA_merged_id.tsv")

# create options for input file dir and output_fasta

# Read input files

file_list <- list.files(path = input_dir, pattern = "*.tsv", full.names = TRUE)
input_dfs <- lapply(file.list, read.delim, simplify = FALSE)

# Stack rows and keep the first occurrence for each merge_id
mirna_df <- bind_rows(input_dfs) %>%
  distinct(merged_id, .keep_all = TRUE) %>%
  select(
    merged_id,
    mature_seq,
    star_seq,
    precursor_seq,
    precursor_coordinate
  )

# Write Merged Table to TSV

# Save the combined data frame as TSV
write_tsv(
  mirna_df,
  file = file.path(results_dir, "combined_all_novel_mirna.tsv"),
  col_names = FALSE
)

# Convert to FASTA format using mature miRNA sequence

mature_mirna <- RNAStringSet(mirna_df$mature_seq)
names(mature_mirna) <- mirna_df$merged_id

# Write FASTA output

writeXStringSet(mature_mirna, filepath = output_fasta)