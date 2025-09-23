# Create a base histologies file for Haydar miRNA project

# Author: Bicna Song

# Load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(rprojroot)
  library(stringr)
})

# Set path to module and results directories
root_dir <- find_root(has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "human-rna-expression")
input_dir <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")

# Load inputs
kw <- readr::read_tsv(file.path(input_dir, "immune_keywords.tsv"), show_col_types = FALSE)
gsva_scores <- readr::read_tsv(
  file.path(results_dir, "gobp-gsva-scores.tsv"),
  show_col_types = FALSE
)

# Normalization helper (lowercase, non-alnum -> space, collapse spaces)
norm_txt <- function(x) {
  x |>
    tolower() |>
    gsub("[^a-z0-9]+", " ", x = _, perl = TRUE) |>
    stringr::str_squish()
}

# Build per-category regex (on normalized keywords)
kw_clean <- kw |>
  mutate(k_norm = norm_txt(keyword)) |>
  filter(k_norm != "")

cat_patterns <- kw_clean |>
  group_by(category) |>
  summarise(
    # turn "alt comp pathway" into "\b(alt\s+comp\s+pathway|...)\b"
    pattern = paste0(
      "\\b(",
      paste(
        unique(stringr::str_replace_all(k_norm, " ", "\\\\s+")),
        collapse = "|"
      ),
      ")\\b"
    ),
    .groups = "drop"
  )

# Match categories to each Pathway (vectorized, no rowwise)
gsva_norm <- gsva_scores |>
  mutate(.path_norm = norm_txt(Pathway))

# For each category pattern, compute a logical vector of matches across all pathways
match_mat <- sapply(
  seq_len(nrow(cat_patterns)),
  function(i) stringr::str_detect(gsva_norm$.path_norm,
                                  stringr::regex(cat_patterns$pattern[i], ignore_case = TRUE))
)
if (!is.matrix(match_mat)) {
  # handle edge case: only one category -> sapply returns vector
  match_mat <- matrix(match_mat, ncol = 1)
}
colnames(match_mat) <- cat_patterns$category

# Collapse matches per row into a semicolon-separated category string (or NA)
cat_col <- apply(match_mat, 1, function(x) {
  hits <- colnames(match_mat)[x]
  if (length(hits) == 0) NA_character_ else paste(unique(hits), collapse = "; ")
})

# Add `category` and filter to immune-related pathways
immune_go <- gsva_norm |>
  mutate(category = cat_col) |>
  select(-.path_norm) |>
  filter(!is.na(category))

# write out
readr::write_tsv(immune_go, file.path(results_dir, "immune_filtered_gsva_scores.tsv"))

# `immune_go` now contains only pathways matching your immune keywords,
# with a `category` column listing the matched categories.
