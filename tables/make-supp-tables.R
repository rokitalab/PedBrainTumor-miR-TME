#!/usr/bin/env Rscript

# Generate the supplementary tables for the manuscript.

# Find the repository root from the current working directory without requiring
# an additional package.
find_repo_root <- function(path = getwd()) {
  path <- normalizePath(path, winslash = "/", mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(path, ".git"))) {
      return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not find the repository root (expected a .git directory).")
    }
    path <- parent
  }
}

root_dir <- find_repo_root()
tables_dir <- file.path(root_dir, "tables")
analysis_dir <- file.path(root_dir, "analyses", "histology-preprocessing")
results_dir <- file.path(analysis_dir, "results")
input_dir <- file.path(analysis_dir, "input")
tables_input_dir <- file.path(tables_dir, "input")
mirna_results_dir <- file.path(
  root_dir, "analyses", "human-mirna-expression", "results"
)
clustering_results_dir <- file.path(
  root_dir, "analyses", "mirna-clustering", "results"
)
human_mirna_go_results_dir <- file.path(
  root_dir, "analyses", "human-mirna-go-enrichment", "results"
)
human_rna_results_dir <- file.path(
  root_dir, "analyses", "human-rna-expression", "results"
)
murine_rna_input_dir <- file.path(
  root_dir, "analyses", "murine-rna-expression", "input"
)

histologies_path <- file.path(results_dir, "histologies.tsv")
metadata_path <- file.path(input_dir, "histologies.tsv")
cbtn_files <- list.files(
  tables_input_dir,
  pattern = "^cbtn_all_.*\\.csv$",
  full.names = TRUE
)
if (length(cbtn_files) != 1L) {
  stop("Expected one cbtn_all_*.csv file in ", tables_input_dir,
       "; found ", length(cbtn_files), ".")
}

read_tsv <- function(path) {
  read.delim(
    path,
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = "NA"
  )
}

histologies <- read_tsv(histologies_path)
metadata <- read_tsv(metadata_path)
cbtn_data <- read.csv(
  cbtn_files[[1L]],
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)
metadata_columns <- c(
  "sample_id", "reported_gender", "race", "tumor_descriptor", "CNS_region"
)

excluded_histologies <- c("HGG", "GNT", "ATRT")
histologies <- histologies[
  is.na(histologies$histology) | !histologies$histology %in% excluded_histologies,
  ,
  drop = FALSE
]

# Omit 1-1855-CC1 because it is a duplicate of sample 1855-T.
duplicate_sample <- !is.na(histologies$external_sample_id) &
  histologies$external_sample_id == "1-1855-CC1"
histologies <- histologies[!duplicate_sample, , drop = FALSE]

# The sample ID is missing from the results file for this known sample.
histologies$sample_id[histologies$external_sample_id == "1516-T"] <- "7316-6908"

missing_columns <- setdiff(metadata_columns, names(metadata))
if (length(missing_columns) > 0L) {
  stop("Missing required columns in input histologies.tsv: ",
       paste(missing_columns, collapse = ", "))
}
if (!"sample_id" %in% names(histologies)) {
  stop("Missing required sample_id column in result histologies.tsv.")
}
cbtn_columns <- c(
  "cbtn_specimen_group_id", "legal_sex", "race", "event_type",
  "cns_integrated_diagnosis"
)
missing_cbtn_columns <- setdiff(cbtn_columns, names(cbtn_data))
if (length(missing_cbtn_columns) > 0L) {
  stop("Missing required columns in CBTN CSV: ",
       paste(missing_cbtn_columns, collapse = ", "))
}

fields_to_append <- setdiff(metadata_columns, "sample_id")
if (any(fields_to_append %in% names(histologies))) {
  stop("One or more metadata columns already exist in result histologies.tsv.")
}

# Prefer the histology input file. Use the CBTN CSV only when a table sample ID
# does not occur in that primary metadata file.
matched_ids <- unique(histologies$sample_id[!is.na(histologies$sample_id)])
primary_ids <- unique(metadata$sample_id[!is.na(metadata$sample_id)])
metadata_matches <- metadata[
  !is.na(metadata$sample_id) & metadata$sample_id %in% matched_ids,
  metadata_columns,
  drop = FALSE
]

secondary_search_ids <- setdiff(matched_ids, primary_ids)
# Map the corresponding CBTN fields to the names used in Table S1.
cbtn_metadata <- data.frame(
  sample_id = cbtn_data$cbtn_specimen_group_id,
  reported_gender = cbtn_data$legal_sex,
  race = cbtn_data$race,
  tumor_descriptor = cbtn_data$event_type,
  CNS_region = NA_character_,
  stringsAsFactors = FALSE
)
cbtn_matches <- cbtn_metadata[
  !is.na(cbtn_metadata$sample_id) &
    cbtn_metadata$sample_id %in% secondary_search_ids,
  metadata_columns,
  drop = FALSE
]
secondary_only_ids <- unique(cbtn_matches$sample_id)
metadata_matches <- rbind(metadata_matches, cbtn_matches)

# The metadata files can contain many records per sample. Confirm that records
# matching this table agree on non-missing values before reducing to one row.
metadata_by_id <- split(metadata_matches, metadata_matches$sample_id)
conflicting_ids <- character()
for (column in fields_to_append) {
  ids_with_conflicts <- names(metadata_by_id)[vapply(
    metadata_by_id,
    function(rows) {
      values <- unique(rows[[column]][!is.na(rows[[column]])])
      length(values) > 1L
    },
    logical(1)
  )]
  conflicting_ids <- union(conflicting_ids, ids_with_conflicts)
}
if (length(conflicting_ids) > 0L) {
  stop("Conflicting requested metadata for sample_id(s): ",
       paste(conflicting_ids, collapse = ", "))
}

metadata <- do.call(rbind, lapply(metadata_by_id, function(rows) {
  all_rows <- rows
  rows <- rows[1L, metadata_columns, drop = FALSE]
  for (column in fields_to_append) {
    values <- unique(all_rows[[column]][!is.na(all_rows[[column]])])
    if (length(values) > 0L) {
      rows[[column]] <- values[[1L]]
    }
  }
  rows
}))
metadata_index <- match(histologies$sample_id, metadata$sample_id)
for (column in fields_to_append) {
  histologies[[column]] <- metadata[[column]][metadata_index]
}

# Use the broad CNS region category only when the detailed primary site is
# missing from the results file.
missing_primary_site <- is.na(histologies$primary_site) |
  !nzchar(trimws(histologies$primary_site))
has_cns_region <- !is.na(histologies$CNS_region) &
  nzchar(trimws(histologies$CNS_region))
use_cns_region <- missing_primary_site & has_cns_region
histologies$primary_site[use_cns_region] <-
  histologies$CNS_region[use_cns_region]

# CBTN's integrated diagnosis specifies Group 4 for this medulloblastoma even
# though its molecular subtype field is unavailable.
cbtn_integrated_diagnosis <- cbtn_data$cns_integrated_diagnosis[
  match(histologies$sample_id, cbtn_data$cbtn_specimen_group_id)
]
missing_molecular_subtype <- is.na(histologies$molecular_subtype) |
  !nzchar(trimws(histologies$molecular_subtype))
group4_diagnosis <- !is.na(cbtn_integrated_diagnosis) &
  cbtn_integrated_diagnosis == "Medulloblastoma, Group 4"
use_group4_subtype <- missing_molecular_subtype & group4_diagnosis &
  histologies$histology == "MB"
histologies$molecular_subtype[use_group4_subtype] <- "MB, Group4"

# Manually supplied metadata for the 2248 tumor/adjacent-normal pair.
pair_rows <- histologies$external_sample_id %in% c("2248-T", "2248-N")
histologies$reported_gender[pair_rows] <- "Male"
histologies$race[pair_rows] <- "White"
histologies$tumor_descriptor[pair_rows] <- "Deceased"
healthy_brain_rows <- !is.na(histologies$external_sample_id) &
  grepl("^B[1-5]$", histologies$external_sample_id)
histologies$tumor_descriptor[healthy_brain_rows] <- "Deceased"
histologies$molecular_subtype[
  histologies$external_sample_id == "2248-T"
] <- "DMG, H3 K28"

# Collapse RNA-seq and miRNA-seq assay records to one row per external sample.
source_row_count <- nrow(histologies)
matched_source_rows <- sum(!is.na(metadata_index))
sample_columns <- setdiff(
  names(histologies),
  c("Bioassay_ID", "experimental_strategy", "pathology_diagnosis")
)
sample_groups <- local({
  sample_keys <- as.character(histologies$external_sample_id)
  missing_external_id <- is.na(sample_keys) | !nzchar(sample_keys)
  sample_keys[missing_external_id] <- paste0(
    "missing_external_sample_id_row_", which(missing_external_id)
  )
  split(
    seq_len(nrow(histologies)),
    factor(sample_keys, levels = unique(sample_keys))
  )
})

histologies <- do.call(rbind, lapply(sample_groups, function(row_indices) {
  assay_rows <- histologies[row_indices, , drop = FALSE]
  sample_row <- assay_rows[1L, sample_columns, drop = FALSE]

  for (column in sample_columns) {
    values <- unique(assay_rows[[column]][!is.na(assay_rows[[column]])])
    if (length(values) > 1L) {
      stop("Conflicting ", column, " values for external_sample_id ",
           assay_rows$external_sample_id[[1L]])
    }
    if (length(values) == 1L) {
      sample_row[[column]] <- values[[1L]]
    }
  }

  sample_row$RNA_seq <- if (any(assay_rows$experimental_strategy == "RNA-Seq")) {
    "Yes"
  } else {
    "No"
  }
  sample_row$miRNA_seq <- if (any(assay_rows$experimental_strategy == "miRNA-Seq")) {
    "Yes"
  } else {
    "No"
  }
  sample_row
}))

output_columns <- c(
  external_sample_id = "Sample ID",
  sample_id = "CBTN Sample ID",
  sample_type = "Sample Type",
  primary_site = "CNS Region",
  cancer_group = "Cancer Group",
  molecular_subtype = "Molecular Subtype",
  histology = "Patient Histology Group",
  reported_gender = "Reported Gender",
  race = "Reported Race",
  tumor_descriptor = "Collection Event",
  RNA_seq = "RNA-Seq?",
  miRNA_seq = "miRNA-Seq?"
)
missing_output_columns <- setdiff(names(output_columns), names(histologies))
if (length(missing_output_columns) > 0L) {
  stop("Missing expected output columns: ",
       paste(missing_output_columns, collapse = ", "))
}
histologies <- histologies[, names(output_columns), drop = FALSE]
names(histologies) <- unname(output_columns)
histologies <- histologies[
  order(histologies[["Patient Histology Group"]], na.last = TRUE),
  ,
  drop = FALSE
]

xlsx_path <- file.path(tables_dir, "TableS1.xlsx")
write_xlsx <- function(sheets, path, na_rep = "") {
  if (is.data.frame(sheets)) {
    sheets <- list(Sheet1 = sheets)
  }
  if (is.null(names(sheets)) || any(!nzchar(names(sheets)))) {
    stop("Every Excel worksheet must have a name.")
  }
  xml_escape <- function(value) {
    value <- gsub("&", "&amp;", value, fixed = TRUE)
    value <- gsub("<", "&lt;", value, fixed = TRUE)
    gsub(">", "&gt;", value, fixed = TRUE)
  }

  excel_column <- function(index) {
    label <- ""
    while (index > 0L) {
      remainder <- (index - 1L) %% 26L
      label <- paste0(LETTERS[remainder + 1L], label)
      index <- (index - 1L) %/% 26L
    }
    label
  }

  cell_xml <- function(value, cell_ref, column_data = NULL) {
    if (is.na(value)) {
      if (!nzchar(na_rep)) {
        return(paste0('<c r="', cell_ref, '"/>'))
      }
      return(paste0(
        '<c r="', cell_ref, '" t="inlineStr"><is><t xml:space="preserve">',
        xml_escape(na_rep), "</t></is></c>"
      ))
    }
    if (is.numeric(column_data)) {
      return(paste0('<c r="', cell_ref, '" t="n"><v>',
                    as.character(value), "</v></c>"))
    }
    if (is.logical(column_data)) {
      return(paste0('<c r="', cell_ref, '" t="b"><v>',
                    as.integer(value), "</v></c>"))
    }
    paste0(
      '<c r="', cell_ref, '" t="inlineStr"><is><t xml:space="preserve">',
      xml_escape(as.character(value)), "</t></is></c>"
    )
  }

  xlsx_dir <- tempfile("supp-tables-")
  dir.create(xlsx_dir)
  dir.create(file.path(xlsx_dir, "_rels"))
  dir.create(file.path(xlsx_dir, "xl", "_rels"), recursive = TRUE)
  dir.create(file.path(xlsx_dir, "xl", "worksheets"), recursive = TRUE)

  write_xml <- function(relative_path, lines) {
    writeLines(lines, file.path(xlsx_dir, relative_path), useBytes = TRUE)
  }

  worksheet_paths <- paste0(
    "xl/worksheets/sheet", seq_along(sheets), ".xml"
  )
  worksheet_overrides <- vapply(seq_along(sheets), function(index) {
    paste0(
      '<Override PartName="/', worksheet_paths[[index]],
      '" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
    )
  }, character(1))
  sheet_entries <- vapply(seq_along(sheets), function(index) {
    paste0(
      '<sheet name="', xml_escape(names(sheets)[[index]]), '" sheetId="',
      index, '" r:id="rId', index, '"/>'
    )
  }, character(1))
  relationship_entries <- vapply(seq_along(sheets), function(index) {
    paste0(
      '<Relationship Id="rId', index,
      '" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet',
      index, '.xml"/>'
    )
  }, character(1))

  write_xml("[Content_Types].xml", c(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">',
    '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>',
    '<Default Extension="xml" ContentType="application/xml"/>',
    '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>',
    worksheet_overrides,
    "</Types>"
  ))
  write_xml("_rels/.rels", c(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">',
    '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>',
    "</Relationships>"
  ))
  write_xml("xl/workbook.xml", c(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">',
    "<sheets>",
    sheet_entries,
    "</sheets>",
    "</workbook>"
  ))
  write_xml("xl/_rels/workbook.xml.rels", c(
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
    '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">',
    relationship_entries,
    "</Relationships>"
  ))

  for (sheet_index in seq_along(sheets)) {
    data <- sheets[[sheet_index]]
    sheet_rows <- vapply(seq_len(nrow(data) + 1L), function(row_index) {
      cells <- vapply(seq_len(ncol(data)), function(column_index) {
        cell_ref <- paste0(excel_column(column_index), row_index)
        if (row_index == 1L) {
          cell_xml(names(data)[[column_index]], cell_ref)
        } else {
          cell_xml(
            data[[column_index]][row_index - 1L],
            cell_ref,
            column_data = data[[column_index]]
          )
        }
      }, character(1))
      paste0('<row r="', row_index, '">', paste(cells, collapse = ""), "</row>")
    }, character(1))

    last_column <- excel_column(ncol(data))
    write_xml(worksheet_paths[[sheet_index]], c(
      '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
      '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">',
      "<sheetData>",
      sheet_rows,
      "</sheetData>",
      paste0('<autoFilter ref="A1:', last_column, nrow(data) + 1L, '"/>'),
      "</worksheet>"
    ))
  }

  previous_dir <- getwd()
  on.exit(setwd(previous_dir), add = TRUE)
  on.exit(unlink(xlsx_dir, recursive = TRUE), add = TRUE)
  path <- file.path(normalizePath(dirname(path)), basename(path))
  setwd(xlsx_dir)
  if (file.exists(path)) {
    invisible(file.remove(path))
  }
  utils::zip(path, c(
    "[Content_Types].xml",
    "_rels/.rels",
    "xl/workbook.xml",
    "xl/_rels/workbook.xml.rels",
    worksheet_paths
  ))
}

write_xlsx(list(S1 = histologies), xlsx_path)

table_s2_files <- c(
  "DIPG or DMG-vs-healthyNormal" = "DESeq2_DIPG or DMG_vs_healthyNormal.csv",
  "EPN-vs-healthyNormal" = "DESeq2_EPN_vs_healthyNormal.csv",
  "LGG-vs-healthyNormal" = "DESeq2_LGG_vs_healthyNormal.csv",
  "MB-vs-healthyNormal" = "DESeq2_MB_vs_healthyNormal.csv",
  "DIPG or DMG-vs-adj-normal" = "DESeq2_DIPG or DMG_paired_full.csv",
  "EPN-vs-adj-normal" = "DESeq2_EPN_paired_full.csv",
  "MB-vs-adj-normal" = "DESeq2_MB_paired_full.csv"
)
table_s2_sheets <- lapply(table_s2_files, function(file_name) {
  file_path <- file.path(mirna_results_dir, file_name)
  if (!file.exists(file_path)) {
    stop("Missing miRNA differential expression file: ", file_path)
  }
  read.csv(
    file_path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
})
expected_s2_columns <- names(table_s2_sheets[[1L]])
if (!all(vapply(table_s2_sheets, function(sheet) {
  identical(names(sheet), expected_s2_columns)
}, logical(1)))) {
  stop("The Table S2 input files do not all have the same columns.")
}
table_s2_path <- file.path(tables_dir, "TableS2.xlsx")
write_xlsx(table_s2_sheets, table_s2_path)

table_s3_files <- c(
  "DIPG or DMG" = "DIPG or DMG-de-mirna-cluster-membership-immune-scores.tsv",
  "Medulloblastoma" = "MB-de-mirna-cluster-membership-immune-scores.tsv"
)
table_s3_sources <- lapply(table_s3_files, function(file_name) {
  file_path <- file.path(clustering_results_dir, file_name)
  if (!file.exists(file_path)) {
    stop("Missing miRNA clustering result file: ", file_path)
  }
  read_tsv(file_path)
})
expected_s3_columns <- names(table_s3_sources[[1L]])
if (!all(vapply(table_s3_sources, function(source) {
  identical(names(source), expected_s3_columns)
}, logical(1)))) {
  stop("The Table S3 input files do not all have the same columns.")
}

table_s3_sheets <- lapply(seq_along(table_s3_sources), function(index) {
  source <- table_s3_sources[[index]]
  names(source)[names(source) == "row_cluster"] <- "miRNA cluster"
  columns_after_first_four <- expected_s3_columns[5:length(expected_s3_columns)]
  for (column in columns_after_first_four) {
    names(source)[names(source) == column] <- paste0(column, "-correlation")
  }
  source
})
names(table_s3_sheets) <- names(table_s3_files)
table_s3_path <- file.path(tables_dir, "TableS3.xlsx")
write_xlsx(table_s3_sheets, table_s3_path, na_rep = "NA")
table_s3_row_count <- sum(vapply(table_s3_sheets, nrow, integer(1)))

message("Prepared Table S1 with ", nrow(histologies), " samples; ",
        matched_source_rows, " of ", source_row_count,
        " source assay rows matched to input metadata.")
message("Wrote ", xlsx_path, ".")
message("Wrote ", table_s2_path, " with ", length(table_s2_sheets),
        " worksheets.")
message("Wrote ", table_s3_path, " with ", length(table_s3_sheets),
        " worksheets and ", table_s3_row_count, " total rows.")
message(length(secondary_only_ids),
        " unique sample_id(s) were found only in the CBTN CSV.")

# Table S4: oncogenic and tumor-suppressive miRNAs.
table_s4_input_path <- file.path(tables_input_dir, "onco-ts-mirs.tsv")
if (!file.exists(table_s4_input_path)) {
  table_s4_input_path <- file.path(tables_input_dir, "onco-ts-mirs.txt")
}
if (!file.exists(table_s4_input_path)) {
  stop("Missing Table S4 input file: ", table_s4_input_path)
}
table_s4 <- read_tsv(table_s4_input_path)
table_s4_path <- file.path(tables_dir, "TableS4.xlsx")
write_xlsx(list(S4 = table_s4), table_s4_path)
message("Wrote ", table_s4_path, " with ", nrow(table_s4), " rows and ",
        ncol(table_s4), " columns.")

# Table S5: full immune-term miRNA target enrichment results.
table_s5_files <- c(
  "DIPG-or-DMG-cluster6" =
    "DIPG-or-DMG-cluster6-mirna-target-go-enr-immune-terms-full.tsv",
  "MB-cluster1" = "MB-cluster1-mirna-target-go-enr-immune-terms-full.tsv"
)
table_s5_sheets <- lapply(table_s5_files, function(file_name) {
  file_path <- file.path(human_mirna_go_results_dir, file_name)
  if (!file.exists(file_path)) {
    stop("Missing Table S5 input file: ", file_path)
  }
  read_tsv(file_path)
})
expected_s5_columns <- names(table_s5_sheets[[1L]])
if (!all(vapply(table_s5_sheets, function(sheet) {
  identical(names(sheet), expected_s5_columns)
}, logical(1)))) {
  stop("The Table S5 input files do not all have the same columns.")
}
table_s5_path <- file.path(tables_dir, "TableS5.xlsx")
write_xlsx(table_s5_sheets, table_s5_path)

# Table S6: RNA differential expression files corresponding to the Table S2
# miRNA comparisons. The paired tumor-versus-adjacent-normal comparisons use
# the paired_full result files, as they do for Table S2.
table_s6_files <- c(
  "DIPG or DMG-vs-healthyNormal" =
    "DESeq2_DIPG or DMG_vs_healthyNormal.csv",
  "EPN-vs-healthyNormal" = "DESeq2_EPN_vs_healthyNormal.csv",
  "LGG-vs-healthyNormal" = "DESeq2_LGG_vs_healthyNormal.csv",
  "MB-vs-healthyNormal" = "DESeq2_MB_vs_healthyNormal.csv",
  "DIPG or DMG-vs-adj-normal" = "DESeq2_DIPG or DMG_paired_full.csv",
  "EPN-vs-adj-normal" = "DESeq2_EPN_paired_full.csv",
  "MB-vs-adj-normal" = "DESeq2_MB_paired_full.csv"
)
table_s6_sheets <- lapply(table_s6_files, function(file_name) {
  file_path <- file.path(human_rna_results_dir, file_name)
  if (!file.exists(file_path)) {
    stop("Missing Table S6 input file: ", file_path)
  }
  read.csv(
    file_path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
})
expected_s6_columns <- names(table_s6_sheets[[1L]])
if (!all(vapply(table_s6_sheets, function(sheet) {
  identical(names(sheet), expected_s6_columns)
}, logical(1)))) {
  stop("The Table S6 input files do not all have the same columns.")
}
table_s6_path <- file.path(tables_dir, "TableS6.xlsx")
write_xlsx(table_s6_sheets, table_s6_path)

# Table S7: one row per sample matched by the Day/Treatment/aliquot key used
# in the murine miRNA-to-RNA correlation analysis. Keep all sequenced samples
# (including samples excluded from particular downstream correlation runs).
table_s7_input_path <- file.path(
  murine_rna_input_dir,
  "Haydar_Mouse_RNA_miRNA_manifest_IDs_assigned.tsv"
)
if (!file.exists(table_s7_input_path)) {
  stop("Missing Table S7 input manifest: ", table_s7_input_path)
}
mouse_manifest <- read_tsv(table_s7_input_path)
required_s7_columns <- c(
  "external_sample_id", "experimental_strategy", "RNA_library",
  "organism", "Day", "Treatment", "external_aliquot_id"
)
missing_s7_columns <- setdiff(required_s7_columns, names(mouse_manifest))
if (length(missing_s7_columns) > 0L) {
  stop("Missing required Table S7 manifest columns: ",
       paste(missing_s7_columns, collapse = ", "))
}
mouse_manifest <- mouse_manifest[
  !is.na(mouse_manifest$Day) & as.character(mouse_manifest$Day) != "28",
  ,
  drop = FALSE
]
mouse_manifest$match_id <- paste(
  as.character(mouse_manifest$Day),
  mouse_manifest$Treatment,
  mouse_manifest$external_aliquot_id,
  sep = "_"
)

miRNA_treatment_tokens <- c(
  "B7H3 CAR" = "B7H3",
  "B7H3 STOP CAR" = "Ctrl",
  "Untreated" = "Untreated"
)
format_s7_mirna_id <- function(sample_id, treatment) {
  if (!grepl("-Tube[0-9]+$", sample_id)) {
    stop("Unexpected miRNA sample ID in Table S7 manifest: ", sample_id)
  }
  token <- unname(miRNA_treatment_tokens[treatment])
  if (length(token) != 1L || is.na(token)) {
    stop("Unexpected miRNA treatment in Table S7 manifest: ", treatment)
  }
  sub(
    "-Tube",
    paste0("-", token, "-Tube"),
    sample_id,
    fixed = TRUE
  )
}
format_s7_rna_id <- function(sample_id) {
  sub("-(\\d+)(-resub)?$", "-Tube\\1\\2", sample_id, perl = TRUE)
}
one_s7_value <- function(values, column, match_id) {
  values <- unique(as.character(values[!is.na(values)]))
  if (length(values) > 1L) {
    stop("Conflicting ", column, " values for Table S7 match_id ", match_id)
  }
  if (length(values) == 0L || !nzchar(values[[1L]])) {
    return("NA")
  }
  values[[1L]]
}

s7_match_ids <- unique(mouse_manifest$match_id)
table_s7_rows <- lapply(s7_match_ids, function(match_id) {
  sample_rows <- mouse_manifest[
    mouse_manifest$match_id == match_id,
    ,
    drop = FALSE
  ]
  rna_rows <- sample_rows[
    sample_rows$experimental_strategy == "RNA-Seq",
    ,
    drop = FALSE
  ]
  mirna_rows <- sample_rows[
    sample_rows$experimental_strategy == "miRNA-Seq",
    ,
    drop = FALSE
  ]
  if (nrow(rna_rows) == 0L && nrow(mirna_rows) == 0L) {
    stop("No RNA-Seq or miRNA-Seq records found for ", match_id)
  }

  rna_ids <- if (nrow(rna_rows) > 0L) {
    unique(vapply(rna_rows$external_sample_id, format_s7_rna_id, character(1)))
  } else {
    "NA"
  }
  mirna_ids <- if (nrow(mirna_rows) > 0L) {
    unique(mapply(
      format_s7_mirna_id,
      mirna_rows$external_sample_id,
      mirna_rows$Treatment,
      USE.NAMES = FALSE
    ))
  } else {
    "NA"
  }
  if (length(rna_ids) != 1L || length(mirna_ids) != 1L) {
    stop("Multiple sample IDs map to Table S7 match_id ", match_id)
  }

  data.frame(
    match_id = match_id,
    RNA_sample_id = rna_ids[[1L]],
    miRNA_sample_id = mirna_ids[[1L]],
    external_aliquot_id = one_s7_value(
      sample_rows$external_aliquot_id, "external_aliquot_id", match_id
    ),
    organism = one_s7_value(sample_rows$organism, "organism", match_id),
    Day = as.integer(one_s7_value(sample_rows$Day, "Day", match_id)),
    Treatment = one_s7_value(sample_rows$Treatment, "Treatment", match_id),
    RNA_library = if (nrow(rna_rows) > 0L) {
      one_s7_value(rna_rows$RNA_library, "RNA_library", match_id)
    } else {
      "NA"
    },
    `miRNA-Seq` = nrow(mirna_rows) > 0L,
    `RNA-Seq` = nrow(rna_rows) > 0L,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
})
table_s7 <- do.call(rbind, table_s7_rows)
for (column in names(table_s7)) {
  if (is.character(table_s7[[column]])) {
    table_s7[[column]] <- gsub(
      "STOP", "Ctrl", table_s7[[column]], fixed = TRUE
    )
  }
}
table_s7_treatment_order <- match(
  table_s7$Treatment,
  c("B7H3 CAR", "B7H3 Ctrl CAR", "Untreated")
)
table_s7_aliquot_number <- as.integer(
  sub("^rep_", "", table_s7$external_aliquot_id)
)
table_s7 <- table_s7[
  order(table_s7$Day, table_s7_treatment_order, table_s7_aliquot_number),
  ,
  drop = FALSE
]
table_s7_path <- file.path(tables_dir, "TableS7.xlsx")
write_xlsx(list("Mouse metadata" = table_s7), table_s7_path)

message("Wrote ", table_s5_path, " with ", length(table_s5_sheets),
        " worksheets.")
message("Wrote ", table_s6_path, " with ", length(table_s6_sheets),
        " worksheets.")
message("Wrote ", table_s7_path, " with ", nrow(table_s7),
        " unique samples; ", sum(table_s7[["RNA-Seq"]] & table_s7[["miRNA-Seq"]]),
        " have both assays.")

# Table S8: begin with the Day 14 and Day 21 miRNA differential expression
# summary. The Day 28 result columns are excluded from this table.
table_s8_input_path <- file.path(
  root_dir,
  "analyses",
  "murine-mirna-differential-expression",
  "results",
  "mirna-differential-expression-deseq2-b7h3-stop-vs-untreated-by-timepoint.tsv"
)
if (!file.exists(table_s8_input_path)) {
  stop("Missing Table S8 input file: ", table_s8_input_path)
}
table_s8_summary <- read_tsv(table_s8_input_path)
table_s8_day28_columns <- grepl("_Day28$", names(table_s8_summary))
table_s8_summary <- table_s8_summary[
  ,
  !table_s8_day28_columns,
  drop = FALSE
]
names(table_s8_summary) <- gsub(
  "stop", "ctrl", names(table_s8_summary), fixed = TRUE
)
if (any(grepl("_Day28$", names(table_s8_summary)))) {
  stop("Day 28 columns remain in the first Table S8 summary file.")
}

table_s8_rna_input_path <- file.path(
  root_dir,
  "analyses",
  "murine-rna-expression",
  "results",
  "rna-differential-expression-deseq2-b7h3-stop-vs-untreated-by-timepoint.tsv"
)
if (!file.exists(table_s8_rna_input_path)) {
  stop("Missing second Table S8 input file: ", table_s8_rna_input_path)
}
table_s8_rna_summary <- read_tsv(table_s8_rna_input_path)
table_s8_rna_day28_columns <- grepl("_Day28$", names(table_s8_rna_summary))
table_s8_rna_summary <- table_s8_rna_summary[
  ,
  !table_s8_rna_day28_columns,
  drop = FALSE
]
names(table_s8_rna_summary) <- gsub(
  "stop", "ctrl", names(table_s8_rna_summary), fixed = TRUE
)
if (any(grepl("_Day28$", names(table_s8_rna_summary)))) {
  stop("Day 28 columns remain in the second Table S8 summary file.")
}

table_s8_path <- file.path(tables_dir, "TableS8.xlsx")
table_s8_sheets <- list(
  "miRNA DE summary" = table_s8_summary,
  "RNA DE summary" = table_s8_rna_summary
)
write_xlsx(table_s8_sheets, table_s8_path)
message("Wrote ", table_s8_path, " with ", length(table_s8_sheets),
        " summary worksheets (miRNA: ", nrow(table_s8_summary), " rows; RNA: ",
        nrow(table_s8_rna_summary), " rows).")

# Table S9: murine miRNA cluster membership.
table_s9_input_path <- file.path(
  root_dir,
  "analyses",
  "murine-mirna-clustering",
  "results",
  "mouse-de-mirna-cluster-membership.tsv"
)
if (!file.exists(table_s9_input_path)) {
  stop("Missing Table S9 input file: ", table_s9_input_path)
}
table_s9 <- read_tsv(table_s9_input_path)
table_s9 <- table_s9[c("mirna_id", "row_cluster", "expr_pattern")]
table_s9_path <- file.path(tables_dir, "TableS9.xlsx")
write_xlsx(list("Cluster membership" = table_s9), table_s9_path)
message("Wrote ", table_s9_path, " with ", nrow(table_s9), " rows and ",
        ncol(table_s9), " columns.")

# Table S10: GO term enrichment for genes upregulated after B7H3 treatment.
table_s10_files <- c(
  "D14_B7H3_up" = "b7h3_up_d14-enriched-go-terms.tsv",
  "D21_B7H3_up" = "b7h3_up_d21-enriched-go-terms.tsv"
)
table_s10_sheets <- lapply(table_s10_files, function(file_name) {
  file_path <- file.path(
    root_dir, "analyses", "murine-rna-expression", "results", file_name
  )
  if (!file.exists(file_path)) {
    stop("Missing Table S10 input file: ", file_path)
  }
  read_tsv(file_path)
})
expected_s10_columns <- names(table_s10_sheets[[1L]])
if (!all(vapply(table_s10_sheets, function(sheet) {
  identical(names(sheet), expected_s10_columns)
}, logical(1)))) {
  stop("The Table S10 input files do not have the same columns.")
}
table_s10_path <- file.path(tables_dir, "TableS10.xlsx")
write_xlsx(table_s10_sheets, table_s10_path)
message("Wrote ", table_s10_path, " with ", length(table_s10_sheets),
        " worksheets.")
