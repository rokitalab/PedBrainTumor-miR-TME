# Correlate DMG cluster 6 miRNA expression with enriched immune-pathway targets
#
# Ryan Corbett
#
# Sep 2026

# This script identifies upregulated DMG cluster 6 miRNAs, finds their enriched
# immune-related GO biological-process terms and DE target genes, calculates
# miRNA/pathway-target expression correlations, and plots significant results.

library(tidyverse)

# Set directory paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

go_analysis_dir <- file.path(root_dir, "analyses", "human-mirna-go-enrichment")
analysis_dir <- file.path(root_dir, "analyses", "human-mirna-target-correlations")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

source(file.path(root_dir, "figures", "theme.R"))

# Set input file paths
dmg_cluster_file <- file.path(root_dir, "analyses", "mirna-clustering", "results",
                              "DIPG or DMG-de-mirna-cluster-membership-immune-scores.tsv")

target_de_files <- c(
  healthy_normal = file.path(go_analysis_dir, "results",
                             "known-mirna-target-gene-differential-expr-DIPG or DMG-healthyNormal.tsv"),
  paired = file.path(go_analysis_dir, "results",
                     "known-mirna-target-gene-differential-expr-DIPG or DMG-paired.tsv")
)

output_file <- file.path(results_dir,
                         "dmg-cluster6-mirna-immune-go-target-differential-expression.tsv")
histologies_file <- file.path(root_dir, "analyses", "histology-preprocessing",
                              "results", "histologies.tsv")
mirna_tpm_file <- file.path(root_dir, "analyses", "immune-deconvolution",
                            "results", "mirna-tpm.rds")
target_tpm_file <- file.path(root_dir, "data",
                             "gene-expression-rsem-tpm-collapsed.all.rds")
sample_means_file <- file.path(
  results_dir,
  "dmg-cluster6-mirna-pathway-target-mean-tpm.tsv"
)
correlations_file <- file.path(
  results_dir,
  "dmg-cluster6-mirna-pathway-target-mean-pearson-correlations.tsv"
)

# Use the same immune-term definition as the cluster 6 GO-enrichment script.
# Exclude MAPK-cascade and cytokinesis terms from the downstream pathway set.
immune_keywords <- paste(
  "inflamm|immunity|T cell|natural killer cell|JAK|immune|leukocyte|lymphocyte",
  "myeloid|macrophage|transforming growth factor beta|defense|cytokine",
  "chemokine|interferon|glycolysis|lactate|antigen|MHC|HLA|antibody|humoral",
  "immunoglobulin|defensin|granzyme|perforin|necrosis factor|NFAT|interleukin",
  sep = "|"
)
excluded_go_terms <- "MAPK cascade|cytokinesis"

# Identify upregulated DMG cluster 6 miRNAs.
cluster6_mirnas <- read_tsv(dmg_cluster_file, show_col_types = FALSE) %>%
  filter(row_cluster == 6, direction == "up") %>%
  pull(miRNA) %>%
  unique()

# Load the per-miRNA GO tables for those miRNAs and retain immune-related terms.
go_files <- list.files(
  file.path(go_analysis_dir, "results"),
  pattern = "-DIPG or DMG-shared-GO-BP_all\\.tsv$",
  full.names = TRUE
)

immune_go_terms <- tibble(go_file = go_files) %>%
  mutate(miRNA = str_remove(basename(go_file), "-DIPG or DMG-shared-GO-BP_all\\.tsv$")) %>%
  filter(miRNA %in% cluster6_mirnas) %>%
  mutate(go_data = map(go_file, ~ read_tsv(.x, show_col_types = FALSE))) %>%
  dplyr::select(-go_file) %>%
  unnest(go_data) %>%
  filter(str_detect(Term, immune_keywords),
         !str_detect(Term, regex(excluded_go_terms, ignore_case = TRUE)),
         nchar(Term) < 75) %>%
  dplyr::select(miRNA, GO.ID, Term, Annotated, Significant, Expected, classicFisher, fdr, Genes)

# `Genes` contains all genes annotated to a GO term.  Intersect it with the
# miRNA-specific DE target lists to obtain the actual target genes contributing
# to each enriched term.  Keep the two DE contrasts separate in the output.
target_de <- imap_dfr(
  target_de_files,
  ~ read_tsv(.x, show_col_types = FALSE) %>% mutate(contrast = .y)
)

term_target_genes <- immune_go_terms %>%
  separate_rows(Genes, sep = ",\\s*") %>%
  dplyr::rename(target_gene_id = Genes) %>%
  inner_join(target_de, by = c("miRNA", "target_gene_id")) %>%
  distinct()

write_tsv(term_target_genes, output_file)

message("Wrote ", nrow(term_target_genes), " miRNA-target-GO-term rows to ", output_file)

# Retain each miRNA-GO term-target relationship once; the target table has one
# row per DE contrast, whereas this expression analysis is contrast-independent.
pathway_targets <- term_target_genes %>%
  dplyr::filter(!str_detect(Term, regex(excluded_go_terms, ignore_case = TRUE))) %>%
  distinct(miRNA, GO.ID, Term, target_gene_id)

# Match DMG tumor miRNA-seq and RNA-seq samples by their external sample ID.
# This prevents relying on separate sort orders for the two TPM matrices.
histologies <- read_tsv(histologies_file, show_col_types = FALSE)
mirna_tpm <- readRDS(mirna_tpm_file)
target_tpm <- readRDS(target_tpm_file)

dmg_mirna_samples <- histologies %>%
  dplyr::filter(experimental_strategy == "miRNA-Seq",
                sample_type == "Tumor",
                histology == "DIPG or DMG",
                external_sample_id %in% colnames(mirna_tpm)) %>%
  dplyr::select(external_sample_id) %>%
  distinct()

dmg_rna_samples <- histologies %>%
  dplyr::filter(experimental_strategy == "RNA-Seq",
                sample_type == "Tumor",
                histology == "DIPG or DMG",
                Bioassay_ID %in% colnames(target_tpm)) %>%
  dplyr::select(external_sample_id, Bioassay_ID) %>%
  distinct()

dmg_samples <- inner_join(dmg_mirna_samples, dmg_rna_samples,
                          by = "external_sample_id")

if (nrow(dmg_samples) < 3) {
  stop("Fewer than three DMG tumors have matched miRNA-seq and RNA-seq TPM data.")
}

# Calculate mean TPM of each miRNA-pathway target set in every DMG tumor.
targets_in_tpm <- pathway_targets %>%
  dplyr::filter(target_gene_id %in% rownames(target_tpm))

if (nrow(targets_in_tpm) == 0) {
  stop("None of the pathway target genes were found in the RNA TPM matrix.")
}

pathway_target_means <- as.data.frame(
  target_tpm[unique(targets_in_tpm$target_gene_id), dmg_samples$Bioassay_ID, drop = FALSE]
) %>%
  rownames_to_column("target_gene_id") %>%
  pivot_longer(-target_gene_id, names_to = "Bioassay_ID", values_to = "target_tpm") %>%
  inner_join(dmg_samples, by = "Bioassay_ID") %>%
  inner_join(targets_in_tpm, by = "target_gene_id") %>%
  group_by(external_sample_id, miRNA, GO.ID, Term) %>%
  summarise(n_target_genes = n(),
            mean_target_tpm = mean(target_tpm, na.rm = TRUE),
            .groups = "drop")

# Add the matching miRNA TPM for each DMG tumor and miRNA-pathway combination.
mirnas_in_tpm <- intersect(unique(pathway_target_means$miRNA), rownames(mirna_tpm))

if (length(mirnas_in_tpm) == 0) {
  stop("None of the pathway miRNAs were found in the miRNA TPM matrix.")
}

pathway_target_means <- pathway_target_means %>%
  dplyr::filter(miRNA %in% mirnas_in_tpm)

mirna_expression <- as.data.frame(mirna_tpm[mirnas_in_tpm,
                                            dmg_samples$external_sample_id,
                                            drop = FALSE]) %>%
  rownames_to_column("miRNA") %>%
  pivot_longer(-miRNA, names_to = "external_sample_id", values_to = "mirna_tpm")

sample_means <- pathway_target_means %>%
  inner_join(mirna_expression, by = c("miRNA", "external_sample_id")) %>%
  dplyr::rename(sample_id = external_sample_id) %>%
  arrange(miRNA, GO.ID, sample_id)

# Compute Pearson correlations across matched DMG tumors. Both the correlation
# and p-value use log2(TPM + 0.01), which stabilizes the TPM scale while
# retaining zero-expression observations.
correlations <- sample_means %>%
  group_by(miRNA, GO.ID, Term) %>%
  summarise(
    n_samples = n(),
    n_target_genes = dplyr::first(n_target_genes),
    pearson_r = if (n() >= 3 &&
                     sd(log2(mirna_tpm + 0.01)) > 0 &&
                     sd(log2(mean_target_tpm + 0.01)) > 0) {
      cor(log2(mirna_tpm + 0.01), log2(mean_target_tpm + 0.01),
          method = "pearson")
    } else {
      NA_real_
    },
    pearson_p = if (n() >= 3 &&
                     sd(log2(mirna_tpm + 0.01)) > 0 &&
                     sd(log2(mean_target_tpm + 0.01)) > 0) {
      cor.test(log2(mirna_tpm + 0.01), log2(mean_target_tpm + 0.01),
               method = "pearson")$p.value
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  group_by(miRNA) %>%
  mutate(pearson_fdr_within_mirna = p.adjust(pearson_p, method = "BH")) %>%
  ungroup() %>%
  arrange(pearson_p, miRNA, Term)

write_tsv(sample_means, sample_means_file)
write_tsv(correlations, correlations_file)

message("Wrote ", nrow(sample_means), " miRNA-pathway sample means to ",
        sample_means_file)
message("Wrote ", nrow(correlations), " miRNA-pathway Pearson correlations to ",
        correlations_file)

# Plot every FDR-significant correlation. Axes show log2 TPM values; Pearson
# statistics remain based on log2(TPM + 0.01).
correlation_plot_dir <- file.path(plot_dir, "correlation-plots")
dir.create(correlation_plot_dir, recursive = TRUE, showWarnings = FALSE)

significant_correlations <- correlations %>%
  dplyr::filter(!is.na(pearson_fdr_within_mirna),
                pearson_fdr_within_mirna < 0.05)

for (i in seq_len(nrow(significant_correlations))) {
  correlation_row <- significant_correlations[i, ]
  plot_data <- sample_means %>%
    dplyr::filter(miRNA == correlation_row$miRNA,
                  GO.ID == correlation_row$GO.ID)

  plot_file_stub <- paste(correlation_row$miRNA, correlation_row$GO.ID,
                          sep = "-") %>%
    str_replace_all("[^A-Za-z0-9_-]", "_")

  correlation_plot <- ggplot(plot_data,
                             aes(x = log2(mirna_tpm + 0.01),
                                 y = log2(mean_target_tpm + 0.01))) +
    geom_point(size = 2.5, alpha = 0.85) +
    geom_smooth(method = "lm", se = FALSE, colour = "darkred") +
    labs(
      title = str_wrap(correlation_row$Term, width = 25),
      subtitle = paste0(correlation_row$miRNA,
                        "\nPearson r = ", round(correlation_row$pearson_r, 2),
                        "; p = ", signif(correlation_row$pearson_p, 3)),
      x = "log2(miRNA TPM)",
      y = "log2(mean target TPM)"
    ) +
    theme_Publication() +
    theme(
      plot.title = element_text(size = rel(0.85)),
      plot.subtitle = element_text(size = rel(0.7)),
      axis.title.x = element_text(size = rel(0.8)),
      axis.title.y = element_text(size = rel(0.8))
    )

  ggsave(file.path(correlation_plot_dir,
                    paste0("dmg-cluster6-mirna-pathway-correlation-",
                           plot_file_stub, ".pdf")),
         correlation_plot, width = 3, height = 3.5)
}

sessionInfo()
