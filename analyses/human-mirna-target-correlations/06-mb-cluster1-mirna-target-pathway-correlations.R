# Correlate MB cluster 1 miRNA expression with enriched immune-pathway targets
#
# This is the medulloblastoma cluster 1 counterpart to
# 05-dmg-cluster6-mirna-target-pathway-correlations.R. It identifies
# upregulated cluster 1 miRNAs, calculates pathway-level target expression in
# matched MB miRNA-seq/RNA-seq samples, and plots FDR-significant correlations.

library(tidyverse)

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
go_analysis_dir <- file.path(root_dir, "analyses", "human-mirna-go-enrichment")
analysis_dir <- file.path(root_dir, "analyses", "human-mirna-target-correlations")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

source(file.path(root_dir, "figures", "theme.R"))

cluster_file <- file.path(root_dir, "analyses", "mirna-clustering", "results",
                          "MB-de-mirna-cluster-membership-immune-scores.tsv")
target_de_files <- c(
  healthy_normal = file.path(go_analysis_dir, "results",
                             "known-mirna-target-gene-differential-expr-MB-healthyNormal.tsv"),
  paired = file.path(go_analysis_dir, "results",
                     "known-mirna-target-gene-differential-expr-MB-paired.tsv")
)
output_file <- file.path(results_dir,
                         "mb-cluster1-mirna-immune-go-target-differential-expression.tsv")
sample_means_file <- file.path(results_dir,
                               "mb-cluster1-mirna-pathway-target-mean-tpm.tsv")
correlations_file <- file.path(results_dir,
                               "mb-cluster1-mirna-pathway-target-mean-pearson-correlations.tsv")
histologies_file <- file.path(root_dir, "analyses", "histology-preprocessing",
                              "results", "histologies.tsv")
mirna_tpm_file <- file.path(root_dir, "analyses", "immune-deconvolution",
                            "results", "mirna-tpm.rds")
target_tpm_file <- file.path(root_dir, "data",
                             "gene-expression-rsem-tpm-collapsed.all.rds")

# Keep the immune GO biological-process definition aligned with the DMG
# counterpart. MAPK-cascade and cytokinesis terms are excluded downstream.
immune_keywords <- paste(
  "inflamm|immunity|T cell|natural killer cell|JAK|immune|leukocyte|lymphocyte",
  "myeloid|macrophage|transforming growth factor beta|defense|cytokine",
  "chemokine|interferon|glycolysis|lactate|antigen|MHC|HLA|antibody|humoral",
  "immunoglobulin|defensin|granzyme|perforin|necrosis factor|NFAT|interleukin",
  sep = "|"
)
excluded_go_terms <- "MAPK cascade|cytokinesis"

# The direction filter retains the miRNAs whose higher MB expression is
# expected to be inversely related to target expression.
cluster1_mirnas <- read_tsv(cluster_file, show_col_types = FALSE) %>%
  filter(row_cluster == 1, direction == "up") %>%
  pull(miRNA) %>%
  unique()

if (length(cluster1_mirnas) == 0) {
  stop("No upregulated MB cluster 1 miRNAs were found in the cluster table.")
}

go_files <- list.files(file.path(go_analysis_dir, "results"),
                       pattern = "-MB-shared-GO-BP_all\\.tsv$", full.names = TRUE)

immune_go_terms <- tibble(go_file = go_files) %>%
  mutate(miRNA = str_remove(basename(go_file), "-MB-shared-GO-BP_all\\.tsv$")) %>%
  filter(miRNA %in% cluster1_mirnas) %>%
  mutate(go_data = map(go_file, ~ read_tsv(.x, show_col_types = FALSE))) %>%
  dplyr::select(-go_file) %>%
  unnest(go_data) %>%
  filter(str_detect(Term, regex(immune_keywords, ignore_case = TRUE)),
         !str_detect(Term, regex(excluded_go_terms, ignore_case = TRUE)),
         nchar(Term) < 75) %>%
  dplyr::select(miRNA, GO.ID, Term, Annotated, Significant, Expected, classicFisher, fdr, Genes)

if (nrow(immune_go_terms) == 0) {
  stop("No immune-related GO terms were found for upregulated MB cluster 1 miRNAs.")
}

# Genes in GO tables are annotations. Intersect them with the per-miRNA DE
# target lists to retain the targets contributing to each enrichment result.
target_de <- imap_dfr(target_de_files, ~ read_tsv(.x, show_col_types = FALSE) %>%
                        mutate(contrast = .y))

term_target_genes <- immune_go_terms %>%
  separate_rows(Genes, sep = ",\\s*") %>%
  dplyr::rename(target_gene_id = Genes) %>%
  inner_join(target_de, by = c("miRNA", "target_gene_id")) %>%
  distinct()

write_tsv(term_target_genes, output_file)
message("Wrote ", nrow(term_target_genes), " miRNA-target-GO-term rows to ", output_file)

pathway_targets <- term_target_genes %>%
  distinct(miRNA, GO.ID, Term, target_gene_id)

if (nrow(pathway_targets) == 0) {
  stop("No MB miRNA-pathway target relationships remained after DE-target matching.")
}

# Match assays by external sample ID rather than relying on independently
# sorted TPM columns.
histologies <- read_tsv(histologies_file, show_col_types = FALSE)
mirna_tpm <- readRDS(mirna_tpm_file)
target_tpm <- readRDS(target_tpm_file)

mb_mirna_samples <- histologies %>%
  filter(experimental_strategy == "miRNA-Seq", sample_type == "Tumor",
         histology == "MB", external_sample_id %in% colnames(mirna_tpm)) %>%
  dplyr::select(external_sample_id) %>% distinct()
mb_rna_samples <- histologies %>%
  filter(experimental_strategy == "RNA-Seq", sample_type == "Tumor",
         histology == "MB", Bioassay_ID %in% colnames(target_tpm)) %>%
  dplyr::select(external_sample_id, Bioassay_ID) %>% distinct()
mb_samples <- inner_join(mb_mirna_samples, mb_rna_samples, by = "external_sample_id")

if (nrow(mb_samples) < 3) {
  stop("Fewer than three MB tumors have matched miRNA-seq and RNA-seq TPM data.")
}

targets_in_tpm <- pathway_targets %>%
  filter(target_gene_id %in% rownames(target_tpm))
if (nrow(targets_in_tpm) == 0) {
  stop("None of the MB pathway target genes were found in the RNA TPM matrix.")
}

pathway_target_means <- as.data.frame(
  target_tpm[unique(targets_in_tpm$target_gene_id), mb_samples$Bioassay_ID, drop = FALSE]
) %>%
  rownames_to_column("target_gene_id") %>%
  pivot_longer(-target_gene_id, names_to = "Bioassay_ID", values_to = "target_tpm") %>%
  inner_join(mb_samples, by = "Bioassay_ID") %>%
  inner_join(targets_in_tpm, by = "target_gene_id") %>%
  group_by(external_sample_id, miRNA, GO.ID, Term) %>%
  summarise(n_target_genes = n(), mean_target_tpm = mean(target_tpm, na.rm = TRUE),
            .groups = "drop")

mirnas_in_tpm <- intersect(unique(pathway_target_means$miRNA), rownames(mirna_tpm))
if (length(mirnas_in_tpm) == 0) {
  stop("None of the MB pathway miRNAs were found in the miRNA TPM matrix.")
}

mirna_expression <- as.data.frame(mirna_tpm[mirnas_in_tpm,
                                             mb_samples$external_sample_id, drop = FALSE]) %>%
  rownames_to_column("miRNA") %>%
  pivot_longer(-miRNA, names_to = "external_sample_id", values_to = "mirna_tpm")

sample_means <- pathway_target_means %>%
  filter(miRNA %in% mirnas_in_tpm) %>%
  inner_join(mirna_expression, by = c("miRNA", "external_sample_id")) %>%
  dplyr::rename(sample_id = external_sample_id) %>%
  arrange(miRNA, GO.ID, sample_id)

correlations <- sample_means %>%
  group_by(miRNA, GO.ID, Term) %>%
  summarise(
    n_samples = n(), n_target_genes = dplyr::first(n_target_genes),
    pearson_r = if (n() >= 3 && sd(log2(mirna_tpm + 0.01)) > 0 &&
                     sd(log2(mean_target_tpm + 0.01)) > 0)
      cor(log2(mirna_tpm + 0.01), log2(mean_target_tpm + 0.01)) else NA_real_,
    pearson_p = if (n() >= 3 && sd(log2(mirna_tpm + 0.01)) > 0 &&
                     sd(log2(mean_target_tpm + 0.01)) > 0)
      cor.test(log2(mirna_tpm + 0.01), log2(mean_target_tpm + 0.01),
               method = "spearman")$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  group_by(miRNA) %>%
  mutate(pearson_fdr_within_mirna = p.adjust(pearson_p, method = "BH")) %>%
  ungroup() %>%
  arrange(pearson_p, miRNA, Term)

write_tsv(sample_means, sample_means_file)
write_tsv(correlations, correlations_file)

correlation_plot_dir <- file.path(plot_dir, "correlation-plots")
dir.create(correlation_plot_dir, recursive = TRUE, showWarnings = FALSE)

significant_correlations <- correlations %>%
  filter(!is.na(pearson_fdr_within_mirna), pearson_fdr_within_mirna < 0.05)

for (i in seq_len(nrow(significant_correlations))) {
  correlation_row <- significant_correlations[i, ]
  plot_data <- sample_means %>%
    filter(miRNA == correlation_row$miRNA, GO.ID == correlation_row$GO.ID)
  plot_stub <- str_replace_all(paste(correlation_row$miRNA, correlation_row$GO.ID,
                                     sep = "-"), "[^A-Za-z0-9_-]", "_")

  plot <- ggplot(plot_data, aes(x = log2(mirna_tpm + 0.01),
                                y = log2(mean_target_tpm + 0.01))) +
    geom_point(size = 2.5, alpha = 0.85) +
    geom_smooth(method = "lm", se = FALSE, colour = "darkred") +
    labs(title = str_wrap(correlation_row$Term, width = 25),
         subtitle = paste0(correlation_row$miRNA, "\\nPearson r = ",
                           round(correlation_row$pearson_r, 2), "; p = ",
                           signif(correlation_row$pearson_p, 3)),
         x = "log2(miRNA TPM)", y = "log2(mean target TPM)") +
    theme_Publication() +
    theme(plot.title = element_text(size = rel(0.85)),
          plot.subtitle = element_text(size = rel(0.7)),
          axis.title = element_text(size = rel(0.8)))

  ggsave(file.path(correlation_plot_dir,
                   paste0("mb-cluster1-mirna-pathway-correlation-", plot_stub, ".pdf")),
         plot, width = 3, height = 3.5)
}

message("Wrote ", nrow(sample_means), " miRNA-pathway sample means to ", sample_means_file)
message("Wrote ", nrow(correlations), " miRNA-pathway Pearson correlations to ", correlations_file)
sessionInfo()
