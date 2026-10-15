# Plot expression of predicted targets of human DMG cluster 6 miRNA orthologs in mouse samples
#
# Ryan Corbett
#
# Oct 2026

# This script performs the following:
# 1. Identifies murine orthologs of human DMG cluster 6 miRNAs that are DE in B7H3 treated versus untreated mice
# 2. Identifies high-confidence miRanda-predicted targets of these miRNAs
# 3. Plots target gene expression in B7H3 treated and untreated mice at Day 14 and Day 21, for all targets
#    and for targets annotated to immune-related GO biological processes

### Setup
suppressPackageStartupMessages({
  library(rprojroot)
  library(tidyverse)
  library(ggplot2)
  library(ggsci)
  library(org.Mm.eg.db)
})

### Paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "human-mirna-expr-in-mouse")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

source(file.path(root_dir, "figures", "theme.R"))

## set file paths
dmg_cluster_file <- file.path(root_dir, "analyses",
                              "mirna-clustering",
                              "results",
                              "DIPG or DMG-de-mirna-cluster-membership-immune-scores.tsv")

mmu_mirna_de_file <- file.path(root_dir, "analyses",
                               "murine-mirna-differential-expression",
                               "results",
                               "mirna-differential-expression-deseq2-b7h3-stop-vs-untreated-by-timepoint.tsv")

mmu_rna_de_file <- file.path(root_dir, "analyses",
                             "murine-rna-expression",
                             "results",
                             "rna-differential-expression-deseq2-b7h3-stop-vs-untreated-by-timepoint.tsv")

mirna_target_file <- file.path(root_dir, "analyses",
                               "murine-novel-mirna-target-prediction",
                               "results",
                               "mouse_miranda_output_parsed_anno.csv")

rna_tpm_file <- file.path(data_dir,
                          "mouse-gene-expression-rsem-tpm-collapsed.rds")

manifest_file <- file.path(root_dir, "analyses",
                           "murine-rna-expression",
                           "input",
                           "Haydar_Mouse_RNA_miRNA_manifest_IDs_assigned.tsv")

# Load hsa DMG miRNA cluster df and filter for cluster 6 miRNAs
dmg_clust6_mirs <- read_tsv(dmg_cluster_file) %>%
  dplyr::filter(row_cluster == 6,
                grepl("hsa", miRNA)) %>%
  # remove hsa prefix from miRNA ID
  dplyr::mutate(miRNA = str_remove(miRNA, "hsa-")) %>%
  pull(miRNA)

# Get cluster 6 orthologs DE in B7H3 versus untreated mice, retaining Day 14 & 21 DE status
de_clust6_orthos_df <- read_tsv(mmu_mirna_de_file) %>%
  dplyr::filter(str_remove(mirna_id, "mmu-") %in% dmg_clust6_mirs,
                expr_pattern != "Not DE") %>%
  dplyr::select(mirna_id, status_Day14, status_Day21, expr_pattern) %>%
  dplyr::rename(mirna_status_Day14 = status_Day14,
                mirna_status_Day21 = status_Day21,
                mirna_expr_pattern = expr_pattern)

# Load mouse DE gene results, and retain B7H3 vs. untreated results at Day 14 & 21
mmu_rna_de_df <- read_tsv(mmu_rna_de_file) %>%
  dplyr::select(gene_symbol,
                matches("b7h3_vs_untr_Day(14|21)"),
                status_Day14, status_Day21) %>%
  dplyr::rename(target_status_Day14 = status_Day14,
                target_status_Day21 = status_Day21)

# Load miRNA target predictions, filter for high confidence targets of DE cluster6 orthologs,
# and append miRNA and target DE results
target_df <- read_csv(mirna_target_file) %>%
  dplyr::filter(Seq1 %in% de_clust6_orthos_df$mirna_id,
                Max_Score >= 145, Max_Energy <= -20) %>%
  # retain the highest scoring binding site per miRNA-gene pair
  dplyr::arrange(desc(Max_Score)) %>%
  distinct(Seq1, GeneSymbol, .keep_all = TRUE) %>%
  dplyr::rename(mirna_id = Seq1) %>%
  dplyr::select(mirna_id, GeneSymbol, GeneID, Max_Score, Max_Energy) %>%
  left_join(de_clust6_orthos_df, by = "mirna_id") %>%
  left_join(mmu_rna_de_df, by = c("GeneSymbol" = "gene_symbol")) %>%
  dplyr::arrange(mirna_id, GeneSymbol)

write_tsv(target_df,
          file.path(results_dir,
                    "hsa-dmg-clust6-mouse-ortholog-predicted-targets.tsv"))

# Load manifest to define RNA-seq sample treatment and timepoint, excluding STOP CAR and Day 28 samples
rna_meta_df <- read_tsv(manifest_file) %>%
  dplyr::filter(experimental_strategy == "RNA-Seq",
                Day %in% c(14, 21),
                Treatment %in% c("B7H3 CAR", "Untreated")) %>%
  distinct(Bioassay_ID, external_sample_id, Day, Treatment) %>%
  dplyr::mutate(time = glue::glue("Day{Day}"),
                treatment = case_when(
                  Treatment == "B7H3 CAR" ~ "CAR",
                  TRUE ~ Treatment
                ))

# Load RNA TPM matrix and subset for retained samples and predicted target genes
rna_tpm_mat <- readRDS(rna_tpm_file)
names(rna_tpm_mat) <- str_remove(names(rna_tpm_mat), "_[^_]+$")

rna_tpm_mat <- rna_tpm_mat[, rna_meta_df$Bioassay_ID]

# Calculate mean log2(TPM + 1) per gene, timepoint, and treatment for all genes
gene_tpm_df <- rna_tpm_mat %>%
  rownames_to_column("GeneSymbol") %>%
  pivot_longer(-GeneSymbol,
               names_to = "Bioassay_ID",
               values_to = "tpm") %>%
  dplyr::mutate(log2_tpm = log2(tpm + 1)) %>%
  left_join(rna_meta_df, by = "Bioassay_ID") %>%
  group_by(GeneSymbol, time, treatment) %>%
  summarise(mean_log2_tpm = mean(log2_tpm),
            .groups = "drop")

target_tpm_df <- gene_tpm_df %>%
  dplyr::filter(GeneSymbol %in% target_df$GeneSymbol)

# Use same treatment colors as cluster 6 ortholog plot (CAR, Ctrl CAR, Untreated)
trt_colors <- setNames(pal_npg()(3)[c(1, 3)],
                       c("CAR", "Untreated"))

# Identify immune-related genes annotated to immune GO biological processes or any descendant term
# (same GO terms as 04-immune-target-expression.R)
immune_go_ids <- c("GO:0002376", # immune system process
                   "GO:0006954", # inflammatory response
                   "GO:0001816", # cytokine production
                   "GO:0019221", # cytokine-mediated signaling pathway
                   "GO:0006952", # defense response
                   "GO:0034097", # response to cytokine
                   "GO:0006909") # phagocytosis

immune_genes <- AnnotationDbi::select(org.Mm.eg.db,
                                      keys = immune_go_ids,
                                      keytype = "GOALL",
                                      columns = "SYMBOL") %>%
  dplyr::filter(ONTOLOGYALL == "BP", !is.na(SYMBOL)) %>%
  pull(SYMBOL) %>%
  unique()

# Function to run paired Wilcoxon signed-rank tests of target expression in B7H3 versus untreated mice
# within each miRNA and timepoint (paired by gene), and plot target expression as violin plots faceted by miRNA
plot_target_expr <- function(mirna_target_pairs, stat_file, plot_file, y_label) {

  # Merge target expression with miRNA-target pairs, and define facet labels
  mirna_target_tpm_df <- mirna_target_pairs %>%
    distinct(mirna_id, GeneSymbol) %>%
    inner_join(target_tpm_df, by = "GeneSymbol", relationship = "many-to-many") %>%
    group_by(mirna_id) %>%
    dplyr::mutate(n_genes = n_distinct(GeneSymbol)) %>%
    ungroup() %>%
    left_join(de_clust6_orthos_df %>% dplyr::select(mirna_id, mirna_expr_pattern),
              by = "mirna_id") %>%
    dplyr::mutate(facet_label = glue::glue("{mirna_id}\n({mirna_expr_pattern}; n = {n_genes})"))

  # order facets numerically by miRNA (strip prefix so mixedorder does not parse hyphens as negative signs)
  facet_levels <- mirna_target_tpm_df %>%
    distinct(mirna_id, facet_label) %>%
    dplyr::slice(gtools::mixedorder(str_remove(mirna_id, "mmu-miR-"))) %>%
    pull(facet_label)

  mirna_target_tpm_df <- mirna_target_tpm_df %>%
    dplyr::mutate(facet_label = factor(facet_label, levels = facet_levels))

  stat_df <- mirna_target_tpm_df %>%
    pivot_wider(names_from = treatment,
                values_from = mean_log2_tpm) %>%
    group_by(mirna_id, facet_label, time) %>%
    summarise(n_genes = n(),
              median_CAR = median(CAR),
              median_Untreated = median(Untreated),
              median_diff_CAR_vs_Untreated = median(CAR - Untreated),
              p_value = wilcox.test(CAR, Untreated, paired = TRUE, exact = FALSE)$p.value,
              .groups = "drop") %>%
    dplyr::mutate(p_adj = p.adjust(p_value, method = "BH"),
                  p_signif = case_when(
                    p_adj < 0.001 ~ "***",
                    p_adj < 0.01 ~ "**",
                    p_adj < 0.05 ~ "*",
                    TRUE ~ ""
                  )) %>%
    dplyr::arrange(facet_label, time)

  write_tsv(stat_df %>% dplyr::select(-facet_label),
            file.path(results_dir, stat_file))

  # position asterisks above violins within each facet
  stat_label_df <- stat_df %>%
    left_join(mirna_target_tpm_df %>%
                group_by(facet_label) %>%
                summarise(y_pos = max(mean_log2_tpm) * 1.08,
                          .groups = "drop"),
              by = "facet_label")

  target_expr_plot <- mirna_target_tpm_df %>%
    ggplot(aes(x = time,
               y = mean_log2_tpm,
               fill = treatment)) +
    geom_violin(alpha = 0.6,
                scale = "width",
                position = position_dodge(width = 0.8)) +
    geom_boxplot(aes(group = interaction(time, treatment)),
                 width = 0.15,
                 outlier.shape = NA,
                 fill = "white",
                 position = position_dodge(width = 0.8),
                 show.legend = FALSE) +
    geom_text(data = stat_label_df,
              aes(x = time, y = y_pos, label = p_signif),
              inherit.aes = FALSE,
              size = 6) +
    scale_fill_manual(values = trt_colors) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.12))) +
    facet_wrap(~facet_label, ncol = 3) +
    labs(x = NULL,
         y = y_label,
         fill = "Treatment",
         caption = "Paired Wilcoxon signed-rank test (CAR vs. Untreated), BH-adjusted: * p < 0.05, ** p < 0.01, *** p < 0.001") +
    theme_Publication() +
    theme(strip.text = element_text(size = 9),
          plot.caption = element_text(size = 8))

  ggsave(file.path(plot_dir, plot_file),
         plot = target_expr_plot,
         width = 10, height = 7)

  return(stat_df)
}

### Plot expression of all predicted targets

all_target_stat_df <- plot_target_expr(target_df,
                                       stat_file = "hsa-dmg-clust6-mouse-ortholog-target-expr-wilcoxon.tsv",
                                       plot_file = "hsa-dmg-clust6-mouse-ortholog-predicted-target-tpm-violin.pdf",
                                       y_label = "Mean target log2(TPM + 1)")

### Plot expression of immune-related predicted targets

immune_target_stat_df <- plot_target_expr(target_df %>%
                                            dplyr::filter(GeneSymbol %in% immune_genes),
                                          stat_file = "hsa-dmg-clust6-mouse-ortholog-immune-target-expr-wilcoxon.tsv",
                                          plot_file = "hsa-dmg-clust6-mouse-ortholog-immune-target-tpm-violin.pdf",
                                          y_label = "Mean immune target log2(TPM + 1)")

print(immune_target_stat_df %>% dplyr::select(-facet_label), n = Inf, width = Inf)

### Wilcoxon rank-sum tests comparing B7H3 vs. untreated expression differences of predicted targets
### to those of all other expressed genes (non-targets), within each miRNA and timepoint

# Calculate per-gene CAR - Untreated difference in mean log2(TPM + 1) for all expressed genes
gene_diff_df <- gene_tpm_df %>%
  pivot_wider(names_from = treatment,
              values_from = mean_log2_tpm) %>%
  dplyr::filter(CAR + Untreated > 0) %>%
  dplyr::mutate(diff_CAR_vs_Untreated = CAR - Untreated)

bg_stat_df <- expand_grid(mirna_id = unique(target_df$mirna_id),
                          time = unique(gene_diff_df$time)) %>%
  dplyr::mutate(res = map2(mirna_id, time, function(mirna, tp) {

    mirna_targets <- target_df$GeneSymbol[target_df$mirna_id == mirna]

    tp_diff_df <- gene_diff_df %>%
      dplyr::filter(time == tp) %>%
      dplyr::mutate(is_target = GeneSymbol %in% mirna_targets)

    target_diffs <- tp_diff_df$diff_CAR_vs_Untreated[tp_diff_df$is_target]
    nontarget_diffs <- tp_diff_df$diff_CAR_vs_Untreated[!tp_diff_df$is_target]

    test_res <- wilcox.test(target_diffs, nontarget_diffs,
                            conf.int = TRUE, exact = FALSE)

    tibble(n_targets = length(target_diffs),
           n_nontargets = length(nontarget_diffs),
           median_diff_targets = median(target_diffs),
           median_diff_nontargets = median(nontarget_diffs),
           # Hodges-Lehmann estimate of the shift in target vs. non-target differences
           location_shift = unname(test_res$estimate),
           p_value = test_res$p.value)

  })) %>%
  unnest(res) %>%
  dplyr::mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  dplyr::slice(gtools::mixedorder(str_remove(mirna_id, "mmu-miR-")))

print(bg_stat_df, width = Inf)

write_tsv(bg_stat_df,
          file.path(results_dir,
                    "hsa-dmg-clust6-mouse-ortholog-target-vs-nontarget-wilcoxon.tsv"))


# print session info
sessionInfo()
