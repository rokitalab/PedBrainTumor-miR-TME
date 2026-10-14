# Assess differential expression of immune-related targets of human DMG cluster 6 miRNA orthologs in mouse samples
#
# Ryan Corbett
#
# Oct 2026

# This script performs the following:
# 1. Identifies high-confidence miRanda-predicted targets of murine orthologs of human DMG cluster 6 miRNAs
#    that are DE in B7H3 treated versus untreated mice
# 2. Filters for targets annotated to immune-related GO biological processes (including all descendant terms)
# 3. Tests whether immune targets are more likely to be up- or downregulated in B7H3 treated versus untreated mice
#    at Day 14 and Day 21:
#    - Binomial test: among DE immune targets, is the proportion downregulated different from 0.5?
#    - Fisher's exact tests: are immune targets more likely to be downregulated (or upregulated) than
#      immune genes that are not predicted targets of the miRNA?

### Setup
suppressPackageStartupMessages({
  library(rprojroot)
  library(tidyverse)
  library(ggplot2)
  library(ggsci)
  library(org.Mm.eg.db)
  library(GO.db)
})

### Paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
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

# Immune-related GO biological processes; genes annotated to these terms or any descendant term are considered immune-related
immune_go_ids <- c("GO:0002376", # immune system process
                   "GO:0006954", # inflammatory response
                   "GO:0001816", # cytokine production
                   "GO:0019221", # cytokine-mediated signaling pathway
                   "GO:0006952", # defense response
                   "GO:0034097", # response to cytokine
                   "GO:0006909") # phagocytosis

timepoints <- c("Day14", "Day21")

# Load hsa DMG miRNA cluster df and filter for cluster 6 miRNAs
dmg_clust6_mirs <- read_tsv(dmg_cluster_file) %>%
  dplyr::filter(row_cluster == 6,
                grepl("hsa", miRNA)) %>%
  # remove hsa prefix from miRNA ID
  dplyr::mutate(miRNA = str_remove(miRNA, "hsa-")) %>%
  pull(miRNA)

# Get cluster 6 orthologs DE in B7H3 versus untreated mice
de_clust6_orthos_df <- read_tsv(mmu_mirna_de_file) %>%
  dplyr::filter(str_remove(mirna_id, "mmu-") %in% dmg_clust6_mirs,
                expr_pattern != "Not DE") %>%
  dplyr::select(mirna_id, status_Day14, status_Day21, expr_pattern) %>%
  dplyr::rename(mirna_status_Day14 = status_Day14,
                mirna_status_Day21 = status_Day21,
                mirna_expr_pattern = expr_pattern)

# Load high-confidence miRanda target predictions for DE cluster 6 orthologs
target_df <- read_csv(mirna_target_file) %>%
  dplyr::filter(Seq1 %in% de_clust6_orthos_df$mirna_id,
                Max_Score >= 145, Max_Energy <= -20) %>%
  # retain the highest scoring binding site per miRNA-gene pair
  dplyr::arrange(desc(Max_Score)) %>%
  distinct(Seq1, GeneSymbol, .keep_all = TRUE) %>%
  dplyr::select(mirna_id = Seq1, GeneSymbol, Max_Score, Max_Energy)

# Identify immune-related genes from GO annotations, recording the immune GO categories each gene is annotated to
immune_gene_df <- AnnotationDbi::select(org.Mm.eg.db,
                                        keys = immune_go_ids,
                                        keytype = "GOALL",
                                        columns = "SYMBOL") %>%
  dplyr::filter(ONTOLOGYALL == "BP", !is.na(SYMBOL)) %>%
  distinct(GOALL, SYMBOL) %>%
  dplyr::mutate(immune_category = AnnotationDbi::Term(GOALL)) %>%
  group_by(GeneSymbol = SYMBOL) %>%
  summarise(immune_categories = paste(sort(unique(immune_category)), collapse = "; "),
            .groups = "drop")

# Load mouse DE gene results and reformat to long format by timepoint
# Gene universe per timepoint = genes with a non-NA adjusted p-value in B7H3 vs. untreated comparison
rna_de_long_df <- read_tsv(mmu_rna_de_file) %>%
  dplyr::select(gene_symbol,
                matches("^(log2FoldChange|padj)_b7h3_vs_untr_Day(14|21)$"),
                status_Day14, status_Day21) %>%
  pivot_longer(-gene_symbol,
               names_to = c(".value", "time"),
               names_pattern = "(.*)_(Day\\d+)$") %>%
  dplyr::rename(GeneSymbol = gene_symbol,
                log2FoldChange = log2FoldChange_b7h3_vs_untr,
                padj = padj_b7h3_vs_untr) %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::mutate(direction = case_when(
    status == "B7H3 upregulated" ~ "Up",
    status == "B7H3 downregulated" ~ "Down",
    TRUE ~ "Not DE"
  ))

immune_de_df <- rna_de_long_df %>%
  dplyr::filter(GeneSymbol %in% immune_gene_df$GeneSymbol)

# Write immune target table with DE results per timepoint
immune_target_df <- target_df %>%
  inner_join(immune_gene_df, by = "GeneSymbol") %>%
  inner_join(immune_de_df %>%
               dplyr::select(GeneSymbol, time, log2FoldChange, padj, direction) %>%
               pivot_wider(names_from = time,
                           values_from = c(log2FoldChange, padj, direction)),
             by = "GeneSymbol") %>%
  left_join(de_clust6_orthos_df, by = "mirna_id") %>%
  dplyr::slice(gtools::mixedorder(str_remove(mirna_id, "mmu-miR-"))) %>%
  dplyr::select(mirna_id, mirna_expr_pattern, mirna_status_Day14, mirna_status_Day21,
                GeneSymbol, Max_Score, Max_Energy,
                ends_with("Day14"), ends_with("Day21"),
                immune_categories)

write_tsv(immune_target_df,
          file.path(results_dir,
                    "hsa-dmg-clust6-mouse-ortholog-immune-targets.tsv"))

### Test whether immune targets are preferentially up- or downregulated in B7H3 treated mice

# Define target sets per miRNA, plus union of targets across all cluster 6 orthologs
mirna_order <- de_clust6_orthos_df$mirna_id[gtools::mixedorder(str_remove(de_clust6_orthos_df$mirna_id, "mmu-miR-"))]

target_sets <- c(split(target_df$GeneSymbol, target_df$mirna_id)[mirna_order],
                 list("All cluster 6 orthologs" = unique(target_df$GeneSymbol)))

test_direction <- function(target_genes, tp) {

  tp_df <- immune_de_df %>%
    dplyr::filter(time == tp) %>%
    dplyr::mutate(is_target = GeneSymbol %in% target_genes)

  tgt <- tp_df %>% dplyr::filter(is_target)
  bg <- tp_df %>% dplyr::filter(!is_target)

  n_up <- sum(tgt$direction == "Up")
  n_down <- sum(tgt$direction == "Down")

  # binomial test of direction among DE immune targets
  binom_p <- if (n_up + n_down > 0) binom.test(n_down, n_up + n_down, p = 0.5)$p.value else NA_real_

  # Fisher's exact tests comparing immune targets to immune non-targets
  fisher_down <- fisher.test(table(factor(tp_df$is_target, levels = c(TRUE, FALSE)),
                                   factor(tp_df$direction == "Down", levels = c(TRUE, FALSE))))
  fisher_up <- fisher.test(table(factor(tp_df$is_target, levels = c(TRUE, FALSE)),
                                 factor(tp_df$direction == "Up", levels = c(TRUE, FALSE))))

  tibble(n_immune_targets = nrow(tgt),
         n_up = n_up,
         n_down = n_down,
         pct_up = 100 * n_up / nrow(tgt),
         pct_down = 100 * n_down / nrow(tgt),
         n_immune_nontargets = nrow(bg),
         pct_up_nontargets = 100 * mean(bg$direction == "Up"),
         pct_down_nontargets = 100 * mean(bg$direction == "Down"),
         binom_p = binom_p,
         down_odds_ratio = unname(fisher_down$estimate),
         down_or_lower = fisher_down$conf.int[1],
         down_or_upper = fisher_down$conf.int[2],
         down_fisher_p = fisher_down$p.value,
         up_odds_ratio = unname(fisher_up$estimate),
         up_or_lower = fisher_up$conf.int[1],
         up_or_upper = fisher_up$conf.int[2],
         up_fisher_p = fisher_up$p.value)
}

stat_df <- expand_grid(mirna_id = names(target_sets),
                       time = timepoints) %>%
  dplyr::mutate(res = map2(mirna_id, time,
                           ~test_direction(target_sets[[.x]], .y))) %>%
  unnest(res) %>%
  left_join(de_clust6_orthos_df %>%
              pivot_longer(c(mirna_status_Day14, mirna_status_Day21),
                           names_to = "time",
                           names_prefix = "mirna_status_",
                           values_to = "mirna_status") %>%
              dplyr::select(mirna_id, time, mirna_status),
            by = c("mirna_id", "time")) %>%
  dplyr::relocate(mirna_status, .after = time) %>%
  # adjust p-values within each test across miRNAs and timepoints
  dplyr::mutate(binom_padj = p.adjust(binom_p, method = "BH"),
                down_fisher_padj = p.adjust(down_fisher_p, method = "BH"),
                up_fisher_padj = p.adjust(up_fisher_p, method = "BH"))

print(stat_df %>%
        dplyr::select(mirna_id, time, mirna_status, n_immune_targets, n_up, n_down,
                      pct_up, pct_up_nontargets, pct_down, pct_down_nontargets,
                      binom_padj, down_odds_ratio, down_fisher_padj,
                      up_odds_ratio, up_fisher_padj),
      n = Inf, width = Inf)

write_tsv(stat_df,
          file.path(results_dir,
                    "hsa-dmg-clust6-mouse-ortholog-immune-target-direction-tests.tsv"))

### Plot percentage of immune targets up- and downregulated versus immune non-targets

signif_label <- function(p) {
  case_when(p < 0.001 ~ "***",
            p < 0.01 ~ "**",
            p < 0.05 ~ "*",
            p < 0.1 ~ "#",
            TRUE ~ "")
}

plot_df <- stat_df %>%
  dplyr::select(mirna_id, time,
                Targets_Up = pct_up, Targets_Down = pct_down,
                `Non-targets_Up` = pct_up_nontargets, `Non-targets_Down` = pct_down_nontargets) %>%
  pivot_longer(-c(mirna_id, time),
               names_to = c("gene_set", "direction"),
               names_sep = "_",
               values_to = "pct") %>%
  dplyr::mutate(mirna_id = factor(mirna_id, levels = rev(names(target_sets))),
                gene_set = factor(gene_set, levels = c("Targets", "Non-targets")),
                direction = factor(direction, levels = c("Up", "Down")))

label_df <- stat_df %>%
  dplyr::select(mirna_id, time, Up = up_fisher_padj, Down = down_fisher_padj) %>%
  pivot_longer(c(Up, Down),
               names_to = "direction",
               values_to = "padj") %>%
  left_join(plot_df %>%
              group_by(mirna_id, time, direction) %>%
              summarise(pct = max(pct), .groups = "drop") %>%
              dplyr::mutate(mirna_id = as.character(mirna_id),
                            direction = as.character(direction)),
            by = c("mirna_id", "time", "direction")) %>%
  dplyr::mutate(label = signif_label(padj),
                mirna_id = factor(mirna_id, levels = rev(names(target_sets))),
                direction = factor(direction, levels = c("Up", "Down")))

immune_target_de_plot <- plot_df %>%
  ggplot(aes(x = pct, y = mirna_id, fill = gene_set)) +
  geom_col(position = position_dodge(width = 0.8),
           width = 0.75) +
  geom_text(data = label_df,
            aes(x = pct, y = mirna_id, label = label),
            inherit.aes = FALSE,
            hjust = -0.3, vjust = 0.75,
            size = 5) +
  scale_fill_manual(values = c("Targets" = pal_npg()(4)[4],
                               "Non-targets" = "grey70")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.2))) +
  facet_grid(time ~ direction, labeller = labeller(direction = c(Up = "Upregulated in B7H3",
                                                                 Down = "Downregulated in B7H3"))) +
  labs(x = "% of immune-related genes",
       y = NULL,
       fill = NULL,
       caption = "Fisher's exact test (immune targets vs. immune non-targets), BH-adjusted: # p < 0.1, * p < 0.05, ** p < 0.01, *** p < 0.001") +
  theme_Publication() +
  theme(plot.caption = element_text(size = 8))

ggsave(file.path(plot_dir,
                 "hsa-dmg-clust6-mouse-ortholog-immune-target-de-barplot.pdf"),
       plot = immune_target_de_plot,
       width = 9, height = 6)

# print session info
sessionInfo()
