# Plot DE miRNA expression across time points
# Author: Bicna Song | Jan 2026

### Setup
suppressPackageStartupMessages({
  library(tidyverse)
  library(rprojroot)
  library(glue)
  library(rstatix)
})

### Paths
root_dir <- find_root(has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "murine-mirna-clustering")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

source(file.path(root_dir, "figures", "theme.R"))

mirna_tpm_file <- file.path(analysis_dir, "results", "murine-mirna-tpm.rds")
de_results_file <- file.path(analysis_dir, "results", "mouse_sig_DE_miRNA_list.csv")
metadata_file <- file.path(root_dir, "analyses", "murine-mirna-differential-expression",
                           "results", "sample-metadata.tsv")

### Load & preprocess
mirna_tpm <- readRDS(mirna_tpm_file)
mirna_tpm <- log2(mirna_tpm + 1)

de_mirnas <- read_csv(de_results_file) %>%
  pull(mirna_id)

mirna_tpm <- mirna_tpm[rownames(mirna_tpm) %in% de_mirnas, , drop = FALSE]

metadata_df <- read_tsv(metadata_file) %>%
  filter(time %in% c("Day14", "Day21"))

### Convert to long format
tpm_long <- mirna_tpm %>%
  as.data.frame() %>%
  rownames_to_column("mirna_id") %>%
  pivot_longer(
    cols = -mirna_id,
    names_to = "id",
    values_to = "log2_TPM"
  ) %>%
  left_join(metadata_df, by = "id")

### Summary statistics
tpm_summary <- tpm_long %>%
  group_by(mirna_id, treatment, time) %>%
  summarise(
    mean_log2_tpm = mean(log2_TPM, na.rm = TRUE),
    sem_log2_tpm  = sd(log2_TPM, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    time = factor(time, levels = c("Day14", "Day21"))
  )

### Extract significance from DE results
sig_df <- read_csv(de_results_file) %>%
  transmute(
    mirna_id,
    expr_pattern,
    # decide which adjusted p-value to use
    padj = case_when(
      expr_pattern %in% c("B7H3 up, Day14", "B7H3 down, Day14") ~
        padj_b7h3_vs_untr_Day14,
      expr_pattern %in% c("B7H3 up, Day21", "B7H3 down, Day21") ~
        padj_b7h3_vs_untr_Day21,
      expr_pattern %in% c("B7H3 up, multiple", "B7H3 down, multiple") ~
        pmin(
          padj_b7h3_vs_untr_Day14,
          padj_b7h3_vs_untr_Day21,
          na.rm = TRUE
        ),
      TRUE ~ NA_real_
    ),
    # assign time point for plotting
    time = case_when(
      str_detect(expr_pattern, "Day14") ~ "Day14",
      str_detect(expr_pattern, "Day21") ~ "Day21",
      str_detect(expr_pattern, "multiple") ~ "Day21",  # place above Day21
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(padj), !is.na(time)) %>%
  mutate(
    p.signif = case_when(
      padj <= 0.001 ~ "***",
      padj <= 0.01  ~ "**",
      padj <= 0.1   ~ "*",
      TRUE          ~ NA_character_
    )
  )

sig_df <- sig_df %>%
  left_join(
    tpm_summary %>%
      group_by(mirna_id, time) %>%
      summarise(
        y_max = max(mean_log2_tpm + sem_log2_tpm, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("mirna_id", "time")
  ) %>%
  mutate(
    y.position = y_max + 0.25
  )

### Plot function
plot_mirna_timecourse <- function(df, mirna, sig_df) {
  ggplot(df, aes(
    x = time,
    y = mean_log2_tpm,
    color = treatment,
    group = treatment
  )) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(
        ymin = mean_log2_tpm - sem_log2_tpm,
        ymax = mean_log2_tpm + sem_log2_tpm
      ),
      width = 0.15,
      linewidth = 0.6
    ) +
    geom_line(
      aes(linetype = treatment),
      linewidth = 1,
      show.legend = FALSE
    ) +
    geom_text(
      data = sig_df %>% filter(mirna_id == mirna),
      aes(
        x = time,
        y = y.position,
        label = p.signif
      ),
      inherit.aes = FALSE,
      color = "black",
      size = 4
    ) +
    scale_linetype_manual(
      values = c(
        "B7H3" = "solid",
        "STOP" = "dashed",
        "untreated" = "longdash"
      )
    ) +
    scale_color_manual(
      values = c(
        "B7H3" = "#E64B35FF",
        "STOP" = "#4DBBD5FF",
        "untreated" = "#00A087FF"
      )
    ) +
    guides(linetype = "none") +
    labs(
      title = mirna,
      x = NULL,
      y = "log2-TPM",
      color = "Treatment"
    ) +
    theme_Publication()
}

### Generate one plot per miRNA
for (mirna in unique(tpm_summary$mirna_id)) {
  p <- plot_mirna_timecourse(
    df = tpm_summary %>% filter(mirna_id == mirna),
    mirna = mirna,
    sig_df = sig_df
  )
  ggsave(
    file.path(plot_dir, glue("{mirna}_log2TPM_Day14_Day21.pdf")),
    p,
    width = 6,
    height = 4
  )
}

### Session info
sessionInfo()
