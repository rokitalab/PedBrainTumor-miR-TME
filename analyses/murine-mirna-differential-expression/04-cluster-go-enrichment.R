# murine-mirna-target GO enrichment for immune-related terms
# Generate immune GO term dot plots for clusters 1, 2, and 5
# Author: Bicna Song | Nov 2025 

### Setup
suppressPackageStartupMessages({
  library(rprojroot)
  library(tidyverse)
  library(ggplot2)
  library(gtools)
  library(forcats)
})

### Paths
root_dir <- find_root(has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "murine-mirna-differential-expression")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

source(file.path(root_dir, "figures", "theme.R"))

cluster_file <- file.path(root_dir, "analyses", "murine-mirna-clustering", "results", "mouse-de-mirna-cluster-membership.tsv")

### Parameters
keywords <- "inflamm|immunity|T cell|natural killer cell|JAK|immune|leukocyte|lymphocyte|myeloid|macrophage|transforming growth factor beta|MAPK cascade|defense|cytokine|chemokine|interferon|glycolysis|lactate|antigen|MHC|HLA|antibody|humoral|immunoglobulin|defensin|granzyme|perforin|necrosis factor|NFAT|interleukin"
plot_keywords <- "MHC|antigen|myeloid|T cell selection|T cell activation|cytotoxicity|T cell proliferation|natural killer|interferon-alpha|interferon-beta|macrophage|interleukin-12 production"

# clusters to plot
target_clusters <- c(1, 2, 5)

### Load cluster membership
cluster_df <- read_tsv(cluster_file) %>%
  filter(grepl("up", expr_pattern))

### Function to build GO enrichment table and plot per cluster
generate_cluster_plot <- function(cluster_id) {
  
  message("Processing cluster ", cluster_id, " ...")
  
  # subset miRNAs for this cluster
  mirnas <- cluster_df %>%
    filter(row_cluster == cluster_id) %>%
    pull(mirna_id)
  
  # list GO enrichment result files for these miRNAs
  go_files <- list.files(results_dir, pattern = "target-go-enrichment", full.names = TRUE)
  cluster_go_files <- go_files[grepl(paste(mirnas, collapse = "|"), go_files)]
  
  if (length(cluster_go_files) == 0) {
    warning("No GO enrichment files found for cluster ", cluster_id)
    return(NULL)
  }
  
  # merge enrichment results
  merged_filtered_go_df <- map_dfr(cluster_go_files, function(file) {
    mirna <- str_remove(basename(file), "-target-go-enrichment.tsv")
    read_tsv(file, show_col_types = FALSE) %>%
      mutate(gene_ratio = Significant / Annotated,
             miRNA = mirna)
  })
  
  # filter for immune-related GO terms
  merged_filtered_go_df <- merged_filtered_go_df %>%
    filter(grepl(keywords, Term), nchar(Term) < 75)
  
  # reduced table for plotting
  plot_df <- merged_filtered_go_df %>%
    filter(grepl(plot_keywords, Term),
           !grepl("CD4|negative|leukocyte|positive regulation", Term)) %>%
    mutate(class = case_when(
      grepl("MHC|antigen", Term)        ~ "MHC biosynthesis/\nantigen presentation",
      grepl("T cell|natural killer", Term) ~ "T cell\nresponses",
      grepl("myeloid|macrophage", Term)   ~ "Myeloid cell\nresponses",
      grepl("interferon|interleukin", Term) ~ "Cytokine\nproduction",
      TRUE ~ "Other"
    )) %>%
    mutate(class = fct_relevel(class, rev(sort(unique(class))))) %>%
    arrange(class) %>%
    mutate(Term = fct_relevel(Term, sort(unique(Term)))) %>%
    mutate(miRNA = fct_relevel(miRNA, rev(mixedsort(unique(miRNA)))))
  
  # -----------------------------------------------------------------
  # Plot
  # -----------------------------------------------------------------
  p <- ggplot(plot_df, aes(x = miRNA, y = Term,
                           size = gene_ratio,
                           colour = -log10(classicFisher))) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(low = "lightpink", high = "darkred") +
    labs(
      x = NULL, y = NULL,
      size = "Gene Ratio",
      colour = "-log10(p)"
    ) +
    facet_grid(class ~ ., scales = "free_y", space = "free_y", switch = "y") +
    theme_Publication() +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 1),
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
    )
  
  p <- p + theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 130))
  
  # save plot
  out_plot <- file.path(plot_dir, sprintf("cluster-%s-immune-go-term-dotplot.pdf", cluster_id))
  ggsave(out_plot, p, width = 14, height = 9)
  
  # save merged and filtered tables
  write_tsv(merged_filtered_go_df, file.path(results_dir, sprintf("cluster%s-mirna-target-go-enr-immune-terms-full.tsv", cluster_id)))
  write_tsv(plot_df, file.path(results_dir, sprintf("cluster%s-mirna-target-go-enr-immune-terms-reduced.tsv", cluster_id)))
  
  message("Saved: ", out_plot)
}

### Run for clusters 1, 2, 5
walk(target_clusters, generate_cluster_plot)

### Session Info

sessionInfo()


