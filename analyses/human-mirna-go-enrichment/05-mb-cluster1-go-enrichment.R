# Plot MB cluster 1 miRNA target GO enrichment for immune related terms
#
# Ryan Corbett
#
# Oct 2025

# This script loads and merges MB cluster 1 miRNA GO enriched terms, filters for pathways of interest, and generates a dot plot

# Load packages
library(tidyverse)
library(ggplot2)
library(gtools)

# set dir paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "human-mirna-go-enrichment")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

source(file.path(root_dir, "figures", "theme.R"))


# set file paths
mb_cluster_file <- file.path(root_dir, "analyses",
                              "mirna-clustering",
                              "results",
                              "MB-de-mirna-cluster-membership-immune-scores.tsv")


# wrangle data
mb_cluster_df <- read_tsv(mb_cluster_file)

# filter for cluster 1 miRNAs
cluster_1_mirnas <- mb_cluster_df %>%
  dplyr::filter(row_cluster == 1) %>%
  pull(miRNA)

# Define key words to use for filtering GO enrichment results
keywords <- "inflamm|immunity|T cell|natural killer cell|JAK|immune|leukocyte|lymphocyte|myeloid|macrophage|transforming growth factor beta|MAPK cascade|defense|cytokine |chemokine|interferon|glycolysis|lactate|antigen|MHC|HLA|antibody|humoral|immunoglobulin|defensin|granzyme|perforin|necrosis factor|NFAT|interleukin"

# List all GO enrichment files and filter for cluster 1 miRNAs (upregulated in MB)
go_files <- list.files(file.path(root_dir, "analyses",
                                 "human-mirna-go-enrichment",
                                 "results"))
mb_go_files <- go_files[grepl("MB", go_files)]
mb_clust1_go_files <- mb_go_files[grepl(paste(cluster_1_mirnas, collapse = "|"), mb_go_files)]

# loop through go enrichment files and filter for terms of interest
for (file in mb_clust1_go_files){
  
  mirna <- str_remove_all(file, "-MB-shared-GO-BP_all.tsv")
  
  if (file == mb_clust1_go_files[1]){
    
    merged_filtered_go_df <- read_tsv(file.path(root_dir, "analyses",
                                                "human-mirna-go-enrichment",
                                                "results", file)) %>%
      dplyr::filter(grepl(keywords, Term),
                    # remove terms with excessive number of characters (for plotting)
                    nchar(Term) < 75) %>%
      dplyr::mutate(gene_ratio = Significant/Annotated) %>%
      dplyr::mutate(miRNA = mirna)
    
  } else {
    
    filtered_go_df <- read_tsv(file.path(root_dir, "analyses",
                                         "human-mirna-go-enrichment",
                                         "results", file)) %>%
      dplyr::filter(grepl(keywords, Term),
                    nchar(Term) < 75) %>%
      dplyr::mutate(gene_ratio = Significant/Annotated) %>%
      dplyr::mutate(miRNA = mirna)
    
    # merge results
    merged_filtered_go_df <- merged_filtered_go_df %>%
      bind_rows(filtered_go_df)
    
  }
  
}

# define key words used to further filter results for plotting
plot_keywords <- "MHC|myeloid|T cell selection|T cell activation|cytotoxicity|T cell proliferation|natural killer|interferon-alpha|interferon-beta|macrophage|interleukin-12 production|interferon|interleukin"

# define plotting df
plot_df <- merged_filtered_go_df %>%
  dplyr::filter(grepl(plot_keywords, Term),
                # remove terms indicating positive/negative regulation to reduce N
                !grepl("CD4|negative|leukocyte|positive regulation", Term)) %>%
  # Define broad classes for Terms to group for plotting
  dplyr::mutate(class = case_when(
    grepl("MHC|antigen", Term) ~ "Antigen\npresentation/\nsignaling",
    grepl("T cell|natural killer", Term) ~ "T cell\nresponses",
    grepl("myeloid|macrophage", Term) ~ "Myeloid cell\nresponses",
    grepl("interferon|interleukin", Term) ~ "Cytokine\nproduction"
  )) %>% 
  # order classes, terms, and miRNAs
  dplyr::mutate(class = fct_relevel(class,
                                    rev(sort(unique(class))))) %>%
  dplyr::arrange(class) %>%
  dplyr::mutate(Term = fct_relevel(Term, sort(unique(Term)))) %>%
  dplyr::mutate(miRNA = fct_relevel(miRNA,
                                    rev(mixedsort(unique(miRNA)))))

# generate dot plot
plot_df %>% 
  ggplot(aes(x = miRNA, y = Term,
             size = gene_ratio,
             colour = -log10(classicFisher))) +
  geom_point(alpha = 0.9) +
  scale_color_gradient(low = "lightpink", high = "darkred") +
  labs(
    x      = NULL,
    size   = "Gene Ratio",
    colour = "-log10(p)",
    y      = NULL
  ) +
  facet_grid(class ~ ., scales = "free_y", 
             space = "free_y", switch = "y") +
  theme_Publication() +
  theme(
    strip.placement = "outside",                     # move strips outside panel
    strip.text.y.left = element_text(angle = 90),    # rotate facet labels
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
  )

# save plot
ggsave(file.path(plot_dir,
                 "MB-cluster-1-immune-go-term-dotplot.pdf"),
       plot = p, width = 11, height = 8)


# save full and reduced merged GO enrichment tables
write_tsv(merged_filtered_go_df,
          file.path(results_dir,
                    "MB-cluster1-mirna-target-go-enr-immune-terms-full.tsv"))

write_tsv(plot_df,
          file.path(results_dir,
                    "MB-cluster1-mirna-target-go-enr-immune-terms-reduced.tsv"))
