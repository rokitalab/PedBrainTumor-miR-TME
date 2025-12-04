# Plot MB cluster 1 miRNA target GO enrichment for immune related terms
#
# Ryan Corbett
#
# Dec 2025

# This script generates immune target-level dot plot indicating:
# 1) Pearson correlation coefficients with targeting cluster 1 miRNAs in MB
# 2) Log2-Fold Change of targets in MB versus normals

# Load packages
library(tidyverse)
library(ggplot2)
library(gtools)
library(ggrepel)

# set dir paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "human-mirna-target-correlations")
input_dir <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

source(file.path(root_dir, "figures", "theme.R"))

# wrangle data
res <- read_tsv(file.path(input_dir,
                          "mb-cluster1-mirna-immune-target-sig-interactions.txt")) %>% 
  dplyr::rename(`Target Gene` = Genes) %>%
  # separate genes with T cell and myeloid cell terms assigned
  separate_rows(category, sep = ";") %>%
  # take max log2FC and min padj from mb vs. healthy and mb vs. adjacent contrasts
  dplyr::mutate(max_log2FC = case_when(
    abs(mb_vs_adj_log2FC) > abs(mb_vs_healthy_log2FC) ~ mb_vs_adj_log2FC,
    TRUE ~ mb_vs_healthy_log2FC)) %>%
  dplyr::mutate(min_padj = case_when(
    mb_vs_adj_padj < mb_vs_healthy_padj ~ mb_vs_adj_padj,
    TRUE ~ mb_vs_healthy_padj
  )) %>%
  # Sort target genes for plotting
  dplyr::mutate(`Target Gene` = fct_relevel(`Target Gene`,
                                            rev(sort(unique(`Target Gene`)))
  )) %>%
  # assign target genes to DE group for plotting 
  dplyr::mutate(de_group = case_when(
    max_log2FC > 0 & min_padj < 0.05 ~ "MB Upregulated",
    max_log2FC < 0 & min_padj < 0.05 ~ "MB Downregulated",
    TRUE ~ "Not DE"
  ))

pdf(NULL)

# Filter for targets associated with T/NK cell processes 
tcell_res <- res %>%
  dplyr::filter(grepl("T Cell|NK Cell", category)) %>%
  droplevels()

# Generate dot plot
ggplot(tcell_res, aes(x = miRNA, y = `Target Gene`,
                      size = mirna_target_pearson_r,
                      fill = de_group)) +
  geom_point(shape = 21,          # circle with fill + border
             colour = "black",    # border color
             stroke = 0.3,        # border thickness
             alpha = 0.9) +
  # add asterisks for sig correlations
  geom_text(
    data = tcell_res %>% dplyr::filter(mirna_target_pearson_r < 0,
                                       mirna_target_pearson_p < 0.05),
    aes(x = miRNA, y = `Target Gene`, label = "*"),
    color = "grey70", 
    size = 5,         
    fontface = "bold",
    vjust = 0.75       
  ) +
  scale_fill_manual(values = c("MB Downregulated" = "red3", 
                               "Not DE" = "whitesmoke", 
                               "MB Upregulated" = "green4"),
                    guide = guide_legend(
                      override.aes = list(size = 5))) +
  labs(
    x      = NULL,
    size   = "miRNA-target\nPearson r",
    fill   = "Target Expression",
    y      = NULL
  ) +
  scale_size_continuous(
    trans = "reverse",
    range = c(2, 6),  
  ) +
  facet_grid(
    category ~ ., 
    scales = "free_y", 
    space  = "free_y", 
    switch = "y",
    labeller = labeller(category = label_wrap_gen(width = 18))
  ) +
  theme_Publication() +
  theme(
    strip.placement   = "outside", 
    strip.text.y.left = element_text(angle = 0),
    axis.text.x       = element_text(angle = 45, vjust = 1, hjust = 1)
  )

# save plot
ggsave(file.path(plot_dir, "mb-cluster1-target-t-cell-dotplot.pdf"),
       height = 6, width = 6)

# Filter for myeloid cell terms 
myeloid_cell_res <- res %>%
  dplyr::filter(grepl("Myeloid", category)) %>%
  droplevels()

# generate dot plot
ggplot(myeloid_cell_res, aes(x = miRNA, y = `Target Gene`,
                             size = mirna_target_pearson_r,
                             fill = de_group)) +
  geom_point(shape = 21,         
             colour = "black",   
             stroke = 0.3,  
             alpha = 0.9) +
  # add asterisks for sig correlations
  geom_text(
    data = myeloid_cell_res %>% dplyr::filter(mirna_target_pearson_r < 0,
                                              mirna_target_pearson_p < 0.05),
    aes(x = miRNA, y = `Target Gene`, label = "*"),
    color = "grey70",   
    size = 5,          
    fontface = "bold",
    vjust = 0.75  
  ) +
  scale_fill_manual(values = c("MB Downregulated" = "red3", 
                               "Not DE" = "whitesmoke", 
                               "MB Upregulated" = "green4"),
                    guide = guide_legend(
                      override.aes = list(size = 5))) +
  labs(
    x      = NULL,
    size   = "miRNA-target\npearson r",
    fill   = "Target Expression",
    y      = NULL
  ) +
  scale_size_continuous(
    trans = "reverse",
    range = c(2, 6),  # adjust min/max point size
  ) +
  facet_grid(
    category ~ ., 
    scales = "free_y", 
    space  = "free_y", 
    switch = "y",
    labeller = labeller(category = function(x) str_wrap(x, width = 18))
  ) +
  theme_Publication() +
  theme(
    strip.placement   = "outside", 
    strip.text.y.left = element_text(angle = 0),
    axis.text.x       = element_text(angle = 45, vjust = 1, hjust = 1)
  )

# Save plot
ggsave(file.path(plot_dir, "mb-cluster1-target-myeloid-cell-dotplot.pdf"),
       height = 10, width = 8)

# filter for other terms related to interleukin/interferon production
other_res <- res %>%
  dplyr::filter(grepl("Macrophage|Interleukin|Interferon", category)) %>%
  droplevels()

# Generate dot plot
ggplot(other_res, aes(x = miRNA, y = `Target Gene`,
                      size = mirna_target_pearson_r,
                      fill = de_group)) +
  geom_point(shape = 21,        
             colour = "black",   
             stroke = 0.3,  
             alpha = 0.9) +
  # add asterisks for sig correlations
  geom_text(
    data = other_res %>% dplyr::filter(mirna_target_pearson_r < 0,
                                       mirna_target_pearson_p < 0.05),
    aes(x = miRNA, y = `Target Gene`, label = "*"),
    color = "grey70",  
    size = 5,        
    fontface = "bold",
    vjust = 0.75  
  ) +
  scale_fill_manual(values = c("MB Downregulated" = "red3", 
                               "Not DE" = "whitesmoke", 
                               "MB Upregulated" = "green4"),
                    guide = guide_legend(
                      override.aes = list(size = 5))) +
  labs(
    x      = NULL,
    size   = "miRNA-target\npearson r",
    fill   = "Target Expression",
    y      = NULL
  ) +
  scale_size_continuous(
    trans = "reverse",
    range = c(2, 6),
  ) +
  facet_grid(
    category ~ ., 
    scales = "free_y", 
    space  = "free_y", 
    switch = "y",
    labeller = labeller(category = function(x) str_wrap(x, width = 18))
  ) +
  theme_Publication() +
  theme(
    strip.placement   = "outside", 
    strip.text.y.left = element_text(angle = 0),
    axis.text.x       = element_text(angle = 45, vjust = 1, hjust = 1)
  )

# Save plot
ggsave(file.path(plot_dir, "mb-cluster1-target-other-term-dotplot.pdf"),
       height = 7, width = 8)

# Print session info
sessionInfo()
