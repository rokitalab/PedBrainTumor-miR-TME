# Plot cluster 6 miRNA target GO enrichment for immune related terms
#
# Ryan Corbett
#
# Oct 2025

# This script loads and merges cluster 6 miRNA GO enriched terms, filters for pathways of interest, and generates a dot plot

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
                          "dipg-dmg-cluster6-mirna-immune-target-sig-interactions.txt")) %>% 
  # filter for genes without assigned terms
  dplyr::filter(!is.na(category_broad) | !is.na(category_narrow)) %>%
  # separate genes with T Cell and Myeloid cell terms assigned
  separate_rows(category_narrow, sep = ",") %>%
  # take max log2FC and min padj from dmg vs. healthy and dmg vs. adjacent contrasts
  dplyr::mutate(max_log2FC = case_when(
    abs(log2FC_tumor_vs_adj_normal) > abs(log2FC_tumor_vs_healthy_normal) ~ log2FC_tumor_vs_adj_normal,
    TRUE ~ log2FC_tumor_vs_healthy_normal)) %>%
  dplyr::mutate(min_padj = case_when(
    padj_tumor_vs_adj_normal < padj_tumor_vs_healthy_normal ~ padj_tumor_vs_adj_normal,
    TRUE ~ padj_tumor_vs_healthy_normal
  )) %>%
  # Sort target genes for plotting
  dplyr::mutate(`Target Gene` = fct_relevel(`Target Gene`,
                                            rev(sort(unique(`Target Gene`)))
  )) %>%
  # assign target genes to DE group for plotting 
  dplyr::mutate(de_group = case_when(
    max_log2FC > 0 & min_padj < 0.05 ~ "DMG Upregulated",
    max_log2FC < 0 & min_padj < 0.05 ~ "DMG Downregulated",
    TRUE ~ "Not DE"
  ))

# Filter for targets associated with T cell processes 
tcell_res <- res %>%
  dplyr::filter(grepl("T Cell|NK Cell", category_narrow)) %>%
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
  scale_fill_manual(values = c("DMG Downregulated" = "red3", 
                               "Not DE" = "whitesmoke", 
                               "DMG Upregulated" = "green4"),
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
    category_narrow ~ ., 
    scales = "free_y", 
    space  = "free_y", 
    switch = "y"
  ) +
  theme_Publication() +
  theme(
    strip.placement   = "outside", 
    strip.text.y.left = element_text(angle = 0),
    axis.text.x       = element_text(angle = 45, vjust = 1, hjust = 1)
  )

# save plot
ggsave(file.path(plot_dir, "t-cell-dotplot.pdf"),
       height = 8, width = 10)

# Filter for myeloid cell processes 
myeloid_cell_res <- res %>%
  dplyr::filter(grepl("Myeloid", category_narrow)) %>%
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
  scale_fill_manual(values = c("DMG Downregulated" = "red3", 
                               "Not DE" = "whitesmoke", 
                               "DMG Upregulated" = "green4"),
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
    category_narrow ~ ., 
    scales = "free_y", 
    space  = "free_y", 
    switch = "y",
    labeller = labeller(category_narrow = function(x) str_wrap(x, width = 30))
  ) +
  theme_Publication() +
  theme(
    strip.placement   = "outside", 
    strip.text.y.left = element_text(angle = 0),
    axis.text.x       = element_text(angle = 45, vjust = 1, hjust = 1)
  )

# Save plot
ggsave(file.path(plot_dir, "myeloid-cell-dotplot.pdf"),
       height = 12, width = 12)

# filter for other terms related to antigen presentation and interleukin/interferon production
other_res <- res %>%
  dplyr::filter(grepl("Antigen|Interleukin|Interferon", category_narrow)) %>%
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
  scale_fill_manual(values = c("DMG Downregulated" = "red3", 
                               "Not DE" = "whitesmoke", 
                               "DMG Upregulated" = "green4"),
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
    category_narrow ~ ., 
    scales = "free_y", 
    space  = "free_y", 
    switch = "y",
    labeller = labeller(category_narrow = function(x) str_wrap(x, width = 30))
  ) +
  theme_Publication() +
  theme(
    strip.placement   = "outside", 
    strip.text.y.left = element_text(angle = 0),
    axis.text.x       = element_text(angle = 45, vjust = 1, hjust = 1)
  )

# Save plot
ggsave(file.path(plot_dir, "other-term-dotplot.pdf"),
       height = 5, width = 8)

# Print session info
sessionInfo()
