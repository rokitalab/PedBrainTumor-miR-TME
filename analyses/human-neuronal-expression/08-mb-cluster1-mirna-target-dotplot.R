# Plot MB cluster 1 miRNA target GO enrichment for neuronal related terms
# Ryan Corbett | Oct 2025
# Updated by Bicna Song | Jan 2026 

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
analysis_dir <- file.path(root_dir, "analyses", "human-neuronal-expression")
input_dir <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

source(file.path(root_dir, "figures", "theme.R"))

# wrangle data
res <- read_tsv(file.path(results_dir,
                          "MB-cluster1-mirna-immune-target-sig-interactions.tsv")) %>% 
  mutate(category_broad = case_when(
    grepl("postsynaptic|scaffolding",
          Term, ignore.case = TRUE) ~
      "postsynaptic\nscaffolding",
    grepl("synaptic|vesicle|trafficking",
          Term, ignore.case = TRUE) &
      !grepl("postsynaptic",
             Term, ignore.case = TRUE) ~
      "synaptic vesicle\ntrafficking",
    grepl("neurotransmitter|receptors",
          Term, ignore.case = TRUE) ~
      "neurotransmitter\nreceptors",
    grepl("cell|adhesion",
          Term, ignore.case = TRUE) ~
      "cell\nadhesion",
    grepl("neuromodulatory",
          Term, ignore.case = TRUE) ~
      "neuromodulatory",
    grepl("glutamatergic",
          Term, ignore.case = TRUE) ~
      "glutamatergic",
    grepl("gabaergic",
          Term, ignore.case = TRUE) ~
      "gabaergic",
    grepl("cholinergic",
          Term, ignore.case = TRUE) ~
      "cholinergic",
    grepl("serotonergic",
          Term, ignore.case = TRUE) ~
      "serotonergic",
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(category_broad)) %>%
  dplyr::rename(`Target Gene` = Genes) %>%
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

# Filter for targets associated with cell adhesion processes 
cell_adhesion_res <- res %>%
  dplyr::filter(grepl("cell|adhesion",
                      category_broad, ignore.case = TRUE)) %>%
  droplevels()

# Generate dot plot
ggplot(cell_adhesion_res, aes(x = miRNA, y = `Target Gene`,
                      size = mirna_target_pearson_r,
                      fill = de_group)) +
  geom_point(shape = 21,          # circle with fill + border
             colour = "black",    # border color
             stroke = 0.3,        # border thickness
             alpha = 0.9) +
  # add asterisks for sig correlations
  geom_text(
    data = cell_adhesion_res %>% dplyr::filter(mirna_target_pearson_r < 0,
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
    category_broad ~ ., 
    scales = "free_y", 
    space  = "free_y", 
    switch = "y",
    labeller = labeller(category_broad = label_wrap_gen(width = 18))
  ) +
  theme_Publication() +
  theme(
    strip.placement   = "outside", 
    strip.text.y.left = element_text(angle = 0),
    axis.text.x       = element_text(angle = 45, vjust = 1, hjust = 1)
  )

# save plot
ggsave(file.path(plot_dir, "MB-cluster1-target-cell-adhesion-dotplot.pdf"),
       height = 26, width = 10)

# Filter for postsynaptic-related processes
postsynaptic_res <- res %>%
  dplyr::filter(grepl("postsynaptic|scaffolding", category_broad, ignore.case = TRUE)) %>%
  droplevels()

# generate dot plot
ggplot(postsynaptic_res, aes(x = miRNA, y = `Target Gene`,
                             size = mirna_target_pearson_r,
                             fill = de_group)) +
  geom_point(shape = 21,         
             colour = "black",   
             stroke = 0.3,  
             alpha = 0.9) +
  # add asterisks for sig correlations
  geom_text(
    data = postsynaptic_res %>% dplyr::filter(mirna_target_pearson_r < 0,
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
    category_broad ~ ., 
    scales = "free_y", 
    space  = "free_y", 
    switch = "y",
    labeller = labeller(category_broad = function(x) str_wrap(x, width = 18))
  ) +
  theme_Publication() +
  theme(
    strip.placement   = "outside", 
    strip.text.y.left = element_text(angle = 0),
    axis.text.x       = element_text(angle = 45, vjust = 1, hjust = 1)
  )

# Save plot
ggsave(file.path(plot_dir, "MB-cluster1-target-postsynaptic-dotplot.pdf"),
       height = 16, width = 10)

# Filter for synaptic-related processes
synaptic_res <- res %>%
  dplyr::filter(grepl("synaptic|vesicle|trafficking", category_broad, ignore.case = TRUE) &
                  !grepl("postsynaptic", category_broad, ignore.case = TRUE)) %>%
  droplevels()

# generate dot plot
ggplot(synaptic_res, aes(x = miRNA, y = `Target Gene`,
                             size = mirna_target_pearson_r,
                             fill = de_group)) +
  geom_point(shape = 21,         
             colour = "black",   
             stroke = 0.3,  
             alpha = 0.9) +
  # add asterisks for sig correlations
  geom_text(
    data = synaptic_res %>% dplyr::filter(mirna_target_pearson_r < 0,
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
    category_broad ~ ., 
    scales = "free_y", 
    space  = "free_y", 
    switch = "y",
    labeller = labeller(category_broad = function(x) str_wrap(x, width = 18))
  ) +
  theme_Publication() +
  theme(
    strip.placement   = "outside", 
    strip.text.y.left = element_text(angle = 0),
    axis.text.x       = element_text(angle = 45, vjust = 1, hjust = 1)
  )

# Save plot
ggsave(file.path(plot_dir, "MB-cluster1-target-synaptic-dotplot.pdf"),
       height = 14, width = 10)

# Filter for neurotransmitter-related processes
neurotransmitter_res <- res %>%
  dplyr::filter(grepl("neurotransmitter|receptors", category_broad, ignore.case = TRUE)) %>%
  droplevels()

# generate dot plot
ggplot(neurotransmitter_res, aes(x = miRNA, y = `Target Gene`,
                         size = mirna_target_pearson_r,
                         fill = de_group)) +
  geom_point(shape = 21,         
             colour = "black",   
             stroke = 0.3,  
             alpha = 0.9) +
  # add asterisks for sig correlations
  geom_text(
    data = neurotransmitter_res %>% dplyr::filter(mirna_target_pearson_r < 0,
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
    category_broad ~ ., 
    scales = "free_y", 
    space  = "free_y", 
    switch = "y",
    labeller = labeller(category_broad = function(x) str_wrap(x, width = 18))
  ) +
  theme_Publication() +
  theme(
    strip.placement   = "outside", 
    strip.text.y.left = element_text(angle = 0),
    axis.text.x       = element_text(angle = 45, vjust = 1, hjust = 1)
  )

# Save plot
ggsave(file.path(plot_dir, "MB-cluster1-target-neurotransmitter-dotplot.pdf"),
       height = 12, width = 10)

# Print session info
sessionInfo()
