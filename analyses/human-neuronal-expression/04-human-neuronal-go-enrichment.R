# Plot cluster-specific miRNA target GO enrichment for neuronal-related terms
# The script has been adpated from /human-mirna-go-enrichment/04-cluster6-go-enrichment.R
# Author: Bicna Song | Jan 2026

### Setup
library(tidyverse)
library(ggplot2)
library(gtools)
library(rprojroot)

### Paths
root_dir <- find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "human-neuronal-expression")
input_dir <- file.path(analysis_dir, "input")
results_dir <- file.path(analysis_dir, "results")
plot_dir <- file.path(analysis_dir, "plots")

source(file.path(root_dir, "figures", "theme.R"))

### Analysis configuration
cluster_configs <- tribble(
  ~tumor,        ~cluster_id, ~direction, ~cluster_label, ~cluster_file,
  "DIPG or DMG", 6,           "up",       "cluster6",
  file.path(results_dir,
            "DIPG_or_DMG_mirna_neuronal_cluster_membership.tsv"),
  "MB",          1,           "up",       "cluster1",
  file.path(results_dir,
            "MB_mirna_neuronal_cluster_membership.tsv")
)

### Load neuronal-related GO BP terms (mapped via GO ID)
## Rationale: GSVA GO BP terms and GO enrichment output use different pathway name formats.
## GO IDs provide a stable identifier to link the two analyses.
neuronal_go_ids <- read_tsv(
  file.path(input_dir, "neuronal_go_term_with_id.tsv"),
  show_col_types = FALSE
) %>%
  filter(
    Pathway == "GOBP_NEURON_CELL_CELL_ADHESION" |
      !grepl("CELL_ADHESION", Pathway)
  )
  
neuronal_go_terms <- read_tsv(
  file.path(results_dir, "neuronal_gobp_terms_from_gsva.tsv"),
  show_col_types = FALSE
) %>%
  left_join(neuronal_go_ids, by = "Pathway") %>%
  filter(!is.na(GO.ID), GO.ID != "") %>%
  distinct(GO.ID) %>%
  pull(GO.ID)

### List all GO enrichment files
all_go_files <- list.files(
  file.path(root_dir, "analyses", "human-mirna-go-enrichment", "results")
)

### Main loop: run each tumor / cluster analysis
for (i in seq_len(nrow(cluster_configs))) {
  
  tumor <- cluster_configs$tumor[i]
  tumor_safe <- str_replace_all(tumor, "\\s+or\\s+", "_or_")
  cluster_id <- cluster_configs$cluster_id[i]
  direction <- cluster_configs$direction[i]
  cluster_lab <- cluster_configs$cluster_label[i]
  cluster_file <- cluster_configs$cluster_file[i]
  
  message("Processing ", tumor,
          " | cluster ", cluster_id,
          " (", direction, ")")
  
  ## Load cluster membership
  cluster_df <- read_tsv(cluster_file, show_col_types = FALSE)
  
  cluster_mirnas <- cluster_df %>%
    filter(row_cluster == cluster_id, direction == direction) %>%
    pull(miRNA)
  
  ## Identify GO enrichment files
  cluster_go_files <- all_go_files %>%
    keep(~ grepl(tumor, .x)) %>%
    keep(~ grepl(paste(cluster_mirnas, collapse = "|"), .x))
  
  ## Load and merge GO enrichment
  merged_filtered_go_df <- map_dfr(cluster_go_files, function(file) {
    
    mirna <- str_remove(
      file,
      paste0("-", tumor, "-shared-GO-BP_all.tsv")
    )
    
    read_tsv(
      file.path(
        root_dir,
        "analyses",
        "human-mirna-go-enrichment",
        "results",
        file
      ),
      show_col_types = FALSE
    ) %>%
      filter(
        !is.na(GO.ID),
        GO.ID %in% neuronal_go_terms
      ) %>%
      mutate(
        gene_ratio = Significant / Annotated,
        miRNA = mirna
      )
  })
  
  ## Define plotting data frame
  ##
  ## - remove regulation-only terms
  ## - exclude purinergic / adrenergic pathways due to
  ##   small numbers of significant enriched pathways
  ## - assign broad neuronal functional classes
  plot_df <- merged_filtered_go_df %>%
    filter(!grepl("negative|positive regulation",
                  Term, ignore.case = TRUE)) %>%
    filter(!grepl("purinergic|adrenergic",
                  Term, ignore.case = TRUE)) %>%
    mutate(class = case_when(
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
    filter(!is.na(class)) %>%
    mutate(
      class = fct_relevel(class, rev(sort(unique(class)))),
      Term  = fct_relevel(Term, sort(unique(Term))),
      miRNA = fct_relevel(miRNA, rev(mixedsort(unique(miRNA))))
    ) %>%
    arrange(class)
  
  ## Generate dot plot
  go_dotplot <- plot_df %>%
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
  
  ## Save outputs
  ggsave(
    file.path(
      plot_dir,
      paste0(tumor_safe, "-", cluster_lab,
             "-neuronal-go-term-dotplot.pdf")
    ),
    go_dotplot,
    width  = 14,
    height = 14
  )
  
  write_tsv(
    merged_filtered_go_df,
    file.path(
      results_dir,
      paste0(tumor_safe, "-", cluster_lab,
             "-mirna-target-go-enr-neuronal-terms-full.tsv")
    )
  )
  
  plot_df_export <- plot_df %>%
    mutate(class = str_replace_all(class, "\n", " "))
  
  write_tsv(
    plot_df_export,
    file.path(
      results_dir,
      paste0(tumor_safe, "-", cluster_lab,
             "-mirna-target-go-enr-neuronal-terms-reduced.tsv")
    )
  )
}

### Session Info
sessionInfo()


