#!/bin/bash

set -e
set -o pipefail

# Run differential expression analysis
Rscript -e "rmarkdown::render('01-mirna-differential-expression.Rmd')"

# Perform GO enrichment analysis on DE miRNA targets
Rscript -e "rmarkdown::render('02-GO_enrichment.Rmd')"

# Generate summary plots for GO enrichment
Rscript -e "rmarkdown::render('03-plot_gsea.Rmd')"

# Generate immune-related GO term dot plots for clusters 1, 2, and 5
Rscript --vanilla 04-cluster-go-enrichment.R

