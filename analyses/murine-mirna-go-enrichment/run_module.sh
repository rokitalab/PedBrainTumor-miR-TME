#!/bin/bash

set -e
set -o pipefail

# Perform GO enrichment analysis on DE miRNA targets
Rscript -e "rmarkdown::render('01-GO_enrichment.Rmd')"

# Generate summary plots for GO enrichment
Rscript -e "rmarkdown::render('02-plot_gsea.Rmd')"

# Generate summary GO enrichment plots by miRNA cluster
Rscript --vanilla 03-cluster-go-enrichment.R

