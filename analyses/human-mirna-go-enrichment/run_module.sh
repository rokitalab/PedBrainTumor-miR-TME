#!/bin/bash

set -e
set -o pipefail

# Run merge DE results script
Rscript -e "rmarkdown::render('01-merge-mirna-rna-differential-expression.Rmd')"

# Run DE miRNA target GO enrichment
Rscript -e "rmarkdown::render('02-GO_enrichment.Rmd')"

# Plot GO enrichment
Rscript -e "rmarkdown::render('03-plot-enriched-terms.Rmd')"

# Plot GO enrichment in DIPG/DMG cluster 6 targets
Rscript --vanilla 04-cluster6-go-enrichment.R

# Plot GO enrichment in MB cluster 1 targets
Rscript --vanilla 05-mb-cluster1-go-enrichment.R