#!/bin/bash
# Module author: Ryan Corbett
# 2026-01

set -e
set -o pipefail

# identify orthologs and expression patterns in mice
R -e "rmarkdown::render('01-summarize-orthology-expression.Rmd')"

# plot DMG cluster 6 DE miRNA expression 
Rscript --vanilla 02-plot-cluster6-ortholog-expr.R

# identify and plot predicted targets of DMG cluster 6 DE miRNA orthologs
Rscript --vanilla 03-plot-cluster6-ortholog-target-expr.R

# test direction of B7H3-associated DE among immune-related targets of DMG cluster 6 DE miRNA orthologs
Rscript --vanilla 04-immune-target-expression.R
