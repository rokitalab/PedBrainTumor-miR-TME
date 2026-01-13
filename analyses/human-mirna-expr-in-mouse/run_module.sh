#!/bin/bash
# Module author: Ryan Corbett
# 2026-01

set -e
set -o pipefail

# identify orthologs and expression patterns in mice
R -e "rmarkdown::render('01-summarize-orthology-expression.Rmd')"

# plot DMG cluster 6 DE miRNA expression 
Rscript --vanilla 02-plot-cluster6-ortholog-expr.R