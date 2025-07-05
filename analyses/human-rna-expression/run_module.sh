#!/bin/bash

set -e
set -o pipefail

# PCA/UMAP
#Rscript -e "rmarkdown::render('01-merge-miRNA-counts.Rmd')"

# Run differential expression
Rscript -e "rmarkdown::render('02-differential-expression.Rmd')"

# Run GO enrichment
Rscript -e "rmarkdown::render('03-GO-enrichment.Rmd')"

# Run GSEA enrichment
Rscript -e "rmarkdown::render('04-GSEA-hallmark.Rmd')"
