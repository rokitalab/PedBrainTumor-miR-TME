#!/bin/bash

set -e
set -o pipefail

# PCA
Rscript -e "rmarkdown::render('01-rna-pca.Rmd')"

# Run differential expression
Rscript -e "rmarkdown::render('02-differential-expression.Rmd')"

# Run GO enrichment
Rscript -e "rmarkdown::render('03-GO-enrichment.Rmd')"

# Run GSEA enrichment
Rscript -e "rmarkdown::render('04-GSEA-hallmark.Rmd')"

# Run GSVA analysis
Rscript -e "rmarkdown::render('05-GSVA-GOBP.Rmd')"

# Run miRNA-GSVA-corrleation analysis
Rscript -e "rmarkdown::render('06-miRNA-GSVA-correlation.Rmd')"
