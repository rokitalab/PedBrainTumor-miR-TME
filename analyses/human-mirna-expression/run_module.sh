#!/bin/bash

set -e
set -o pipefail

# Run merge miRNA counts script
Rscript -e "rmarkdown::render('01-merge-miRNA-counts.Rmd')"

# Run PCA and UMAP script
Rscript -e "rmarkdown::render('02-mirna-pca-umap.Rmd')"

# Run miRNA DE analysis script
Rscript -e "rmarkdown::render('03-mirna-differential-expression-analysis.Rmd')"
