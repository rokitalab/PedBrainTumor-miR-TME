#!/bin/bash

set -e
set -o pipefail

# Run merge miRNA counts script
Rscript -e "rmarkdown::render('01-merge-miRNA-counts.Rmd')"

# Run PCA batch comparison script
Rscript -e "rmarkdown::render('02-pca-mirna-batch-comparison.Rmd')"

# Run PCA by histology script
Rscript -e "rmarkdown::render('03-pca-by-histology.Rmd')"

# Run miRNA DE analysis script
Rscript -e "rmarkdown::render('04-mirna-differential-expression-analysis.Rmd')"

# Run miRNA UMAP script
Rscript -e "rmarkdown::render('05-umap-mirna.Rmd')"
