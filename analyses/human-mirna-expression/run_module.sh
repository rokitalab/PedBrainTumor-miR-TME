#!/bin/bash

set -e
set -o pipefail

# Run merge miRNA counts script
Rscript -e "rmarkdown::render('01-merge-miRNA-counts.Rmd')"

# Run PCA batch comparison script
Rscript -e "rmarkdown::render('02-pca-mirna-batch-comparison.Rmd')"