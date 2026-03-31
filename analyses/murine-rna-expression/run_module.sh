#!/bin/bash

set -e
set -o pipefail

# Run differential expression
Rscript -e "rmarkdown::render('01-murine-rna-differential-expression-analysis.Rmd')"

# Run GO enrichment 
Rscript -e "rmarkdown::render('02-GO-enrichment.Rmd')"

# Identify downregulated miRNAS targeting immune genes
Rscript -e "rmarkdown::render('03-mirna-target-correlations.Rmd')"