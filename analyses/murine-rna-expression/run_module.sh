#!/bin/bash

set -e
set -o pipefail

# Run differential expression
Rscript -e "rmarkdown::render('01-murine-rna-differential-expression-analysis.Rmd')"
