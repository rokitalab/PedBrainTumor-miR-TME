#!/bin/bash

set -e
set -o pipefail

# Run differential expression analysis
Rscript -e "rmarkdown::render('01-mirna-differential-expression.Rmd')"
