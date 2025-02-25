#!/bin/bash

set -e
set -o pipefail

# Run DE analysis script
Rscript -e "rmarkdown::render('01-miRNA-differential-expression.Rmd')"