#!/bin/bash

set -e
set -o pipefail

# Run merge DE results script
Rscript -e "rmarkdown::render('01-merge-mirna-rna-differential-expression.Rmd')"