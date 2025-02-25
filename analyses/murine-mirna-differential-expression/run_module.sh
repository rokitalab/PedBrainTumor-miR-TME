#!/bin/bash

set -e
set -o pipefail

# Run DE analyses
Rscript -e "rmarkdown::render('01-mirna-differential-expression.Rmd')"