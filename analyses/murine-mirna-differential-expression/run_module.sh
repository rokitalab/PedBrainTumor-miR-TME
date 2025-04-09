#!/bin/bash

set -e
set -o pipefail

# Run DE analyses
Rscript -e "rmarkdown::render('01-mirna-differential-expression.Rmd')"

# Run GO enrichment on DE miRNA targets
Rscript -e "rmarkdown::render('02-GO_enrichment.Rmd')"