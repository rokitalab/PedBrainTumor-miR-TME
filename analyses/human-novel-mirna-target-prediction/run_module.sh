#!/bin/bash

set -e
set -o pipefail

# Run extracting 3′ UTR sequences from GENCODE v39
Rscript -e "rmarkdown::render('00-extract-3utr-gencode-v39.Rmd')"
