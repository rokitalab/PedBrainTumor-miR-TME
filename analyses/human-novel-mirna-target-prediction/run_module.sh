#!/bin/bash

set -e
set -o pipefail

# Run extracting 3′ UTR sequences from GENCODE v39
Rscript -e "rmarkdown::render('00-extract-3utr-gencode-v39')"

# Run convert novel miRNA table to FASTA format
Rscript -e "rmarkdown::render('01-convert-all-novel-miRNA-to-fasta.Rmd')"
