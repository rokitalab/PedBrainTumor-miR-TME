#!/bin/bash

set -e
set -o pipefail

# Run extracting 3′ UTR sequences from GENCODE v39
Rscript -e "rmarkdown::render('00-extract-3utr-gencode-v39.Rmd')"

# Run convert novel miRNA table to FASTA format
Rscript -e "rmarkdown::render('01-convert-all-novel-miRNA-to-fasta.Rmd')"

# Run extract DE Novel miRNA Sequences from FASTA
Rscript -e "rmarkdown::render('02-extract-de-novel-miRNAs-to-fasta.Rmd')"

# Run miranda for novel miRNA target prediction
bash  03-run-miranda.sh 

