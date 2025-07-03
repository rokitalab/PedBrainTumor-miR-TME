#!/bin/bash

set -e
set -o pipefail

# Run extracting 3′ UTR sequences from GENCODE v39
Rscript 00-extract-3utr-gencode-v39.R

# Run convert novel miRNA table to FASTA format
Rscript -e "rmarkdown::render('01-convert-all-novel-miRNA-to-fasta.Rmd')"

# Run miranda for novel miRNA target prediction
bash  02-run-miranda.sh

# Run python script to parse miranda output file
python 03-parse-miranda-output.py 

# Run annotate miranda output with Gene ID and Gene Symbol
Rscript -e "rmarkdown::render('04-annotate-parsed-miranda-output.Rmd')"

