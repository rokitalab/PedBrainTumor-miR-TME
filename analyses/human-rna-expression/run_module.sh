#!/bin/bash

set -e
set -o pipefail

# Run merge RNA counts script
Rscript -e "rmarkdown::render('01-merge-and-annotate-RNAseq-STAR-counts.Rmd')"