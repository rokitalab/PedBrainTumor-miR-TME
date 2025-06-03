#!/bin/bash

set -e
set -o pipefail

# Run RNA-seq PCA script
Rscript -e "rmarkdown::render('rna_seq_pca.Rmd')"

# Run RNA-seq UMAP script
Rscript -e "rmarkdown::render('rna_seq_umap.Rmd')"

