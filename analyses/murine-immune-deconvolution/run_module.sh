#!/bin/bash
# -------------------------------------------------------------------------
# Script: run_module.sh
# Author: Bicna Song
# Date: 2025-10
# Description:
#   Run immune deconvolution of mouse bulk RNA-seq data using xCell2.
# -------------------------------------------------------------------------

set -euo pipefail

# Run XCell2 workflow
Rscript --vanilla 01-run-xcell2.R

# Run xCell2 Scaled Enrichment Heatmaps by Reference
Rscript --vanilla 02-xcell2-heatmap.R

# Run xCell2 plot
R -e "rmarkdown::render('03-xcell2-plot.Rmd')"
