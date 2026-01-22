#!/bin/bash
# -------------------------------------------------------------------------
# Script: run_module.sh
# Author: Bicna Song
# Date: 2025-12
# Description:
#   Assess neuronal gene expression in tumor vs. normal brain samples across DIPG/DMG, MB, EPN, LGG
# -------------------------------------------------------------------------

set -euo pipefail

# Run neuronal expression workflow
Rscript -e "rmarkdown::render('01-human-neuronal-expression.Rmd')"
echo "[$(date)] Neuronal expression analysis completed successfully."

# Run correlation analysis
Rscript -e "rmarkdown::render('02-human-neuronal-mirna-target-correlation.Rmd')"
echo "[$(date)] miRNA–neuronal pathway correlation analysis completed successfully."

# Run clustering and generate heatmap
Rscript --vanilla 03-human-neuronal-heatmap.R
echo "[$(date)] Clustering and heatmap generation completed successfully."
