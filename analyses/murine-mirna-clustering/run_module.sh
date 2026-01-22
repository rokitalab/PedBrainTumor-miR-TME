#!/bin/bash
# Module author: Bicna Song
# 2026-01

set -e
set -o pipefail

# Run clustering and generate heatmaps
Rscript --vanilla 01-murine-de-mirna-clustering.R

# Plots DE miRNA expression across time points
Rscript --vanilla 02-plot-de-mirna.R


