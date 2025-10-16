#!/bin/bash
# Module author: Bicna Song
# 2025-10

set -e
set -o pipefail

# Run clustering and generate heatmaps
Rscript --vanilla 01-murine-de-mirna-clustering.R