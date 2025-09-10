#!/bin/bash
# Module author: Ryan Corbett
# 2025-09

set -e
set -o pipefail

# Run clustering and generate heatmaps
Rscript --vanilla 01-de-mirna-clustering.R