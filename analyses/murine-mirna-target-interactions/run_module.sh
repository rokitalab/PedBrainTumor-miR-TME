#!/bin/bash
# -------------------------------------------------------------------------
# Script: run_module.sh
# Author: Bicna Song
# Date: 2025-11
# Description:
#   Run murine miRNA–target interaction analyses and generate cluster-specific
#   dot plots for immune-related target genes
# -------------------------------------------------------------------------

set -euo pipefail

# Generate miRNA–target interaction tables
Rscript --vanilla 01-murine-mirna-target-interactions.R
echo "[$(date)] miRNA–target interaction table generation completed successfully."

# Generate cluster-specific dot plots
for c in 1 2 3 4 5; do
  echo "[$(date)] Generating dot plots for cluster ${c}..."
  Rscript --vanilla 02-murine-mirna-target-dotplot.R \
    results/cluster${c}-mirna-immune-target-interactions.tsv
  echo "[$(date)] Dot plots for cluster ${c} completed successfully."
done

echo "[$(date)] All murine miRNA–target interaction analyses completed successfully."
