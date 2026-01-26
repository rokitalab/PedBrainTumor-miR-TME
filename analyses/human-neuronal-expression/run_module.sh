#!/bin/bash
# -------------------------------------------------------------------------
# Script: run_module.sh
# Author: Bicna Song
# Date: Jan 2026
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

# Plot cluster-specific miRNA target GO enrichment
Rscript --vanilla 04-human-neuronal-go-enrichment.R
echo "[$(date)] Plotting cluster-specific miRNA target GO enrichment dot plot completed successfully."

# Identify miRNA–target interactions for DIPG/DMG cluster 6
Rscript --vanilla 05-dmg-cluster6-mirna-target-interactions.R
echo "[$(date)] DIPG/DMG cluster 6 miRNA–target interaction analysis completed successfully."

# Plot DIPG/DMG cluster 6 miRNA target dot plots
Rscript --vanilla 06-dmg-cluster6-mirna-target-dotplot.R
echo "[$(date)] DIPG/DMG cluster 6 miRNA target dot plots generated successfully."

# Identify miRNA–target interactions for MB cluster 1
Rscript --vanilla 07-mb-cluster1-mirna-target-interactions.R
echo "[$(date)] MB cluster 1 miRNA–target interaction analysis completed successfully."

# Plot MB cluster 1 miRNA target dot plots
Rscript --vanilla 08-mb-cluster1-mirna-target-dotplot.R
echo "[$(date)] MB cluster 1 miRNA target dot plots generated successfully."

echo "[$(date)] Human neuronal expression module finished successfully."


