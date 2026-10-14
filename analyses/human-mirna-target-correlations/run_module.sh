#!/bin/bash
# Module author: Ryan Corbett
# 2025-10

set -e
set -o pipefail

# Identify significant miRNA - immune target assocations among DMG cluster6 miRNAs
Rscript --vanilla 01-dmg-cluster6-mirna-target-interactions.R

# Plot significant associations in DMG by immune cell process
Rscript --vanilla 02-dmg-cluster6-mirna-target-dotplot.R

# Identify significant miRNA - immune target assocations among MB cluster1 miRNAs
Rscript --vanilla 03-mb-cluster1-mirna-target-interactions.R

# Plot significant associations in MB by immune cell process
Rscript --vanilla 04-mb-cluster1-mirna-target-dotplot.R

# Correlate DMG cluster 6 miRNA and immune-pathway target expression
Rscript --vanilla 05-dmg-cluster6-mirna-target-pathway-correlations.R

# Correlate MB cluster 1 miRNA and immune-pathway target expression
Rscript --vanilla 06-mb-cluster1-mirna-target-pathway-correlations.R
