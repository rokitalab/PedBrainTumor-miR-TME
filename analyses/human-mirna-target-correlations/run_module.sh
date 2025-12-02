#!/bin/bash
# Module author: Ryan Corbett
# 2025-10

set -e
set -o pipefail

# Identify significant miRNA - immune target assocations among DMG cluster6 miRNAs
Rscript --vanilla 01-dmg-cluster6-mirna-target-interactions.R

# Plot significant associations by immune cell process
Rscript --vanilla 02-dmg-cluster6-mirna-target-dotplot.R