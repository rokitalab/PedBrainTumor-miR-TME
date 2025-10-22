#!/bin/bash
# -------------------------------------------------------------------------
# Script: run_module.sh
# Author: Bicna Song
# Date: 2025-10
# Description:
#   Run miranda for mouse DE miRNAs target prediction.
# -------------------------------------------------------------------------

set -euo pipefail

# Run miranda
Rscript -e "rmarkdown::render('01-run-miranda-mouse.Rmd')"
echo "[$(date)] miranda analysis completed successfully."

