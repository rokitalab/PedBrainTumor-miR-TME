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
Rscript --vanilla 01-human-neuronal-expression.Rmd
echo "[$(date)] Neuronal expression analysis completed successfully."

