# Human mouse miRNA orthology analyses 

This module identifies mouse orthologs to human miRNAs that are DE in DMG versus normal brain, and assesses their expression patterns in B7H3 treated versus untreated mice. 

## Usage

`bash run_module.sh`

## Folder contents

1. `01-summarize-orthology-expression.Rmd`: identifies mouse orthologs and assessses DMG cluster-specific expression patterns in B7H3 treated mice
2. `02-plot-cluster6-ortholog-expr.R`: Plots expression of human DMG cluster6 orthologs that are DE in B7H3 treated versus untreated mice. 

## Analysis module directory structure

```
.
├── 01-summarize-orthology-expression.Rmd
├── 02-plot-cluster6-ortholog-expr.R
├── plots
│   ├── hsa-dmg-clust6-mouse-ortholog-de-tpm.pdf
│   ├── hsa-dmg-de-mirna-mmu-ortho-barplot.pdf
│   └── human-dmg-de-mirna-orthology-mouse-expr-enrichment-plot.pdf
├── README.md
└── run_module.sh
```