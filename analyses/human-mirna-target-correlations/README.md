# Human DE miRNA - target associations

This module identifies significant associations between DE miRNAs in clusters of interest and predicted immune-related targets as follows:

1) Calculates pearson correlation coefficients between DE miRNA and target expression
2) Assesses whether predicted targets exhibit differential expression in the opposite direction of predicted targeting miRNAs

## Usage

`bash run_module.sh`

## Folder contents

1. `01-dmg-cluster6-mirna-target-interactions.R`: identify significant associations between DMG cluster6 miRNAs and predicted immune-related targets
2. `02-dmg-cluster6-mirna-target-dotplot.R`: Plot significant associations between cluster6 miRNAs and predicted immune-related targets

## Analysis module directory structure

```
.
├── 01-dmg-cluster6-mirna-target-interactions.R
├── 02-dmg-cluster6-mirna-target-dotplot.R
├── dmg-cluster6-mirna-target-dotplot.R
├── input
│   └── dipg-dmg-cluster6-mirna-immune-target-sig-interactions.txt
├── plots
│   ├── dmg-cluster6-target-myeloid-cell-dotplot.pdf
│   ├── dmg-cluster6-target-other-term-dotplot.pdf
│   └── dmg-cluster6-target-t-cell-dotplot.pdf
├── README.md
├── results
│   ├── dmg-cluster6-mirna-immune-target-interactions.tsv
│   └── dmg-cluster6-mirna-immune-target-sig-interactions.tsv
└── run_module.sh
```