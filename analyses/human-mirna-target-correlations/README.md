# Human DE miRNA - target associations

This module identifies significant associations between DE miRNAs in clusters of interest and predicted immune-related targets as follows:

1) Calculates pearson correlation coefficients between DE miRNA and target expression
2) Assesses whether predicted targets exhibit differential expression in the opposite direction of predicted targeting miRNAs

## Usage

`bash run_module.sh`

## Folder contents

1. `01-dmg-cluster6-mirna-target-interactions.R`: identify significant associations between DMG cluster6 miRNAs and predicted immune-related targets
2. `02-dmg-cluster6-mirna-target-dotplot.R`: Plot significant associations between DMG cluster6 miRNAs and predicted immune-related targets
3. `03-mb-cluster1-mirna-target-interactions.R`: identify significant associations between MB cluster1 miRNAs and predicted immune-related targets
4. `04-mb-cluster1-mirna-target-dotplot.R`: Plot significant associations between MB cluster1 miRNAs and predicted immune-related targets

### Input files

- `dipg-dmg-cluster6-mirna-immune-target-sig-interactions.txt`: reviewed and annotated target GO terms
- `mb-cluster1-mirna-immune-target-sig-interactions.txt`: reviewed and annotated target GO terms

## Analysis module directory structure

```
.
├── 01-dmg-cluster6-mirna-target-interactions.R
├── 02-dmg-cluster6-mirna-target-dotplot.R
├── 03-mb-cluster1-mirna-target-interactions.R
├── 04-mb-cluster1-mirna-target-dotplot.R
├── input
│   ├── dipg-dmg-cluster6-mirna-immune-target-sig-interactions.txt
│   └── mb-cluster1-mirna-immune-target-sig-interactions.txt
├── plots
│   ├── dmg-cluster6-target-myeloid-cell-dotplot.pdf
│   ├── dmg-cluster6-target-other-term-dotplot.pdf
│   ├── dmg-cluster6-target-t-cell-dotplot.pdf
│   ├── mb-cluster1-target-myeloid-cell-dotplot.pdf
│   ├── mb-cluster1-target-other-term-dotplot.pdf
│   └── mb-cluster1-target-t-cell-dotplot.pdf
├── README.md
├── results
│   ├── dmg-cluster6-mirna-immune-target-interactions.tsv
│   ├── dmg-cluster6-mirna-immune-target-sig-interactions.tsv
│   ├── mb-cluster1-mirna-immune-target-interactions.tsv
│   └── mb-cluster1-mirna-immune-target-sig-interactions.tsv
└── run_module.sh
```