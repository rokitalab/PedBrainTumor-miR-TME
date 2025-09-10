# DE miRNA clustering

This module performs hierarchical clustering of differentially expressed miRNAs in DIPG/DMG and medulloblastoma, and appends immune cell fraction correlations to heatmaps

## Usage

`bash run_module.sh`

## Folder contents

1. `01-de-mirna-clustering.R`: perform clustering and generate heatmap with immune cell fraction correlation annotations

## Analysis module directory structure

```
.
├── 01-de-mirna-clustering.R
├── plots
│   ├── DIPG or DMG-de-mirna-heatmap.pdf
│   └── MB-de-mirna-heatmap.pdf
├── README.md
├── results
│   ├── DIPG or DMG-de-mirna-cluster-membership-immune-scores.tsv
│   └── MB-de-mirna-cluster-membership-immune-scores.tsv
└── run_module.sh
```