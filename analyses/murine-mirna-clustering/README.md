# DE miRNA clustering

## Usage

`bash run_module.sh`

## Folder contents

1. `01-murine-de-mirna-clustering.R`: clusters DE miRNAs from CAR-T versus untreated mice by expression

## Analysis module directory structure

```
.
├── 01-murine-de-mirna-clustering.R
├── README.md
├── plots
│   ├── de-mirna-heatmap-Day14.pdf
│   ├── de-mirna-heatmap-Day21.pdf
│   └── de-mirna-heatmap.pdf
├── results
│   ├── mouse_sig_DE_miRNA_list.csv
│   ├── mouse-de-mirna-cluster-membership.tsv
│   └── murine-mirna-tpm.rds
└── run_module.sh
```