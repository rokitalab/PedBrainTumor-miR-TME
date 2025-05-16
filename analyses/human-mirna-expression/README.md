# Human miRNA expression analysis

## Usage

`bash run_module.sh`

## Folder contents

1. `01-merge-miRNA-counts.Rmd`: merge merges two batches of miRNA expression count data, standardizes novel miRNA IDs using mature and precurosr sequence information, and saves the merged output for downstream analysis.
2. `02-pca-mirna-batch-comparison.Rmd`: compares PCA plots of miRNA expression data before and after batch correction.

## Analysis module directory structure

```
.
├── 01-merge-miRNA-counts.Rmd
├── 01-merge-miRNA-counts.html
├── 02-pca-mirna-batch-comparison.Rmd
├── 02-pca-mirna-batch-comparison.html
├── README.md
├── results
│   ├── 30-931106737-all_novel_miRNA_merged_id.tsv
│   ├── 30-1075661268-all_novel_miRNA_merged_id.tsv
│   └── merged-miRNA-expression.tsv
├── plots
│   ├── mirna-pca-corrected.pdf
│   └── mirna-pca-uncorrected.pdf
└── run_module.sh
```