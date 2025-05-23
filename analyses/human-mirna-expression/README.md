# Human miRNA expression analysis

## Usage

`bash run_module.sh`

## Folder contents

1. `01-merge-miRNA-counts.Rmd`: merges two batches of miRNA expression count data, standardizes novel miRNA IDs using mature and precurosr sequence information, and saves the merged output for downstream analysis.
2. `02-pca-mirna-batch-comparison.Rmd`: compares PCA plots of miRNA expression data before and after batch correction.
3. `03-pca-by-histology.Rmd`: generates PCA plots of miRNA expression data to evaluate batch effects and biological variation across samples.

## Analysis module directory structure

```
.
├── 01-merge-miRNA-counts.Rmd
├── 01-merge-miRNA-counts.html
├── 02-pca-mirna-batch-comparison.Rmd
├── 02-pca-mirna-batch-comparison.html
├── 03-pca-by-histology.Rmd
├── 03-pca-by-histology.html
├── README.md
├── results
│   ├── 30-931106737-all_novel_miRNA_merged_id.tsv
│   ├── 30-1075661268-all_novel_miRNA_merged_id.tsv
│   └── merged-miRNA-expression.tsv
├── plots
│   ├── mirna-pca-corrected.pdf
│   ├── mirna-pca-uncorrected_filtered_n50.pdf
│   ├── mirna-pca-uncorrected.pdf
│   ├── mirna-pca-tumor-only-uncorrected.pdf
│   ├── mirna-pca-tumor-only-corrected.pdf
│   ├── mirna-pca-tumor-only-corrected-combat.pdf
│   ├── pca_ATRT_vs_controls.pdf
│   ├── pca_DIPG_vs_controls.pdf
│   ├── pca_Ependymoma_vs_controls.pdf
│   ├── pca_HGG_vs_controls.pdf
│   ├── pca_LGG_vs_controls.pdf
│   └── pca_Medulloblastoma_vs_controls.pdf
└── run_module.sh
```