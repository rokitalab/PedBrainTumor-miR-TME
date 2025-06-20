# Human miRNA expression analysis

## Usage

`bash run_module.sh`

## Folder contents

1. `01-merge-miRNA-counts.Rmd`: merges two batches of miRNA expression count data, standardizes novel miRNA IDs using mature and precurosr sequence information, and saves the merged output for downstream analysis.
2. `02-mirna-pca-umap.Rmd`: generates PCA (batch comparison, histology, tumor-only) and UMAP plots using top 500 most variable miRNAs across samples.
3. `04-mirna-differential-expression-analysis.Rmd`: performs differential expression analysis of miRNA-seq data using DESeq2.

## Analysis module directory structure

```
.
├── 01-merge-miRNA-counts.Rmd
├── 01-merge-miRNA-counts.html
├── 02-mirna-pca-umap.Rmd
├── 02-mirna-pca-umap.html
├── 04-mirna-differential-expression-analysis.Rmd
├── 04-mirna-differential-expression-analysis.html
├── README.md
├── results
│   ├── 30-931106737-all_novel_miRNA_merged_id.tsv
│   ├── 30-1075661268-all_novel_miRNA_merged_id.tsv
│   ├── merged-miRNA-expression.tsv
│   ├── DE_summary_counts.csv
│   ├── DESeq2_ATRT_vs_adjNormal.csv
│   ├── DESeq2_ATRT_vs_adjNormal_sig.csv
│   ├── DESeq2_ATRT_vs_healthyNormal.csv
│   ├── DESeq2_ATRT_vs_healthyNormal_sig.csv
│   ├── DESeq2_DIPG_vs_adjNormal.csv
│   ├── DESeq2_DIPG_vs_adjNormal_sig.csv
│   ├── DESeq2_DIPG_vs_healthyNormal.csv
│   ├── DESeq2_DIPG_vs_healthyNormal_sig.csv
│   ├── DESeq2_Ependymoma_vs_adjNormal.csv
│   ├── DESeq2_Ependymoma_vs_adjNormal_sig.csv
│   ├── DESeq2_Ependymoma_vs_healthyNormal.csv
│   ├── DESeq2_Ependymoma_vs_healthyNormal_sig.csv
│   ├── DESeq2_HGG_vs_adjNormal.csv
│   ├── DESeq2_HGG_vs_adjNormal_sig.csv
│   ├── DESeq2_HGG_vs_healthyNormal.csv
│   ├── DESeq2_HGG_vs_healthyNormal_sig.csv
│   ├── DESeq2_LGG_vs_healthyNormal.csv
│   ├── DESeq2_LGG_vs_healthyNormal_sig.csv
│   ├── DESeq2_Medulloblastoma_vs_adjNormal.csv
│   ├── DESeq2_Medulloblastoma_vs_adjNormal_sig.csv
│   ├── DESeq2_Medulloblastoma_vs_healthyNormal.csv
│   └── DESeq2_Medulloblastoma_vs_healthyNormal_sig.csv
├── plots
│   ├── mirna-pca-uncorrected.pdf
│   ├── mirna-pca-corrected.pdf
│   ├── mirna-pca-uncorrected-filtered.pdf
│   ├── mirna-pca-uncorrected-filtered-outlier-removed.pdf
│   ├── mirna-pca-tumor-only-corrected.pdf
│   ├── mirna-pca-tumor-only-uncorrected.pdf
│   ├── mirna-pca-ATRT-vs-controls.pdf
│   ├── mirna-pca-DIPG or DMG-vs-controls.pdf
│   ├── mirna-pca-EPN-vs-controls.pdf
│   ├── mirna-pca-GNT-vs-controls.pdf
│   ├── mirna-pca-HGG-vs-controls.pdf
│   ├── mirna-pca-LGG-vs-controls.pdf
│   ├── mirna-pca-MB-vs-controls.pdf
│   ├── mirna-umap-condition.pdf
│   ├── mirna-umap-sample-type.pdf
│   ├── volcano_ATRT_vs_adjNormal.pdf
│   ├── volcano_ATRT_vs_healthyNormal.pdf
│   ├── volcano_DIPG_vs_adjNormal.pdf
│   ├── volcano_DIPG_vs_healthyNormal.pdf
│   ├── volcano_Ependymoma_vs_adjNormal.pdf
│   ├── volcano_Ependymoma_vs_healthyNormal.pdf
│   ├── volcano_HGG_vs_adjNormal.pdf
│   ├── volcano_HGG_vs_healthyNormal.pdf
│   ├── volcano_LGG_vs_healthyNormal.pdf
│   ├── volcano_Medulloblastoma_vs_adjNormal.pdf
│   └── volcano_Medulloblastoma_vs_healthyNormal.pdf
└── run_module.sh
```