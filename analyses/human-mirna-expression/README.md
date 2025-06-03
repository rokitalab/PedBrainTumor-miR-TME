# Human miRNA expression analysis

## Usage

`bash run_module.sh`

## Folder contents

1. `01-merge-miRNA-counts.Rmd`: merges two batches of miRNA expression count data, standardizes novel miRNA IDs using mature and precurosr sequence information, and saves the merged output for downstream analysis.
2. `02-pca-mirna-batch-comparison.Rmd`: compares PCA plots of miRNA expression data before and after batch correction.
3. `03-pca-by-histology.Rmd`: generates PCA plots of miRNA expression data to evaluate batch effects and biological variation across samples.
4. `04-mirna-differential-expression-analysis.Rmd`: performs differential expression analysis of miRNA-seq data using DESeq2.
5. `05-umap-mirna.Rmd`: investigates whether miRNA-seq samples show batch effects or biological clustering based on histology and sample type. The goal is to visualize global expression patterns using UMAP.

## Analysis module directory structure

```
.
├── 01-merge-miRNA-counts.Rmd
├── 01-merge-miRNA-counts.html
├── 02-pca-mirna-batch-comparison.Rmd
├── 02-pca-mirna-batch-comparison.html
├── 03-pca-by-histology.Rmd
├── 03-pca-by-histology.html
├── 04-mirna-differential-expression-analysis.Rmd
├── 04-mirna-differential-expression-analysis.html
├── 05-umap-mirna.Rmd
├── 05-umap-mirna.html
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
│   ├── mirna-pca-uncorrected_filtered_n50.pdf
│   ├── mirna-pca-corrected.pdf
│   ├── mirna-pca-tumor-only-uncorrected.pdf
│   ├── mirna-pca-tumor-only-corrected-PC1-PC2.pdf
│   ├── mirna-pca-tumor-only-corrected-PC2-PC3.pdf
│   ├── mirna-pca-tumor-only-corrected-PC3-PC4.pdf
│   ├── mirna-pca-tumor-only-corrected-combat.pdf
│   ├── mirna-pca-tumor-only-corrected-combat.pdf
│   ├── mirna_umap.pdf
│   ├── mirna_umap_sample_type.pdf
│   ├── mirna_umap_filtered.pdf
│   ├── mirna_umap_filtered_sample_type.pdf
│   ├── pca_ATRT_vs_controls.pdf
│   ├── pca_DIPG_vs_controls.pdf
│   ├── pca_Ependymoma_vs_controls.pdf
│   ├── pca_HGG_vs_controls.pdf
│   ├── pca_LGG_vs_controls.pdf
│   ├── pca_Medulloblastoma_vs_controls.pdf
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