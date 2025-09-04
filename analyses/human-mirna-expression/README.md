# Human miRNA expression analysis

## Usage

`bash run_module.sh`

## Folder contents

1. `01-merge-miRNA-counts.Rmd`: merges two batches of miRNA expression count data, standardizes novel miRNA IDs using mature and precurosr sequence information, and saves the merged output for downstream analysis.
2. `02-mirna-pca-umap.Rmd`: generates PCA (batch comparison, histology, tumor-only) and UMAP plots using top 500 most variable miRNAs across samples.
3. `03-mirna-differential-expression-analysis.Rmd`: performs differential expression analysis of miRNA-seq data using DESeq2.
4. `04-mirna-venn-upset.Rmd`: generates Venn diagrams and UpSet plots to visualize shared and unique dysregulated miRNAs across subtypes and comparisons.

## Analysis module directory structure

```
.
├── 01-merge-miRNA-counts.Rmd
├── 01-merge-miRNA-counts.html
├── 02-mirna-pca-umap.Rmd
├── 02-mirna-pca-umap.html
├── 03-mirna-differential-expression-analysis.Rmd
├── 03-mirna-differential-expression-analysis.html
├── 04-mirna-venn-upset.Rmd
├── 04-mirna-venn-upset.html
├── README.md
├── plots
│   ├── mirna-pca-DIPG or DMG-vs-controls.pdf
│   ├── mirna-pca-EPN-vs-controls.pdf
│   ├── mirna-pca-LGG-vs-controls.pdf
│   ├── mirna-pca-MB-vs-controls.pdf
│   ├── mirna-pca-corrected.pdf
│   ├── mirna-pca-tumor-only-corrected.pdf
│   ├── mirna-pca-tumor-only-uncorrected.pdf
│   ├── mirna-pca-uncorrected-filtered-outlier-removed.pdf
│   ├── mirna-pca-uncorrected-filtered.pdf
│   ├── mirna-pca-uncorrected.pdf
│   ├── mirna-umap-condition.pdf
│   ├── mirna-umap-sample-type.pdf
│   ├── upset_downregulated_miRNAs.pdf
│   ├── upset_downregulated_miRNAs_tumor_vs_adjacent.pdf
│   ├── upset_downregulated_miRNAs_tumor_vs_healthynormal.pdf
│   ├── upset_upregulated_miRNAs.pdf
│   ├── upset_upregulated_miRNAs_tumor_vs_adjacent.pdf
│   ├── upset_upregulated_miRNAs_tumor_vs_healthynormal.pdf
│   ├── venn_down_paired_vs_healthy_DIPG or DMG.pdf
│   ├── venn_down_paired_vs_healthy_EPN.pdf
│   ├── venn_down_paired_vs_healthy_MB.pdf
│   ├── venn_up_paired_vs_healthy_DIPG or DMG.pdf
│   ├── venn_up_paired_vs_healthy_EPN.pdf
│   ├── venn_up_paired_vs_healthy_MB.pdf
│   ├── volcano_DIPG or DMG_paired.pdf
│   ├── volcano_DIPG or DMG_vs_adjNormal.pdf
│   ├── volcano_DIPG or DMG_vs_healthyNormal.pdf
│   ├── volcano_EPN_paired.pdf
│   ├── volcano_EPN_vs_adjNormal.pdf
│   ├── volcano_EPN_vs_healthyNormal.pdf
│   ├── volcano_LGG_vs_healthyNormal.pdf
│   ├── volcano_MB_paired.pdf
│   ├── volcano_MB_vs_adjNormal.pdf
│   └── volcano_MB_vs_healthyNormal.pdf
├── results
│   ├── 30-1075661268-all_novel_miRNA_merged_id.tsv
│   ├── 30-931106737-all_novel_miRNA_merged_id.tsv
│   ├── DESeq2_DIPG or DMG_paired_full.csv
│   ├── DESeq2_DIPG or DMG_paired_sig.csv
│   ├── DESeq2_DIPG or DMG_vs_adjNormal.csv
│   ├── DESeq2_DIPG or DMG_vs_adjNormal_sig.csv
│   ├── DESeq2_DIPG or DMG_vs_healthyNormal.csv
│   ├── DESeq2_DIPG or DMG_vs_healthyNormal_sig.csv
│   ├── DESeq2_EPN_paired_full.csv
│   ├── DESeq2_EPN_paired_sig.csv
│   ├── DESeq2_EPN_vs_adjNormal.csv
│   ├── DESeq2_EPN_vs_adjNormal_sig.csv
│   ├── DESeq2_EPN_vs_healthyNormal.csv
│   ├── DESeq2_EPN_vs_healthyNormal_sig.csv
│   ├── DESeq2_LGG_vs_healthyNormal.csv
│   ├── DESeq2_LGG_vs_healthyNormal_sig.csv
│   ├── DESeq2_MB_paired_full.csv
│   ├── DESeq2_MB_paired_sig.csv
│   ├── DESeq2_MB_vs_adjNormal.csv
│   ├── DESeq2_MB_vs_adjNormal_sig.csv
│   ├── DESeq2_MB_vs_healthyNormal.csv
│   ├── DESeq2_MB_vs_healthyNormal_sig.csv
│   ├── DE_summary_counts.csv
│   ├── DE_summary_counts_paired.csv
│   ├── DIPG or DMG_sig_DE_miRNA_list.csv
│   ├── EPN_sig_DE_miRNA_list.csv
│   ├── MB_sig_DE_miRNA_list.csv
│   ├── merged-miRNA-expression.tsv
│   ├── venn_table_down_paired_vs_healthy_DIPG or DMG.csv
│   ├── venn_table_down_paired_vs_healthy_EPN.csv
│   ├── venn_table_down_paired_vs_healthy_MB.csv
│   ├── venn_table_up_paired_vs_healthy_DIPG or DMG.csv
│   ├── venn_table_up_paired_vs_healthy_EPN.csv
│   └── venn_table_up_paired_vs_healthy_MB.csv
└── run_module.sh
```