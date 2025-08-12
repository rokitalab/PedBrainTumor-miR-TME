# Human RNA expression analysis

## Usage

`bash run_module.sh`

## Folder contents

1. `01-rna-pca-umap.Rmd`: generates PCA (batch comparison, histology, tumor-only) and UMAP plots using top 500 most variable RNAs across samples.
2. `02-differential-expression.Rmd`: Run DESeq2 differential expression (DE) analyses on tumors versus matched normal and healthy normals.
3. `03-GO-enrichment.Rmd`: Run GO enrichment using TopGO on tumor DEGs.
4. `04-GSEA-hallmark.Rmd`: Run GSEA on RNA-seq DE results across multiple histologies to identify enriched Hallmark pathways.
5. `05-GSVA-GOBP.Rmd`: Perform Gene Set Variation Analysis (GSVA) on RSEM-derived TPM expression for multiple histologies, using GO Biological Process terms, and assess differential pathway activity.
6. `06-miRNA-GSVA-correlation.Rmd`: Correlate DE miRNA expression with GSVA pathway scores in DIPG/DMG tumors; generate full and significant-pairs heatmaps with stable annotations.

## Analysis module directory structure

```
.
├── 01-rna-pca.Rmd
├── 01-rna-pca.html
├── 02-rna-differential-expression-analysis.Rmd
├── 02-rna-differential-expression-analysis.html
├── 03-GO-enrichment.Rmd
├── 03-GO-enrichment.html
├── 04-GSEA-hallmark.Rmd
├── 04-GSEA-hallmark.html
├── 05-GSVA-GOBP.Rmd
├── 05-GSVA-GOBP.html
├── 06-miRNA-GSVA-correlation.Rmd
├── 06-miRNA-GSVA-correlation.html
├── README.md
├── input
│   ├── DIPG or DMG_sig_DE_miRNA_list.csv
│   └── mirna-tpm.rds
├── plots
│   ├── DIPG or DMG-downregulated-genes-healthyNormal-go-term-enrichment-dotplot.pdf
│   ├── DIPG or DMG-downregulated-genes-paired-go-term-enrichment-dotplot.pdf
│   ├── DIPG or DMG-upregulated-genes-healthyNormal-go-term-enrichment-dotplot.pdf
│   ├── DIPG or DMG-upregulated-genes-paired-go-term-enrichment-dotplot.pdf
│   ├── DIPG or DMG_healthyNormal_GSVA_DE_heatmap.pdf
│   ├── DIPG or DMG_healthyNormal_HALLMARK_ANGIOGENESIS_GSEA.pdf
│   ├── DIPG or DMG_healthyNormal_HALLMARK_E2F_TARGETS_GSEA.pdf
│   ├── DIPG or DMG_healthyNormal_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_GSEA.pdf
│   ├── DIPG or DMG_healthyNormal_HALLMARK_G2M_CHECKPOINT_GSEA.pdf
│   ├── DIPG or DMG_healthyNormal_HALLMARK_MITOTIC_SPINDLE_GSEA.pdf
│   ├── DIPG or DMG_paired_GSVA_DE_heatmap.pdf
│   ├── DIPG or DMG_paired_HALLMARK_E2F_TARGETS_GSEA.pdf
│   ├── DIPG or DMG_paired_HALLMARK_G2M_CHECKPOINT_GSEA.pdf
│   ├── DIPG or DMG_paired_HALLMARK_INTERFERON_ALPHA_RESPONSE_GSEA.pdf
│   ├── DIPG or DMG_paired_HALLMARK_INTERFERON_GAMMA_RESPONSE_GSEA.pdf
│   ├── DIPG or DMG_paired_HALLMARK_MYC_TARGETS_V1_GSEA.pdf
│   ├── EPN-downregulated-genes-healthyNormal-go-term-enrichment-dotplot.pdf
│   ├── EPN-downregulated-genes-paired-go-term-enrichment-dotplot.pdf
│   ├── EPN-upregulated-genes-healthyNormal-go-term-enrichment-dotplot.pdf
│   ├── EPN-upregulated-genes-paired-go-term-enrichment-dotplot.pdf
│   ├── EPN_healthyNormal_HALLMARK_E2F_TARGETS_GSEA.pdf
│   ├── EPN_healthyNormal_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_GSEA.pdf
│   ├── EPN_healthyNormal_HALLMARK_G2M_CHECKPOINT_GSEA.pdf
│   ├── EPN_healthyNormal_HALLMARK_HYPOXIA_GSEA.pdf
│   ├── EPN_healthyNormal_HALLMARK_UNFOLDED_PROTEIN_RESPONSE_GSEA.pdf
│   ├── EPN_paired_HALLMARK_ANGIOGENESIS_GSEA.pdf
│   ├── EPN_paired_HALLMARK_E2F_TARGETS_GSEA.pdf
│   ├── EPN_paired_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_GSEA.pdf
│   ├── EPN_paired_HALLMARK_G2M_CHECKPOINT_GSEA.pdf
│   ├── EPN_paired_HALLMARK_TNFA_SIGNALING_VIA_NFKB_GSEA.pdf
│   ├── LGG-downregulated-genes-healthyNormal-go-term-enrichment-dotplot.pdf
│   ├── LGG-upregulated-genes-healthyNormal-go-term-enrichment-dotplot.pdf
│   ├── LGG_healthyNormal_HALLMARK_COAGULATION_GSEA.pdf
│   ├── LGG_healthyNormal_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_GSEA.pdf
│   ├── LGG_healthyNormal_HALLMARK_INTERFERON_GAMMA_RESPONSE_GSEA.pdf
│   ├── LGG_healthyNormal_HALLMARK_TGF_BETA_SIGNALING_GSEA.pdf
│   ├── LGG_healthyNormal_HALLMARK_UV_RESPONSE_DN_GSEA.pdf
│   ├── MB-downregulated-genes-healthyNormal-go-term-enrichment-dotplot.pdf
│   ├── MB-downregulated-genes-paired-go-term-enrichment-dotplot.pdf
│   ├── MB-upregulated-genes-healthyNormal-go-term-enrichment-dotplot.pdf
│   ├── MB-upregulated-genes-paired-go-term-enrichment-dotplot.pdf
│   ├── MB_healthyNormal_GSVA_DE_heatmap.pdf
│   ├── MB_healthyNormal_HALLMARK_E2F_TARGETS_GSEA.pdf
│   ├── MB_healthyNormal_HALLMARK_G2M_CHECKPOINT_GSEA.pdf
│   ├── MB_healthyNormal_HALLMARK_MITOTIC_SPINDLE_GSEA.pdf
│   ├── MB_healthyNormal_HALLMARK_MYC_TARGETS_V1_GSEA.pdf
│   ├── MB_healthyNormal_HALLMARK_MYC_TARGETS_V2_GSEA.pdf
│   ├── MB_paired_GSVA_DE_heatmap.pdf
│   ├── MB_paired_HALLMARK_DNA_REPAIR_GSEA.pdf
│   ├── MB_paired_HALLMARK_E2F_TARGETS_GSEA.pdf
│   ├── MB_paired_HALLMARK_G2M_CHECKPOINT_GSEA.pdf
│   ├── MB_paired_HALLMARK_MTORC1_SIGNALING_GSEA.pdf
│   ├── MB_paired_HALLMARK_MYC_TARGETS_V1_GSEA.pdf
│   ├── heatmap_all.pdf
│   ├── heatmap_sig_cluster.pdf
│   ├── rna-pca-DIPG or DMG-vs-controls.pdf
│   ├── rna-pca-EPN-vs-controls.pdf
│   ├── rna-pca-LGG-vs-controls.pdf
│   ├── rna-pca-MB-vs-controls.pdf
│   ├── rna-pca.pdf
│   ├── silhouette_curves_ht2.pdf
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
│   ├── DIPG or DMG-downregulated-genes-healthyNormal-enriched-go-terms.tsv
│   ├── DIPG or DMG-downregulated-genes-paired-enriched-go-terms.tsv
│   ├── DIPG or DMG-upregulated-genes-healthyNormal-enriched-go-terms.tsv
│   ├── DIPG or DMG-upregulated-genes-paired-enriched-go-terms.tsv
│   ├── DIPG or DMG_healthyNormal_GSEA_results.tsv
│   ├── DIPG or DMG_healthyNormal_GSVA_DE_results.tsv
│   ├── DIPG or DMG_paired_GSEA_results.tsv
│   ├── DIPG or DMG_paired_GSVA_DE_results.tsv
│   ├── EPN-downregulated-genes-healthyNormal-enriched-go-terms.tsv
│   ├── EPN-downregulated-genes-paired-enriched-go-terms.tsv
│   ├── EPN-upregulated-genes-healthyNormal-enriched-go-terms.tsv
│   ├── EPN-upregulated-genes-paired-enriched-go-terms.tsv
│   ├── EPN_healthyNormal_GSEA_results.tsv
│   ├── EPN_paired_GSEA_results.tsv
│   ├── LGG-downregulated-genes-healthyNormal-enriched-go-terms.tsv
│   ├── LGG-upregulated-genes-healthyNormal-enriched-go-terms.tsv
│   ├── LGG_healthyNormal_GSEA_results.tsv
│   ├── MB-downregulated-genes-healthyNormal-enriched-go-terms.tsv
│   ├── MB-downregulated-genes-paired-enriched-go-terms.tsv
│   ├── MB-upregulated-genes-healthyNormal-enriched-go-terms.tsv
│   ├── MB-upregulated-genes-paired-enriched-go-terms.tsv
│   ├── MB_healthyNormal_GSEA_results.tsv
│   ├── MB_healthyNormal_GSVA_DE_results.tsv
│   ├── MB_paired_GSEA_results.tsv
│   ├── MB_paired_GSVA_DE_results.tsv
│   ├── gobp-gsva-scores-with-target-genes.tsv
│   └── gobp-gsva-scores.tsv
└── run_module.sh
```