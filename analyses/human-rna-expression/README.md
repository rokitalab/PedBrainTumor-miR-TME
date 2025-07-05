# Human RNA expression analysis

## Usage

`bash run_module.sh`

## Folder contents

2. `02-differential-expression.Rmd`: Run DESeq2 differential expression (DE) analyses on tumors versus matched normal and healthy normals.
3. `03-GO-enrichment.Rmd`: Run GO enrichment using TopGO on tumor DEGs.
4. `04-GSEA-hallmark.Rmd`: Run GSEA on RNA-seq DE results across multiple histologies to identify enriched Hallmark pathways.

## Analysis module directory structure

```
.
├── 02-differential-expression.Rmd
├── 02-differential-expression.html
├── 03-GO-enrichment.Rmd
├── 03-GO-enrichment.html
├── 04-GSEA-hallmark.Rmd
├── 04-GSEA-hallmark.html
├── README.md
├── plots
│   ├── DIPG or DMG-downregulated-genes-healthyNormal-go-term-enrichment-dotplot.pdf
│   ├── DIPG or DMG-downregulated-genes-paired-go-term-enrichment-dotplot.pdf
│   ├── DIPG or DMG-upregulated-genes-healthyNormal-go-term-enrichment-dotplot.pdf
│   ├── DIPG or DMG-upregulated-genes-paired-go-term-enrichment-dotplot.pdf
│   ├── DIPG or DMG_healthyNormal_HALLMARK_ANGIOGENESIS_GSEA.pdf
│   ├── DIPG or DMG_healthyNormal_HALLMARK_E2F_TARGETS_GSEA.pdf
│   ├── DIPG or DMG_healthyNormal_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_GSEA.pdf
│   ├── DIPG or DMG_healthyNormal_HALLMARK_G2M_CHECKPOINT_GSEA.pdf
│   ├── DIPG or DMG_healthyNormal_HALLMARK_MITOTIC_SPINDLE_GSEA.pdf
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
│   ├── HGG_healthyNormal_HALLMARK_ALLOGRAFT_REJECTION_GSEA.pdf
│   ├── HGG_healthyNormal_HALLMARK_E2F_TARGETS_GSEA.pdf
│   ├── HGG_healthyNormal_HALLMARK_G2M_CHECKPOINT_GSEA.pdf
│   ├── HGG_healthyNormal_HALLMARK_IL6_JAK_STAT3_SIGNALING_GSEA.pdf
│   ├── HGG_healthyNormal_HALLMARK_MITOTIC_SPINDLE_GSEA.pdf
│   ├── LGG-downregulated-genes-healthyNormal-go-term-enrichment-dotplot.pdf
│   ├── LGG-upregulated-genes-healthyNormal-go-term-enrichment-dotplot.pdf
│   ├── LGG_healthyNormal_HALLMARK_COAGULATION_GSEA.pdf
│   ├── LGG_healthyNormal_HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION_GSEA.pdf
│   ├── LGG_healthyNormal_HALLMARK_INTERFERON_GAMMA_RESPONSE_GSEA.pdf
│   ├── LGG_healthyNormal_HALLMARK_TGF_BETA_SIGNALING_GSEA.pdf
│   ├── LGG_healthyNormal_HALLMARK_UV_RESPONSE_DN_GSEA.pdf
│   ├── MB-downregulated-genes-paired-go-term-enrichment-dotplot.pdf
│   ├── MB-upregulated-genes-healthyNormal-go-term-enrichment-dotplot.pdf
│   ├── MB-upregulated-genes-paired-go-term-enrichment-dotplot.pdf
│   ├── MB_healthyNormal_HALLMARK_E2F_TARGETS_GSEA.pdf
│   ├── MB_healthyNormal_HALLMARK_G2M_CHECKPOINT_GSEA.pdf
│   ├── MB_healthyNormal_HALLMARK_MITOTIC_SPINDLE_GSEA.pdf
│   ├── MB_healthyNormal_HALLMARK_MYC_TARGETS_V1_GSEA.pdf
│   ├── MB_healthyNormal_HALLMARK_MYC_TARGETS_V2_GSEA.pdf
│   ├── MB_paired_HALLMARK_E2F_TARGETS_GSEA.pdf
│   ├── MB_paired_HALLMARK_G2M_CHECKPOINT_GSEA.pdf
│   ├── MB_paired_HALLMARK_MITOTIC_SPINDLE_GSEA.pdf
│   ├── MB_paired_HALLMARK_MTORC1_SIGNALING_GSEA.pdf
│   ├── MB_paired_HALLMARK_MYC_TARGETS_V1_GSEA.pdf
│   ├── volcano_ATRT_vs_adjNormal.pdf
│   ├── volcano_ATRT_vs_healthyNormal.pdf
│   ├── volcano_DIPG or DMG_paired.pdf
│   ├── volcano_DIPG or DMG_vs_adjNormal.pdf
│   ├── volcano_DIPG or DMG_vs_healthyNormal.pdf
│   ├── volcano_EPN_paired.pdf
│   ├── volcano_EPN_vs_adjNormal.pdf
│   ├── volcano_EPN_vs_healthyNormal.pdf
│   ├── volcano_HGG_paired.pdf
│   ├── volcano_HGG_vs_adjNormal.pdf
│   ├── volcano_HGG_vs_healthyNormal.pdf
│   ├── volcano_LGG_vs_healthyNormal.pdf
│   ├── volcano_MB_paired.pdf
│   ├── volcano_MB_vs_adjNormal.pdf
│   └── volcano_MB_vs_healthyNormal.pdf
├── results
│   ├── ATRT_healthyNormal_GSEA_results.tsv
│   ├── DESeq2_ATRT_vs_adjNormal.csv
│   ├── DESeq2_ATRT_vs_adjNormal_sig.csv
│   ├── DESeq2_ATRT_vs_healthyNormal.csv
│   ├── DESeq2_ATRT_vs_healthyNormal_sig.csv
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
│   ├── DESeq2_HGG_paired_full.csv
│   ├── DESeq2_HGG_paired_sig.csv
│   ├── DESeq2_HGG_vs_adjNormal.csv
│   ├── DESeq2_HGG_vs_adjNormal_sig.csv
│   ├── DESeq2_HGG_vs_healthyNormal.csv
│   ├── DESeq2_HGG_vs_healthyNormal_sig.csv
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
│   ├── DIPG or DMG_paired_GSEA_results.tsv
│   ├── EPN-downregulated-genes-healthyNormal-enriched-go-terms.tsv
│   ├── EPN-downregulated-genes-paired-enriched-go-terms.tsv
│   ├── EPN-upregulated-genes-healthyNormal-enriched-go-terms.tsv
│   ├── EPN-upregulated-genes-paired-enriched-go-terms.tsv
│   ├── EPN_healthyNormal_GSEA_results.tsv
│   ├── EPN_paired_GSEA_results.tsv
│   ├── HGG-downregulated-genes-healthyNormal-enriched-go-terms.tsv
│   ├── HGG-upregulated-genes-healthyNormal-enriched-go-terms.tsv
│   ├── HGG-upregulated-genes-paired-enriched-go-terms.tsv
│   ├── HGG_healthyNormal_GSEA_results.tsv
│   ├── LGG-downregulated-genes-healthyNormal-enriched-go-terms.tsv
│   ├── LGG-upregulated-genes-healthyNormal-enriched-go-terms.tsv
│   ├── LGG_healthyNormal_GSEA_results.tsv
│   ├── MB-downregulated-genes-healthyNormal-enriched-go-terms.tsv
│   ├── MB-downregulated-genes-paired-enriched-go-terms.tsv
│   ├── MB-upregulated-genes-healthyNormal-enriched-go-terms.tsv
│   ├── MB-upregulated-genes-paired-enriched-go-terms.tsv
│   ├── MB_healthyNormal_GSEA_results.tsv
│   └── MB_paired_GSEA_results.tsv
└── run_module.sh
```