# Human miRNA expression analysis

## Usage

`bash run_module.sh`

## Folder contents

2. `02-differential-expression.Rmd`: Run DESeq2 differential expression analyses on tumors versus matched normal and healthy normals. 

## Analysis module directory structure

```
.
├── 02-differential-expression.Rmd
├── plots
│   ├── volcano_ATRT_vs_adjNormal.pdf
│   ├── volcano_ATRT_vs_healthyNormal.pdf
│   ├── volcano_DIPG or DMG_paired.pdf
│   ├── volcano_DIPG or DMG_vs_adjNormal.pdf
│   ├── volcano_DIPG or DMG_vs_healthyNormal.pdf
│   ├── volcano_EPN_paired.pdf
│   ├── volcano_EPN_vs_adjNormal.pdf
│   ├── volcano_EPN_vs_healthyNormal.pdf
│   ├── volcano_HGG_paired.pdf
│   ├── volcano_HGG_vs_adjNormal.pdf
│   ├── volcano_HGG_vs_healthyNormal.pdf
│   ├── volcano_LGG_vs_healthyNormal.pdf
│   ├── volcano_MB_paired.pdf
│   ├── volcano_MB_vs_adjNormal.pdf
│   └── volcano_MB_vs_healthyNormal.pdf
├── results
│   ├── DE_summary_counts_paired.csv
│   ├── DE_summary_counts.csv
│   ├── DESeq2_ATRT_vs_adjNormal_sig.csv
│   ├── DESeq2_ATRT_vs_adjNormal.csv
│   ├── DESeq2_ATRT_vs_healthyNormal_sig.csv
│   ├── DESeq2_ATRT_vs_healthyNormal.csv
│   ├── DESeq2_DIPG or DMG_paired_full.csv
│   ├── DESeq2_DIPG or DMG_paired_sig.csv
│   ├── DESeq2_DIPG or DMG_vs_adjNormal_sig.csv
│   ├── DESeq2_DIPG or DMG_vs_adjNormal.csv
│   ├── DESeq2_DIPG or DMG_vs_healthyNormal_sig.csv
│   ├── DESeq2_DIPG or DMG_vs_healthyNormal.csv
│   ├── DESeq2_EPN_paired_full.csv
│   ├── DESeq2_EPN_paired_sig.csv
│   ├── DESeq2_EPN_vs_adjNormal_sig.csv
│   ├── DESeq2_EPN_vs_adjNormal.csv
│   ├── DESeq2_EPN_vs_healthyNormal_sig.csv
│   ├── DESeq2_EPN_vs_healthyNormal.csv
│   ├── DESeq2_HGG_paired_full.csv
│   ├── DESeq2_HGG_paired_sig.csv
│   ├── DESeq2_HGG_vs_adjNormal_sig.csv
│   ├── DESeq2_HGG_vs_adjNormal.csv
│   ├── DESeq2_HGG_vs_healthyNormal_sig.csv
│   ├── DESeq2_HGG_vs_healthyNormal.csv
│   ├── DESeq2_LGG_vs_healthyNormal_sig.csv
│   ├── DESeq2_LGG_vs_healthyNormal.csv
│   ├── DESeq2_MB_paired_full.csv
│   ├── DESeq2_MB_paired_sig.csv
│   ├── DESeq2_MB_vs_adjNormal_sig.csv
│   ├── DESeq2_MB_vs_adjNormal.csv
│   ├── DESeq2_MB_vs_healthyNormal_sig.csv
│   ├── DESeq2_MB_vs_healthyNormal.csv
└── run_module.sh
```