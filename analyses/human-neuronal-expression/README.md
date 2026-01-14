# Human neuronal expression analyses

## Usage

`bash run_module.sh`

## Folder contents

1. `01-human-neuronal-expression.Rmd`: performs neuronal gene–level comparisons in tumor vs. normal brain samples across DIPG/DMG, MB, EPN, and LGG cohorts.
2. `02-human-neuronal-mirna-target-correlation.Rmd`: performs correlation analysis between DE miRNAs and GSVA scores of neuronal-related GO terms.
3. `03-human-neuronal-heatmap.R`: regenerates the miRNA clustering heatmap with neuronal GO BP correlation coefficients as annotations.

## Analysis module directory structure

```
.
├── 01-human-neuronal-expression.Rmd
├── 01-human-neuronal-expression.html
├── 02-human-neuronal-mirna-target-correlation.Rmd
├── 02-human-neuronal-mirna-target-correlation.html
├── 03-human-neuronal-heatmap.R
├── 03-human-neuronal-heatmap.html
├── README.md
├── plots
│   ├── DIPG_or_DMG_mirna_neuronal_gobp_heatmap.pdf
│   ├── neuronal_DIPG_or_DMG_healthyNormal.pdf
│   ├── neuronal_DIPG_or_DMG_paired.pdf
│   ├── neuronal_EPN_healthyNormal.pdf
│   ├── neuronal_EPN_paired.pdf
│   ├── neuronal_LGG_healthyNormal.pdf
│   ├── neuronal_MB_healthyNormal.pdf
│   └── neuronal_MB_paired.pdf
├── results
│   ├── DIPG_or_DMG_mirna_neuronal_cluster_memebership.tsv
│   ├── mirna_neuronal_gobp_spearman.tsv
│   ├── mirna_neuronal_gobp_spearman_FDR01.tsv
│   └── neuronal_gobp_terms_from_gsva.tsv
└── run_module.sh
```

