# Immune deconvolution analyses

## Usage

`bash run_module.sh`

## Folder contents

1. `01-run-xcell2.R`: Performs immune deconvolution on mouse bulk RNA-seq data using the xCell2 framework.
2. `02-xcell2-heatmap.R`: Generates scaled xCell2 enrichment heatmaps for each reference dataset.

## Analysis module directory structure

```
.
├── 01-run-xcell2.R
├── 02-xcell2-heatmap.R
├── README.md
├── plots
│   ├── xCell2_scaled_heatmap_ImmGenData.pdf
│   ├── xCell2_scaled_heatmap_MouseRNAseqData.pdf
│   └── xCell2_scaled_heatmap_TabulaMurisBlood.pdf
├── results
│   ├── ImmGenData_xCell2_results.tsv
│   ├── MouseRNAseqData_xCell2_results.tsv
│   ├── TabulaMurisBlood_xCell2_results.tsv
│   ├── xCell2_all_references_merged.tsv
│   └── xCell2_scaled_within_reference.tsv
└── run_module.sh
```