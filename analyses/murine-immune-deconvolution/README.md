# Immune deconvolution analyses

## Usage

`bash run_module.sh`

## Folder contents

1. `01-run-xcell2.R`: Performs immune deconvolution on mouse bulk RNA-seq data using the xCell2 framework.
2. `02-xcell2-heatmap.R`: Generates scaled xCell2 enrichment heatmaps for each reference dataset.
3. `03-xcell2-plot.Rmd`: Generates summary plots of xCell2 immune cell enrichment scores across treatments, time points, and reference datasets.

## Analysis module directory structure

```
.
├── 01-run-xcell2.R
├── 02-xcell2-heatmap.R
├── 03-xcell2-plot.Rmd
├── 03-xcell2-plot.html
├── README.md
├── plots
│   ├── xCell2_scaled_heatmap_ImmGenData.pdf
│   ├── xCell2_scaled_heatmap_MouseRNAseqData.pdf
│   ├── xCell2_scaled_heatmap_TabulaMurisBlood.pdf
│   ├── xcell2_EnrichmentScore_ImmGenData_Day14.pdf
│   ├── xcell2_EnrichmentScore_ImmGenData_Day21.pdf
│   ├── xcell2_EnrichmentScore_MouseRNAseqData_Day14.pdf
│   ├── xcell2_EnrichmentScore_MouseRNAseqData_Day21.pdf
│   ├── xcell2_EnrichmentScore_byTreatment_Day14.pdf
│   └── xcell2_EnrichmentScore_byTreatment_Day21.pdf
├── results
│   ├── ImmGenData_xCell2_results.tsv
│   ├── MouseRNAseqData_xCell2_results.tsv
│   ├── TabulaMurisBlood_xCell2_results.tsv
│   ├── xCell2_all_references_merged.tsv
│   └── xCell2_scaled_within_reference.tsv
└── run_module.sh
```