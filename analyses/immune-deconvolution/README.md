# Immune deconvolution analyses

This module runs immune deconvolution of tumor and normal bulk RNA-seq data to obtain immune cell fractions and scores, and utlizes output to
1. Compare immune cell fractions across tumor histologies and between tumor and normal samples
2. Assess correlations between DE miRNA expression and immune cell fractions and scores. 

## Usage

`bash run_module.sh`

## Folder contents

1. `01-run-immune-deconvolution.R`: runs immune deconvolution of bulk-RNA seq data using quantiseq and Xcell methods
2. `02-immune-deconv-summary.Rmd`; plot immune cell fraction and scores by histology and sample type (Tumor, adj and healthy normal)
3. `03-mirna-quantiseq-fraction-correlations.Rmd`; calculate correlations between DE miRNA TPM and quantiseq-derived immune cell fractions.
4. `04-mirna-xcell-score-correlations.Rmd`; calculate correlations between DE miRNA TPM and Xcell-derived cell scores.
5. `05-mirna-xcell-score-correlations.Rmd`; calculate correlations between DE miRNA TPM and T cell marker gene expression.
6. `06-plot-cell-fraction-heatmap.R`; plot cell fraction heatmaps

## Analysis module directory structure

```
.
├── 01-run-immune-deconvolution.R
├── 02-immune-deconv-summary.html
├── 02-immune-deconv-summary.Rmd
├── 03-mirna-quantiseq-fraction-correlations.html
├── 03-mirna-quantiseq-fraction-correlations.Rmd
├── 04-mirna-xcell-score-correlations.html
├── 04-mirna-xcell-score-correlations.Rmd
├── 05-tcell-marker-gene-correlations.Rmd
├── 06-plot-cell-fraction-heatmap.R
├── plots
│   ├── B cell-proportions-byHist-xcell.pdf
│   ├── DIPG or DMG
│   ├── EPN
│   ├── Innate immune cell-proportions-byHist-xcell.pdf
│   ├── LGG
│   ├── MB
│   ├── Other-proportions-byHist-xcell.pdf
│   ├── Quantiseq_scaled_heatmap_clustered.pdf
│   ├── quantiseq-cell-proportions-byHist.pdf
│   ├── quantiseq-m2-m1-macrophage-ratio-byHist.pdf
│   ├── Score-proportions-byHist-xcell.pdf
│   ├── T cell-proportions-byHist-xcell.pdf
│   ├── XCell_scaled_heatmap_clustered.pdf
│   └── xcell-m2-m1-macrophage-ratio-byHist.pdf
├── README.md
├── results
│   ├── DIPG or DMG-de-mirna-quantiseq-fraction-correlations.tsv
│   ├── DIPG or DMG-de-mirna-tcell-marker-gene-correlations.tsv
│   ├── DIPG or DMG-de-mirna-xcell-fraction-correlations.tsv
│   ├── MB-de-mirna-quantiseq-fraction-correlations.tsv
│   ├── MB-de-mirna-tcell-marker-gene-correlations.tsv
│   ├── MB-de-mirna-xcell-fraction-correlations.tsv
│   ├── mirna-tpm.rds
│   ├── quantiseq_output.rds
│   ├── xcell_output.rds
│   └── xCell_scaled_fraction_mat.tsv
└── run_module.sh
```