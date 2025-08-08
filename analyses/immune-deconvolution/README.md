# Immune deconvolution analyses

This module runs immune deconvolution of tumor and normal bulk RNA-seq data to obtain immune cell fractions and scores, and utlizes output to
1. Compare immune cell fractions across tumor histologies and between tumor and normal samples
2. Assess correlations between DE miRNA expression and immune cell fractions and scores. 

## Usage

`bash run_module.sh`

## Folder contents

1. `01-run-immune-deconvolution.R`: runs immune deconvolution of bulk-RNA seq data using quantiseq and Xcell methods
2. `02-immune-deconv-summary.Rmd`; plot immune cell fraction and scores by histology and sample type (Tumor, adj and healthy normal)

## Analysis module directory structure

```
├── 01-run-immune-deconvolution.R
├── 02-immune-deconv-summary.html
├── 02-immune-deconv-summary.Rmd
├── plots
│   ├── B cell-proportions-byHist-xcell.pdf
│   ├── DMG
│   │   ├── DMG-B cell-proportions-tumor-v-adj-healthyN-xcell.pdf
│   │   ├── DMG-B cell-proportions-tumor-v-adjN-xcell.pdf
│   │   ├── DMG-cell-proportions-tumor-v-adj-healthyN-quantiseq.pdf
│   │   ├── DMG-cell-proportions-tumor-v-adjN-quantiseq.pdf
│   │   ├── DMG-Innate immune cell-proportions-tumor-v-adj-healthyN-xcell.pdf
│   │   ├── DMG-Innate immune cell-proportions-tumor-v-adjN-xcell.pdf
│   │   ├── DMG-Other-proportions-tumor-v-adj-healthyN-xcell.pdf
│   │   ├── DMG-Other-proportions-tumor-v-adjN-xcell.pdf
│   │   ├── DMG-Score-proportions-tumor-v-adj-healthyN-xcell.pdf
│   │   ├── DMG-Score-proportions-tumor-v-adjN-xcell.pdf
│   │   ├── DMG-T cell-proportions-tumor-v-adj-healthyN-xcell.pdf
│   │   └── DMG-T cell-proportions-tumor-v-adjN-xcell.pdf
│   ├── EPN
│   │   ├── EPN-B cell-proportions-tumor-v-adj-healthyN-xcell.pdf
│   │   ├── EPN-B cell-proportions-tumor-v-adjN-xcell.pdf
│   │   ├── EPN-cell-proportions-tumor-v-adj-healthyN-quantiseq.pdf
│   │   ├── EPN-cell-proportions-tumor-v-adjN-quantiseq.pdf
│   │   ├── EPN-Innate immune cell-proportions-tumor-v-adj-healthyN-xcell.pdf
│   │   ├── EPN-Innate immune cell-proportions-tumor-v-adjN-xcell.pdf
│   │   ├── EPN-Other-proportions-tumor-v-adj-healthyN-xcell.pdf
│   │   ├── EPN-Other-proportions-tumor-v-adjN-xcell.pdf
│   │   ├── EPN-Score-proportions-tumor-v-adj-healthyN-xcell.pdf
│   │   ├── EPN-Score-proportions-tumor-v-adjN-xcell.pdf
│   │   ├── EPN-T cell-proportions-tumor-v-adj-healthyN-xcell.pdf
│   │   └── EPN-T cell-proportions-tumor-v-adjN-xcell.pdf
│   ├── Innate immune cell-proportions-byHist-xcell.pdf
│   ├── LGG
│   │   ├── LGG-B cell-proportions-tumor-v-adj-healthyN-xcell.pdf
│   │   ├── LGG-cell-proportions-tumor-v-adj-healthyN-quantiseq.pdf
│   │   ├── LGG-Innate immune cell-proportions-tumor-v-adj-healthyN-xcell.pdf
│   │   ├── LGG-Other-proportions-tumor-v-adj-healthyN-xcell.pdf
│   │   ├── LGG-Score-proportions-tumor-v-adj-healthyN-xcell.pdf
│   │   └── LGG-T cell-proportions-tumor-v-adj-healthyN-xcell.pdf
│   ├── MB
│   │   ├── MB-B cell-proportions-tumor-v-adj-healthyN-xcell.pdf
│   │   ├── MB-B cell-proportions-tumor-v-adjN-xcell.pdf
│   │   ├── MB-cell-proportions-tumor-v-adj-healthyN-quantiseq.pdf
│   │   ├── MB-cell-proportions-tumor-v-adjN-quantiseq.pdf
│   │   ├── MB-Innate immune cell-proportions-tumor-v-adj-healthyN-xcell.pdf
│   │   ├── MB-Innate immune cell-proportions-tumor-v-adjN-xcell.pdf
│   │   ├── MB-Other-proportions-tumor-v-adj-healthyN-xcell.pdf
│   │   ├── MB-Other-proportions-tumor-v-adjN-xcell.pdf
│   │   ├── MB-Score-proportions-tumor-v-adj-healthyN-xcell.pdf
│   │   ├── MB-Score-proportions-tumor-v-adjN-xcell.pdf
│   │   ├── MB-T cell-proportions-tumor-v-adj-healthyN-xcell.pdf
│   │   └── MB-T cell-proportions-tumor-v-adjN-xcell.pdf
│   ├── Other-proportions-byHist-xcell.pdf
│   ├── quantiseq-cell-proportions-byHist.pdf
│   ├── Score-proportions-byHist-xcell.pdf
│   └── T cell-proportions-byHist-xcell.pdf
├── README.md
├── results
│   ├── quantiseq_output.rds
│   └── xcell_output.rds
└── run_module.sh
```