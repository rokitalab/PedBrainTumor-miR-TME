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
├── plots
│   ├── B cell-proportions-byHist-xcell.pdf
│   ├── DIPG or DMG
│   │   ├── DIPG or DMG-B cell-proportions-tumor-v-adj-healthyN-xcell.pdf
│   │   ├── DIPG or DMG-B cell-proportions-tumor-v-adjN-xcell.pdf
│   │   ├── DIPG or DMG-cell-proportions-tumor-v-adj-healthyN-quantiseq.pdf
│   │   ├── DIPG or DMG-cell-proportions-tumor-v-adjN-quantiseq.pdf
│   │   ├── DIPG or DMG-Innate immune cell-proportions-tumor-v-adj-healthyN-xcell.pdf
│   │   ├── DIPG or DMG-Innate immune cell-proportions-tumor-v-adjN-xcell.pdf
│   │   ├── DIPG or DMG-Other-proportions-tumor-v-adj-healthyN-xcell.pdf
│   │   ├── DIPG or DMG-Other-proportions-tumor-v-adjN-xcell.pdf
│   │   ├── DIPG or DMG-Score-proportions-tumor-v-adj-healthyN-xcell.pdf
│   │   ├── DIPG or DMG-Score-proportions-tumor-v-adjN-xcell.pdf
│   │   ├── DIPG or DMG-T cell-proportions-tumor-v-adj-healthyN-xcell.pdf
│   │   ├── DIPG or DMG-T cell-proportions-tumor-v-adjN-xcell.pdf
│   │   ├── hsa-miR-10b-5p-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-10b-5p-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-1269a-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-1269a-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-135b-5p-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-135b-5p-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-137-3p-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-137-3p-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-139-5p-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-139-5p-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-145-5p-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-145-5p-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-203a-3p-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-203a-3p-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-206-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-206-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-212-5p-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-212-5p-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-216a-5p-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-216a-5p-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-216b-5p-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-216b-5p-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-217-5p-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-217-5p-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-222-3p-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-222-3p-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-301b-3p-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-301b-3p-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-3059-5p-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-3059-5p-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-3065-5p-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-3065-5p-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-320c-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-320c-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-3615-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-3615-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-378d-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-378d-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-450b-5p-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-450b-5p-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-455-5p-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-455-5p-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-5010-3p-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-5010-3p-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-542-3p-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-542-3p-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-582-5p-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-582-5p-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-7-5p-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-7-5p-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-769-5p-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-769-5p-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-770-3p-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-770-3p-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-888-5p-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-888-5p-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-891a-5p-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-891a-5p-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-892a-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-892a-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-92b-3p-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── hsa-miR-92b-3p-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── NovelmiRNA-1734-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   ├── NovelmiRNA-1734-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
│   │   ├── NovelmiRNA-529-tpm-vs-quantiseq-cell-proportions-DIPG or DMG.pdf
│   │   └── NovelmiRNA-529-tpm-vs-xcell-cell-proportions-DIPG or DMG.pdf
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
│   │   ├── hsa-miR-106b-3p-tpm-vs-quantiseq-cell-proportions-MB.pdf
│   │   ├── hsa-miR-106b-3p-tpm-vs-xcell-cell-proportions-MB.pdf
│   │   ├── hsa-miR-130b-3p-tpm-vs-quantiseq-cell-proportions-MB.pdf
│   │   ├── hsa-miR-130b-3p-tpm-vs-xcell-cell-proportions-MB.pdf
│   │   ├── hsa-miR-203b-3p-tpm-vs-quantiseq-cell-proportions-MB.pdf
│   │   ├── hsa-miR-203b-3p-tpm-vs-xcell-cell-proportions-MB.pdf
│   │   ├── hsa-miR-216b-5p-tpm-vs-quantiseq-cell-proportions-MB.pdf
│   │   ├── hsa-miR-216b-5p-tpm-vs-xcell-cell-proportions-MB.pdf
│   │   ├── hsa-miR-217-5p-tpm-vs-quantiseq-cell-proportions-MB.pdf
│   │   ├── hsa-miR-217-5p-tpm-vs-xcell-cell-proportions-MB.pdf
│   │   ├── hsa-miR-378d-tpm-vs-quantiseq-cell-proportions-MB.pdf
│   │   ├── hsa-miR-378d-tpm-vs-xcell-cell-proportions-MB.pdf
│   │   ├── hsa-miR-421-tpm-vs-quantiseq-cell-proportions-MB.pdf
│   │   ├── hsa-miR-421-tpm-vs-xcell-cell-proportions-MB.pdf
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
│   │   ├── MB-T cell-proportions-tumor-v-adjN-xcell.pdf
│   │   ├── NovelmiRNA-1734-tpm-vs-quantiseq-cell-proportions-MB.pdf
│   │   ├── NovelmiRNA-1734-tpm-vs-xcell-cell-proportions-MB.pdf
│   │   ├── NovelmiRNA-1819-tpm-vs-quantiseq-cell-proportions-MB.pdf
│   │   ├── NovelmiRNA-1819-tpm-vs-xcell-cell-proportions-MB.pdf
│   │   ├── NovelmiRNA-529-tpm-vs-quantiseq-cell-proportions-MB.pdf
│   │   ├── NovelmiRNA-529-tpm-vs-xcell-cell-proportions-MB.pdf
│   │   ├── NovelmiRNA-78-tpm-vs-quantiseq-cell-proportions-MB.pdf
│   │   └── NovelmiRNA-78-tpm-vs-xcell-cell-proportions-MB.pdf
│   ├── Other-proportions-byHist-xcell.pdf
│   ├── quantiseq-cell-proportions-byHist.pdf
│   ├── quantiseq-m2-m1-macrophage-ratio-byHist.pdf
│   ├── Score-proportions-byHist-xcell.pdf
│   ├── T cell-proportions-byHist-xcell.pdf
│   └── xcell-m2-m1-macrophage-ratio-byHist.pdf
├── README.md
├── results
│   ├── DIPG or DMG-de-mirna-quantiseq-fraction-correlations.tsv
│   ├── MB-de-mirna-quantiseq-fraction-correlations.tsv
│   ├── mirna-tpm.rds
│   ├── quantiseq_output.rds
│   └── xcell_output.rds
└── run_module.sh
```