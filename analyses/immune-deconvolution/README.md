# Immune deconvolution analyses

This module runs immune deconvolution of tumor and normal bulk RNA-seq data to obtain immune cell fractions and scores, and utlizes output to
1. Compare immune cell fractions across tumor histologies and between tumor and normal samples
2. Assess correlations between DE miRNA expression and immune cell fractions and scores. 

## Usage

`bash run_module.sh`

## Folder contents

1. `01-run-immune-deconvolution.R`: runs immune deconvolution of bulk-RNA seq data using quantiseq and Xcell methods

## Analysis module directory structure

```
.
├── 01-run-immune-deconvolution.R
├── results
│   ├── quantiseq_output.rds
│   └── xcell_output.rds
└── run_module.sh
```