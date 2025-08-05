# Human miRNA targets GO enrichment analysis

This module identifies and annotates mRNA targets of DE miRNAs across tumor histologies and control contrasts, performs GO enrichment on DE target sets, and visualizes enriched terms.

## Usage

`bash run_module.sh`

## Folder contents

1. `01-merge-mirna-rna-differential-expression.Rmd`; This script identifies and annotates mRNA targets of DE miRNAs across tumor histologies and control contrasts, then saves the combined interaction tables for downstream analysis.
2. `02-GO_enrichment.Rmd`; runs gene ontology enrichment analyses on DE miRNA target genes.
3. `03-plot-enriched-terms.Rmd`; plots signficantly enriched GO terms for each DE miRNA target list. 

## Input files

miRTarBase files were pulled from [this database](https://awi.cuhk.edu.cn/~miRTarBase/miRTarBase_2025/), and included predicted miRNA targets separated by strong evidence (SE) versus weak evidence (WE) and experimental evidence type (W = western blot, R = reported assay, Clip = CLIP-seq).

##Analysis module directory structure

```
.
├── 01-merge-mirna-rna-differential-expression.Rmd
├── 02-GO_enrichment.Rmd
├── 03-plot-enriched-terms.Rmd
├── README.md
├── input
│   ├── miRTarBase_SE_R.csv
│   ├── miRTarBase_SE_W.csv
│   ├── miRTarBase_SE_WR.csv
│   ├── miRTarBase_WE_Clip.tsv.gz
│   └── miRTarBase_WE_Other.csv
├── plots
├── results
└── run_module.sh
```