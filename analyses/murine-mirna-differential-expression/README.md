# Analysis of miRNA expression in CAR-T versus untreated mice 

This analysis module performs miRNA differential expression in mice treated with B7H3 or control CAR-T versus untreated mice over three time points (14, 21, and 28 days post-treatment)

## Usage

`bash run_module.sh`

## Folder contents

1. `01-mirna-differential-expression.Rmd`; Performs miRNA differential expression analysis of CAR-T versus untreated mice at each time point.

##Analysis module directory structure

```
.
├── 01-mirna-differential-expression.Rmd
├── README.md
├── plots
│   ├── mirna-pca-plot-all.pdf
│   ├── mirna-pca-plot-by-timepoint.pdf
│   ├── mirna-pca-plot-by-treatment.pdf
│   ├── murine-mirna-deseq2-b7h3-versus-untreated-Day14.pdf
│   ├── murine-mirna-deseq2-b7h3-versus-untreated-Day21.pdf
│   └── murine-mirna-deseq2-b7h3-versus-untreated-Day28.pdf
├── results
│   ├── mirna-differential-expression-deseq2-b7h3-stop-vs-untreated-by-timepoint.tsv
│   └── sample-metadata.tsv
└── run_module.sh
```