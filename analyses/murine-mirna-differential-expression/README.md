# Analysis of miRNA expression in CAR-T versus untreated miRNA 

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
│   ├── NovelmiRNA-251-expr.pdf
│   ├── NovelmiRNA-372.pdf
│   ├── mirna-pca-plot-all.pdf
│   ├── mirna-pca-plot-by-timepoint.pdf
│   ├── mirna-pca-plot-by-treatment.pdf
│   ├── mmu-let-7d-5p-expr.pdf
│   ├── mmu-miR-146a-5p-expr.pdf
│   ├── mmu-miR-1964-3p-expr.pdf
│   ├── mmu-miR-30c-5p.pdf
│   ├── mmu-miR-344-3p-expr.pdf
│   ├── mmu-miR-99b-5p-expr.pdf
│   ├── murine-mirna-deseq2-b7h3-versus-untreated-Day14.pdf
│   ├── murine-mirna-deseq2-b7h3-versus-untreated-Day21.pdf
│   └── murine-mirna-deseq2-b7h3-versus-untreated-Day28.pdf
├── results
│   ├── mirna-differential-expression-deseq2-b7h3-stop-vs-untreated-by-timepoint.tsv
│   └── sample-metadata.tsv
└── run_module.sh
```