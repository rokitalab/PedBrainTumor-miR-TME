# Analysis of miRNA expression in CAR-T versus untreated mice 

This analysis module performs miRNA differential expression in mice treated with B7H3 or control CAR-T versus untreated mice over three time points (14, 21, and 28 days post-treatment)

## Usage

`bash run_module.sh`

## Folder contents

1. `01-mirna-differential-expression.Rmd`; Performs miRNA differential expression analysis of CAR-T versus untreated mice at each time point.
2. `02-GO_enrichment.Rmd`; perform GO enrichment analysis on DE miRNA targets

##Analysis module directory structure

```
.
├── 01-mirna-differential-expression.Rmd
├── 01-mirna-differential-expression.html
├── 02-GO_enrichment.Rmd
├── 02-GO_enrichment.html
├── 03-plot_gsea.Rmd
├── 03-plot_gsea.html
├── README.md
├── plots
│   ├── mirna-pca-plot-all.pdf
│   ├── mirna-pca-plot-by-timepoint.pdf
│   ├── mirna-pca-plot-by-treatment.pdf
│   ├── mmu-miR-1198-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-128-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-146a-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-151-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-155-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-185-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-204-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-30c-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-423-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-466i-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-679-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── murine-mirna-deseq2-b7h3-versus-untreated-Day14.pdf
│   ├── murine-mirna-deseq2-b7h3-versus-untreated-Day21.pdf
│   └── murine-mirna-deseq2-b7h3-versus-untreated-Day28.pdf
├── results
│   ├── mirna-differential-expression-deseq2-b7h3-stop-vs-untreated-by-timepoint.tsv
│   ├── mmu-miR-1198-5p-target-go-enrichment.tsv
│   ├── mmu-miR-128-3p-target-go-enrichment.tsv
│   ├── mmu-miR-146a-5p-target-go-enrichment.tsv
│   ├── mmu-miR-151-5p-target-go-enrichment.tsv
│   ├── mmu-miR-155-5p-target-go-enrichment.tsv
│   ├── mmu-miR-185-5p-target-go-enrichment.tsv
│   ├── mmu-miR-1964-3p-target-go-enrichment.tsv
│   ├── mmu-miR-1981-5p-target-go-enrichment.tsv
│   ├── mmu-miR-204-5p-target-go-enrichment.tsv
│   ├── mmu-miR-30b-5p-target-go-enrichment.tsv
│   ├── mmu-miR-30c-5p-target-go-enrichment.tsv
│   ├── mmu-miR-423-3p-target-go-enrichment.tsv
│   ├── mmu-miR-466i-3p-target-go-enrichment.tsv
│   ├── mmu-miR-679-5p-target-go-enrichment.tsv
│   ├── mouse-mirna-target-predictions-miRTarBase.tsv.gz
│   └── sample-metadata.tsv
└── run_module.sh
```