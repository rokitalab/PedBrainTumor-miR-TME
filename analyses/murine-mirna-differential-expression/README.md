# Analysis of miRNA expression in CAR-T versus untreated mice 

This analysis module performs miRNA differential expression in mice treated with B7H3 or control CAR-T versus untreated mice over three time points (14, 21, and 28 days post-treatment)

## Usage

`bash run_module.sh`

## Folder contents

1. `01-mirna-differential-expression.Rmd`; Performs miRNA differential expression analysis of CAR-T versus untreated mice at each time point
2. `02-GO_enrichment.Rmd`; perform GO enrichment analysis on DE miRNA targets
3. `03-plot_gsea.Rmd`; plots DE miRNA target gene GO enrichment

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
│   ├── NovelmiRNA-1265-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── NovelmiRNA-1309-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── NovelmiRNA-325-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mirna-pca-plot-all.pdf
│   ├── mirna-pca-plot-by-timepoint.pdf
│   ├── mirna-pca-plot-by-treatment.pdf
│   ├── mmu-let-7a-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-let-7b-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-let-7c-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-let-7d-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-let-7j-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-106b-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-107-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-10a-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-10b-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-1191a-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-1198-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-128-2-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-128-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-133b-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-135b-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-138-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-146a-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-148a-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-151-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-155-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-17-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-17-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-181a-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-1843a-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-185-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-191-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-1943-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-1947-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-1964-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-196b-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-1981-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-204-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-210-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-21a-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-22-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-23a-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-25-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-26a-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-296-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-299a-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-30b-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-30c-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-30e-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-32-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-330-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-344-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-361-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-362-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-376b-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-376c-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-380-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-384-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-421-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-423-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-467d-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-490-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-500-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-501-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-503-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-5126-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-532-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-668-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-669c-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-671-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-672-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-674-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-6958-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-7043-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-7224-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-760-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-7689-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-770-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-8120-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-9-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-92a-1-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-92a-2-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-93-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── murine-mirna-deseq2-b7h3-versus-untreated-Day14.pdf
│   ├── murine-mirna-deseq2-b7h3-versus-untreated-Day21.pdf
│   └── murine-mirna-deseq2-b7h3-versus-untreated-Day28.pdf
├── results
│   ├── NovelmiRNA-1265-target-go-enrichment.tsv
│   ├── NovelmiRNA-1309-target-go-enrichment.tsv
│   ├── NovelmiRNA-325-target-go-enrichment.tsv
│   ├── mirna-differential-expression-deseq2-b7h3-stop-vs-untreated-by-timepoint.tsv
│   ├── mmu-let-7a-5p-target-go-enrichment.tsv
│   ├── mmu-let-7b-3p-target-go-enrichment.tsv
│   ├── mmu-let-7c-5p-target-go-enrichment.tsv
│   ├── mmu-let-7d-5p-target-go-enrichment.tsv
│   ├── mmu-let-7j-target-go-enrichment.tsv
│   ├── mmu-miR-106b-5p-target-go-enrichment.tsv
│   ├── mmu-miR-107-3p-target-go-enrichment.tsv
│   ├── mmu-miR-10a-5p-target-go-enrichment.tsv
│   ├── mmu-miR-10b-5p-target-go-enrichment.tsv
│   ├── mmu-miR-1191a-target-go-enrichment.tsv
│   ├── mmu-miR-1198-5p-target-go-enrichment.tsv
│   ├── mmu-miR-126a-3p-target-go-enrichment.tsv
│   ├── mmu-miR-128-2-3p-target-go-enrichment.tsv
│   ├── mmu-miR-128-3p-target-go-enrichment.tsv
│   ├── mmu-miR-133b-3p-target-go-enrichment.tsv
│   ├── mmu-miR-135b-5p-target-go-enrichment.tsv
│   ├── mmu-miR-138-5p-target-go-enrichment.tsv
│   ├── mmu-miR-146a-5p-target-go-enrichment.tsv
│   ├── mmu-miR-148a-3p-target-go-enrichment.tsv
│   ├── mmu-miR-151-5p-target-go-enrichment.tsv
│   ├── mmu-miR-155-5p-target-go-enrichment.tsv
│   ├── mmu-miR-17-3p-target-go-enrichment.tsv
│   ├── mmu-miR-17-5p-target-go-enrichment.tsv
│   ├── mmu-miR-181a-5p-target-go-enrichment.tsv
│   ├── mmu-miR-1843a-5p-target-go-enrichment.tsv
│   ├── mmu-miR-185-5p-target-go-enrichment.tsv
│   ├── mmu-miR-191-5p-target-go-enrichment.tsv
│   ├── mmu-miR-1943-5p-target-go-enrichment.tsv
│   ├── mmu-miR-1947-5p-target-go-enrichment.tsv
│   ├── mmu-miR-1964-3p-target-go-enrichment.tsv
│   ├── mmu-miR-196b-5p-target-go-enrichment.tsv
│   ├── mmu-miR-1981-5p-target-go-enrichment.tsv
│   ├── mmu-miR-204-5p-target-go-enrichment.tsv
│   ├── mmu-miR-210-3p-target-go-enrichment.tsv
│   ├── mmu-miR-21a-5p-target-go-enrichment.tsv
│   ├── mmu-miR-22-3p-target-go-enrichment.tsv
│   ├── mmu-miR-23a-3p-target-go-enrichment.tsv
│   ├── mmu-miR-25-3p-target-go-enrichment.tsv
│   ├── mmu-miR-26a-5p-target-go-enrichment.tsv
│   ├── mmu-miR-296-3p-target-go-enrichment.tsv
│   ├── mmu-miR-299a-3p-target-go-enrichment.tsv
│   ├── mmu-miR-30b-5p-target-go-enrichment.tsv
│   ├── mmu-miR-30c-5p-target-go-enrichment.tsv
│   ├── mmu-miR-30e-5p-target-go-enrichment.tsv
│   ├── mmu-miR-32-5p-target-go-enrichment.tsv
│   ├── mmu-miR-330-5p-target-go-enrichment.tsv
│   ├── mmu-miR-344-3p-target-go-enrichment.tsv
│   ├── mmu-miR-361-3p-target-go-enrichment.tsv
│   ├── mmu-miR-362-5p-target-go-enrichment.tsv
│   ├── mmu-miR-376b-3p-target-go-enrichment.tsv
│   ├── mmu-miR-376c-3p-target-go-enrichment.tsv
│   ├── mmu-miR-380-3p-target-go-enrichment.tsv
│   ├── mmu-miR-384-3p-target-go-enrichment.tsv
│   ├── mmu-miR-410-3p-target-go-enrichment.tsv
│   ├── mmu-miR-421-3p-target-go-enrichment.tsv
│   ├── mmu-miR-423-3p-target-go-enrichment.tsv
│   ├── mmu-miR-467d-5p-target-go-enrichment.tsv
│   ├── mmu-miR-487b-3p-target-go-enrichment.tsv
│   ├── mmu-miR-490-3p-target-go-enrichment.tsv
│   ├── mmu-miR-500-3p-target-go-enrichment.tsv
│   ├── mmu-miR-501-3p-target-go-enrichment.tsv
│   ├── mmu-miR-503-3p-target-go-enrichment.tsv
│   ├── mmu-miR-5126-target-go-enrichment.tsv
│   ├── mmu-miR-532-5p-target-go-enrichment.tsv
│   ├── mmu-miR-668-3p-target-go-enrichment.tsv
│   ├── mmu-miR-669c-5p-target-go-enrichment.tsv
│   ├── mmu-miR-671-3p-target-go-enrichment.tsv
│   ├── mmu-miR-672-5p-target-go-enrichment.tsv
│   ├── mmu-miR-674-3p-target-go-enrichment.tsv
│   ├── mmu-miR-6958-3p-target-go-enrichment.tsv
│   ├── mmu-miR-7043-3p-target-go-enrichment.tsv
│   ├── mmu-miR-7224-3p-target-go-enrichment.tsv
│   ├── mmu-miR-760-3p-target-go-enrichment.tsv
│   ├── mmu-miR-7689-3p-target-go-enrichment.tsv
│   ├── mmu-miR-770-3p-target-go-enrichment.tsv
│   ├── mmu-miR-8120-target-go-enrichment.tsv
│   ├── mmu-miR-9-5p-target-go-enrichment.tsv
│   ├── mmu-miR-92a-1-3p-target-go-enrichment.tsv
│   ├── mmu-miR-92a-2-3p-target-go-enrichment.tsv
│   ├── mmu-miR-93-5p-target-go-enrichment.tsv
│   └── sample-metadata.tsv
└── run_module.sh
```