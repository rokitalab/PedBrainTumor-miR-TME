# GO enrichment analyses of DE miRNA targets in CAR-T versus untreated mice 

This analysis module performs GO enrichment analyses of targets of DE miRNAs in B7H3 CAR T treated versus untreated mice

## Usage

`bash run_module.sh`

## Folder contents

1. `01-GO_enrichment.Rmd`: Performs GO enrichment analysis on DE miRNA targets
2. `02-plot_gsea.Rmd`: Plots DE miRNA target gene GO enrichment results
3. `03-cluster-go-enrichment.R`: Generates immune-related GO term dot plots for murine DE miRNA clusters

##Analysis module directory structure

```
.
├── 01-GO_enrichment.Rmd
├── 01-GO_enrichment.html
├── 02-plot_gsea.Rmd
├── 02-plot_gsea.html
├── 03-cluster-go-enrichment.R
├── README.md
├── plots
│   ├── NovelmiRNA-1265-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── NovelmiRNA-1309-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── NovelmiRNA-325-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── cluster-1-immune-go-term-dotplot.pdf
│   ├── cluster-2-immune-go-term-dotplot.pdf
│   ├── cluster-3-immune-go-term-dotplot.pdf
│   ├── cluster-4-immune-go-term-dotplot.pdf
│   ├── cluster-5-immune-go-term-dotplot.pdf
│   ├── mmu-let-7a-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-let-7b-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-let-7c-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-let-7d-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-let-7j-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-106b-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-107-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-10a-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-10b-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-1191a-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-1198-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-128-2-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-128-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-133b-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-135b-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-138-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-146a-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-148a-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-151-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-155-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-17-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-17-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-181a-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-1843a-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-185-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-191-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-1943-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-1947-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-1964-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-196b-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-1981-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-204-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-210-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-21a-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-22-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-23a-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-25-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-26a-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-296-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-299a-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-30b-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-30c-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-30e-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-32-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-326-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-330-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-344-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-361-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-362-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-376b-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-376c-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-380-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-384-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-421-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-423-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-467d-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-490-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-500-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-501-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-503-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-5126-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-532-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-543-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-668-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-669c-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-671-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-672-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-6958-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-7043-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-7224-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-760-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-7689-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-770-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-8120-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-9-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-92a-1-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-92a-2-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   └── mmu-miR-93-5p-targets-immune-go-term-enrichment-dotplot.pdf
├── results
│   ├── NovelmiRNA-1265-target-go-enrichment.tsv
│   ├── NovelmiRNA-1309-target-go-enrichment.tsv
│   ├── NovelmiRNA-325-target-go-enrichment.tsv
│   ├── cluster1-mirna-target-go-enr-immune-terms-full.tsv
│   ├── cluster1-mirna-target-go-enr-immune-terms-reduced.tsv
│   ├── cluster2-mirna-target-go-enr-immune-terms-full.tsv
│   ├── cluster2-mirna-target-go-enr-immune-terms-reduced.tsv
│   ├── cluster3-mirna-target-go-enr-immune-terms-full.tsv
│   ├── cluster3-mirna-target-go-enr-immune-terms-reduced.tsv
│   ├── cluster4-mirna-target-go-enr-immune-terms-full.tsv
│   ├── cluster4-mirna-target-go-enr-immune-terms-reduced.tsv
│   ├── cluster5-mirna-target-go-enr-immune-terms-full.tsv
│   ├── cluster5-mirna-target-go-enr-immune-terms-reduced.tsv
│   ├── mmu-let-7a-5p-target-go-enrichment.tsv
│   ├── mmu-let-7b-3p-target-go-enrichment.tsv
│   ├── mmu-let-7c-5p-target-go-enrichment.tsv
│   ├── mmu-let-7d-5p-target-go-enrichment.tsv
│   ├── mmu-let-7j-target-go-enrichment.tsv
│   ├── mmu-miR-106b-5p-target-go-enrichment.tsv
│   ├── mmu-miR-107-3p-target-go-enrichment.tsv
│   ├── mmu-miR-10a-5p-target-go-enrichment.tsv
│   ├── mmu-miR-10b-5p-target-go-enrichment.tsv
│   ├── mmu-miR-1191a-target-go-enrichment.tsv
│   ├── mmu-miR-1198-5p-target-go-enrichment.tsv
│   ├── mmu-miR-126a-3p-target-go-enrichment.tsv
│   ├── mmu-miR-128-2-3p-target-go-enrichment.tsv
│   ├── mmu-miR-128-3p-target-go-enrichment.tsv
│   ├── mmu-miR-133b-3p-target-go-enrichment.tsv
│   ├── mmu-miR-135b-5p-target-go-enrichment.tsv
│   ├── mmu-miR-138-5p-target-go-enrichment.tsv
│   ├── mmu-miR-146a-5p-target-go-enrichment.tsv
│   ├── mmu-miR-148a-3p-target-go-enrichment.tsv
│   ├── mmu-miR-151-5p-target-go-enrichment.tsv
│   ├── mmu-miR-155-5p-target-go-enrichment.tsv
│   ├── mmu-miR-17-3p-target-go-enrichment.tsv
│   ├── mmu-miR-17-5p-target-go-enrichment.tsv
│   ├── mmu-miR-181a-5p-target-go-enrichment.tsv
│   ├── mmu-miR-1843a-5p-target-go-enrichment.tsv
│   ├── mmu-miR-185-5p-target-go-enrichment.tsv
│   ├── mmu-miR-191-5p-target-go-enrichment.tsv
│   ├── mmu-miR-1943-5p-target-go-enrichment.tsv
│   ├── mmu-miR-1947-5p-target-go-enrichment.tsv
│   ├── mmu-miR-1964-3p-target-go-enrichment.tsv
│   ├── mmu-miR-196b-5p-target-go-enrichment.tsv
│   ├── mmu-miR-1981-5p-target-go-enrichment.tsv
│   ├── mmu-miR-204-5p-target-go-enrichment.tsv
│   ├── mmu-miR-210-3p-target-go-enrichment.tsv
│   ├── mmu-miR-21a-5p-target-go-enrichment.tsv
│   ├── mmu-miR-22-3p-target-go-enrichment.tsv
│   ├── mmu-miR-23a-3p-target-go-enrichment.tsv
│   ├── mmu-miR-25-3p-target-go-enrichment.tsv
│   ├── mmu-miR-26a-5p-target-go-enrichment.tsv
│   ├── mmu-miR-296-3p-target-go-enrichment.tsv
│   ├── mmu-miR-299a-3p-target-go-enrichment.tsv
│   ├── mmu-miR-30b-5p-target-go-enrichment.tsv
│   ├── mmu-miR-30c-5p-target-go-enrichment.tsv
│   ├── mmu-miR-30e-5p-target-go-enrichment.tsv
│   ├── mmu-miR-32-5p-target-go-enrichment.tsv
│   ├── mmu-miR-326-3p-target-go-enrichment.tsv
│   ├── mmu-miR-330-5p-target-go-enrichment.tsv
│   ├── mmu-miR-344-3p-target-go-enrichment.tsv
│   ├── mmu-miR-361-3p-target-go-enrichment.tsv
│   ├── mmu-miR-362-5p-target-go-enrichment.tsv
│   ├── mmu-miR-365-1-3p-target-go-enrichment.tsv
│   ├── mmu-miR-365-2-3p-target-go-enrichment.tsv
│   ├── mmu-miR-376b-3p-target-go-enrichment.tsv
│   ├── mmu-miR-376c-3p-target-go-enrichment.tsv
│   ├── mmu-miR-380-3p-target-go-enrichment.tsv
│   ├── mmu-miR-384-3p-target-go-enrichment.tsv
│   ├── mmu-miR-410-3p-target-go-enrichment.tsv
│   ├── mmu-miR-421-3p-target-go-enrichment.tsv
│   ├── mmu-miR-423-3p-target-go-enrichment.tsv
│   ├── mmu-miR-467d-5p-target-go-enrichment.tsv
│   ├── mmu-miR-487b-3p-target-go-enrichment.tsv
│   ├── mmu-miR-490-3p-target-go-enrichment.tsv
│   ├── mmu-miR-500-3p-target-go-enrichment.tsv
│   ├── mmu-miR-501-3p-target-go-enrichment.tsv
│   ├── mmu-miR-503-3p-target-go-enrichment.tsv
│   ├── mmu-miR-5126-target-go-enrichment.tsv
│   ├── mmu-miR-532-5p-target-go-enrichment.tsv
│   ├── mmu-miR-543-3p-target-go-enrichment.tsv
│   ├── mmu-miR-668-3p-target-go-enrichment.tsv
│   ├── mmu-miR-669c-5p-target-go-enrichment.tsv
│   ├── mmu-miR-671-3p-target-go-enrichment.tsv
│   ├── mmu-miR-672-5p-target-go-enrichment.tsv
│   ├── mmu-miR-6958-3p-target-go-enrichment.tsv
│   ├── mmu-miR-7043-3p-target-go-enrichment.tsv
│   ├── mmu-miR-7224-3p-target-go-enrichment.tsv
│   ├── mmu-miR-760-3p-target-go-enrichment.tsv
│   ├── mmu-miR-7689-3p-target-go-enrichment.tsv
│   ├── mmu-miR-770-3p-target-go-enrichment.tsv
│   ├── mmu-miR-8120-target-go-enrichment.tsv
│   ├── mmu-miR-9-5p-target-go-enrichment.tsv
│   ├── mmu-miR-92a-1-3p-target-go-enrichment.tsv
│   ├── mmu-miR-92a-2-3p-target-go-enrichment.tsv
│   └── mmu-miR-93-5p-target-go-enrichment.tsv
└── run_module.sh
```