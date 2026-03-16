# Murine RNA expression analysis

## Usage

`bash run_module.sh`

## Folder contents

1. `01-murine-rna-differential-expression-analysis.Rmd`: performs differential expression analysis of murine RNA-seq data using DESeq2.
2. `02-GO-enrichment.Rmd`: run GO enrichment analyses on differentially expressed genes, and plot significant results.
3. `03-mirna-target-correlations.Rmd`: identify downregulated miRNAs target upregulated immune genes in enriched terms identified above.

## Analysis module directory structure

```
.
├── 01-murine-rna-differential-expression-analysis.html
├── 01-murine-rna-differential-expression-analysis.Rmd
├── 02-GO-enrichment.html
├── 02-GO-enrichment.Rmd
├── 03-mirna-target-correlations.Rmd
├── input
│   └── Haydar_Mouse_RNA_miRNA_manifest_IDs_assigned.tsv
├── plots
│   ├── b7h3_up_d14-go-term-enrichment-dotplot.pdf
│   ├── b7h3_up_d21-go-term-enrichment-dotplot.pdf
│   ├── mmu-miR-30b-5p-Alox5-tpm-dotplot.pdf
│   ├── mmu-miR-421-3p-Tal1-tpm-dotplot.pdf
│   ├── mmu-miR-467d-5p-Tnfrsf9-tpm-dotplot.pdf
│   ├── murine-rna-deseq2-b7h3-versus-untreated-Day14.pdf
│   ├── murine-rna-deseq2-b7h3-versus-untreated-Day21.pdf
│   ├── rna-pca-plot-all.pdf
│   ├── rna-pca-plot-by-timepoint.pdf
│   └── rna-pca-plot-by-treatment.pdf
├── README.md
├── results
│   ├── b7h3_up_d14-enriched-go-terms.tsv
│   ├── b7h3_up_d21-enriched-go-terms.tsv
│   └── rna-differential-expression-deseq2-b7h3-stop-vs-untreated-by-timepoint.tsv
└── run_module.sh
```