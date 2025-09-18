# Murine RNA expression analysis

## Usage

`bash run_module.sh`

## Folder contents

1. `01-murine-rna-differential-expression-analysis.Rmd`: performs differential expression analysis of murine RNA-seq data using DESeq2.

## Analysis module directory structure

```
.
├── 01-murine-rna-differential-expression-analysis.Rmd
├── 01-murine-rna-differential-expression-analysis.html
├── README.md
├── input
│   └── Haydar_Mouse_RNA_miRNA_manifest_IDs_assigned.tsv
├── plots
│   ├── murine-rna-deseq2-b7h3-versus-untreated-Day14.pdf
│   ├── murine-rna-deseq2-b7h3-versus-untreated-Day21.pdf
│   ├── rna-pca-plot-all.pdf
│   └── rna-pca-plot-by-timepoint.pdf
├── results
│   └── rna-differential-expression-deseq2-b7h3-stop-vs-untreated-by-timepoint.tsv
└── run_module.sh
```