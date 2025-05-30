# Human RNA expression analysis

## Usage

`bash run_module.sh`

## Folder contents

1. `01-merge-and-annotate-RNAseq-STAR-counts.Rmd`: merges gzipped STAR ReadsPerGene.out.tab.gz files and annotates each row with gene symbols from the Gencode v39 GTF.

## Analysis module directory structure

```
.
├── 01-merge-and-annotate-RNAseq-STAR-counts.Rmd
├── 01-merge-and-annotate-RNAseq-STAR-counts.html
├── README.md
├── input
│   ├── rnaseq_star_gene_count
│   └── Haydar-RNAseq-manifest_IDs_assigned.tsv
├── results
│   └── merged-RNA-rawcounts-annotated.tsv
└── run_module.sh
```