# Merge miRNA and RNA differential expression results in DIPG/DMG versus normal pons

## Usage

`bash run_module.sh`

## Folder contents

1. `01-merge-mirna-rna-differential-expression.Rmd`; merges DESeq2 results from miRNA and DE miRNA target genes in DIPG/DMG versus normal pons.

##Analysis module directory structure

```
.
├── 01-merge-mirna-rna-differential-expression.Rmd
├── README.md
├── plots
├── results
│   └── mirna-target-gene-differential-expr-dipg-dmg-versus-normal.tsv
└── run_module.sh
```