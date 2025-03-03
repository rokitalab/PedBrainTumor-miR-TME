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
├── input
│   ├── miRTarBase_SE_R.csv
│   ├── miRTarBase_SE_W.csv
│   ├── miRTarBase_SE_WR.csv
│   ├── miRTarBase_WE_Clip.tsv.gz
│   └── miRTarBase_WE_Other.csv
├── results
│   ├── human-mirna-target-predictions-miRTarBase.tsv.gz
│   ├── known-mirna-target-gene-differential-expr-dipg-dmg-versus-normal.tsv
│   └── novel-mirna-target-gene-differential-expr-dipg-dmg-versus-normal.tsv
└── run_module.sh
```