# Predict Target Genes of Novel miRNAs

## Usage

`bash run_module.sh`

## Folder contents

1. `01-run-miranda-human.Rmd`; Extracts human 3' UTR sequences from GENCODE v39, runs miRanda for differentially expressed miRNAs, parses and annotates predicted targets with gene IDs and symbols.

## Input files

## Analysis module directory structure

```
.
├── 01-run-miranda-human.Rmd
├── README.md
├── results
└── run_module.sh
```