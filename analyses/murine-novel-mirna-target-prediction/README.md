# Murine miRNA target prediction analysis

## Usage

`bash run_module.sh`

## Folder contents

1. `01-run-miranda-mouse.Rmd`: Extracts mouse 3′UTR sequences from Ensembl GTF, runs miRanda for differentially expressed miRNAs, parses and annotates predicted targets with gene IDs and symbols.

## Analysis module directory structure

```
.
├── 01-run-miranda-mouse.Rmd
├── README.md
├── input
│   ├── 30-963755216-all_miRNA.fa
│   └── mouse_sig_DE_miRNA_list.csv
├── results
│   ├── mouse_DE_miRNA.fa
│   └── mouse_miranda_output_parsed_anno.csv
└── run_module.sh

```