# Murine miRNA target prediction analysis

## Usage

`bash run_module.sh`

## Folder contents

1. `01-run-miranda-mouse.Rmd`: Extracts mouse 3′UTR sequences from Ensembl GTF, runs miRanda for differentially expressed miRNAs, parses and annotates predicted targets with gene IDs and symbols.

Note: `mouse_miranda_output.txt` and `Mus_musculus.GRCm39.3utr.fa` were not uploaded to the repository due to file size.

## Analysis module directory structure

```
.
├── 01-run-miranda-mouse.Rmd
├── 01-run-miranda-mouse.html
├── README.md
├── results
│   ├── Mus_musculus.GRCm39.3utr.fa
│   ├── mouse_DE_miRNA.fa
│   ├── mouse_miranda_output.txt
│   └── mouse_miranda_output_parsed_anno.csv
└── run_module.sh
```