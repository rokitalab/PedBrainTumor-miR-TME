# Predict Target Genes of Novel miRNAs

## Usage

`bash run_module.sh`

## Folder contents

1. `00-extract-3utr-gencode-v39.Rmd`; Extracts 3' UTR sequences from GENCODE v39 and outputs a FASTA file.
2. `01_convert_all_novel_miRNA_to_fasta.Rmd`; Converts a plain-text table of novel miRNA sequences into FASTA format.
3. `02-extract-de-novel-miRNAs-to-fasta.Rmd`; Extracts DE Novel miRNA Sequences from FASTA

## Input files

## Analysis module directory structure

```
.
├── 00-extract-3utr-gencode-v39.Rmd
├── 01_convert_all_novel_miRNA_to_fasta.Rmd
├── 02-extract-de-novel-miRNAs-to-fasta.Rmd
├── README.md
├── results
│   ├── gencode.v39.3utr.fa
│   ├── all_novel_miRNA_mature.fa
│   ├── novel_de_miRNAs_mature.fa
└── run_module.sh
```