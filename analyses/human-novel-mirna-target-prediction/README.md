# Predict Target Genes of Novel miRNAs

## Usage

`bash run_module.sh`

## Folder contents

1. `00-extract-3utr-gencode-v39.Rmd`; Extracts 3' UTR sequences from GENCODE v39 and outputs a FASTA file.
2. `01_convert_all_novel_miRNA_to_fasta.Rmd`; Converts a plain-text table of novel miRNA sequences into FASTA format.
3. `02-run-miranda.sh`; Runs miRanda using novel DE miRNAs and 3′ UTR sequences extracted from GENCODE v39.
4. `03-parse-miranda-output.py`; Parses the raw miRanda output file to extract predicted miRNA-target interactions and saves them in a clean CSV format.
5. `04-annotate-parsed-miranda-output.Rmd`; Annotates the parsed miRanda output file with Ensembl Gene IDs and gene symbols using a GTF annotation file.

## Input files

## Analysis module directory structure

```
.
├── 00-extract-3utr-gencode-v39.Rmd
├── 00-extract-3utr-gencode-v39.html
├── 01-convert-all-novel-miRNA-to-fasta.Rmd
├── 01-convert-all-novel-miRNA-to-fasta.html 
├── 02-run-miranda.sh
├── 03-parse-miranda-output.py
├── 04-annotate-parsed-miranda-output.Rmd
├── 04-annotate-parsed-miranda-output.html
├── README.md
├── results
│   ├── comined_all_novel_mirna.fa
│   └── comined_all_novel_mirna.tsv
└── run_module.sh
```