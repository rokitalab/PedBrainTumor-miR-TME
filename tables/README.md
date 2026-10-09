# Write manuscript tables

Module authors: Ryan Corbett (@rjcorb)

The purpose of this module is to generate supplementary tables for the manuscript.

## Usage
### script to run analysis
<br>**Run the R script from this directory to generate the final supplementary tables below**
```
Rscript --vanilla make-supp-tables.R
```

The script uses inputs in `input/` and existing results from the human and murine analysis modules under `../analyses/`. R and the `zip` command are required.

## Output
`TableS1.xlsx` Human sample histology, clinical metadata, and RNA-seq/miRNA-seq assay availability <br>
`TableS2.xlsx` Human miRNA differential expression results comparing tumors with healthy or paired adjacent normal tissue <br>
`TableS3.xlsx` Human miRNA cluster membership and correlations with immune scores for DIPG/DMG and medulloblastoma <br>
`TableS4.xlsx` Oncogenic and tumor-suppressive miRNAs <br>
`TableS5.xlsx` Immune-related GO term enrichment results for miRNA targets in DIPG/DMG cluster 6 and medulloblastoma cluster 1 <br>
`TableS6.xlsx` Human RNA differential expression results comparing tumors with healthy or paired adjacent normal tissue <br>
`TableS7.xlsx` Murine sample metadata and matched RNA-seq/miRNA-seq assay availability, excluding Day 28 <br>
`TableS8.xlsx` Murine miRNA and RNA differential expression summaries for B7H3 CAR and control CAR versus untreated samples at Days 14 and 21 <br>
`TableS9.xlsx` Murine miRNA cluster membership and expression patterns <br>
`TableS10.xlsx` GO term enrichment results for genes upregulated after B7H3 treatment at Days 14 and 21 <br>

## Folder content
* `make-supp-tables.R` script to generate supplementary tables in xlsx format for the manuscript
* `input/` metadata and oncogenic/tumor-suppressive miRNA reference files

## Directory structure
```
.
├── input
│   ├── cbtn_all_20251105164014_trt_removed.csv
│   ├── histologies.tsv
│   └── onco-ts-mirs.txt
├── make-supp-tables.R
├── README.md
├── TableS1.xlsx
├── TableS2.xlsx
├── TableS3.xlsx
├── TableS4.xlsx
├── TableS5.xlsx
├── TableS6.xlsx
├── TableS7.xlsx
├── TableS8.xlsx
├── TableS9.xlsx
└── TableS10.xlsx
```
