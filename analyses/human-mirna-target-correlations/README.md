# Human DE miRNA - target associations

This module identifies significant associations between DE miRNAs in clusters of interest and predicted immune-related targets as follows:

1) Calculates pearson correlation coefficients between DE miRNA and target expression
2) Assesses whether predicted targets exhibit differential expression in the opposite direction of predicted targeting miRNAs

## Usage

`bash run_module.sh`

## Folder contents

1. `01-dmg-cluster6-mirna-target-interactions.R`: identify significant associations between DMG cluster6 miRNAs and predicted immune-related targets
2. `02-dmg-cluster6-mirna-target-dotplot.R`: Plot significant associations between DMG cluster6 miRNAs and predicted immune-related targets
3. `03-mb-cluster1-mirna-target-interactions.R`: identify significant associations between MB cluster1 miRNAs and predicted immune-related targets
4. `04-mb-cluster1-mirna-target-dotplot.R`: Plot significant associations between MB cluster1 miRNAs and predicted immune-related targets
5. `05-dmg-cluster6-mirna-target-pathway-correlations.R`: identify DMG cluster 6 miRNA targets in enriched immune GO terms, calculate miRNA/pathway-target expression correlations, and plot all FDR-significant correlations
6. `06-mb-cluster1-mirna-target-pathway-correlations.R`: MB cluster 1 counterpart, using matched MB miRNA-seq/RNA-seq samples and MB differential-expression contrasts

### Input files

- `dipg-dmg-cluster6-mirna-immune-target-sig-interactions.txt`: reviewed and annotated target GO terms
- `mb-cluster1-mirna-immune-target-sig-interactions.txt`: reviewed and annotated target GO terms

## Analysis module directory structure

```
.
├── 01-dmg-cluster6-mirna-target-interactions.R
├── 02-dmg-cluster6-mirna-target-dotplot.R
├── 03-mb-cluster1-mirna-target-interactions.R
├── 04-mb-cluster1-mirna-target-dotplot.R
├── 05-dmg-cluster6-mirna-target-pathway-correlations.R
├── 06-mb-cluster1-mirna-target-pathway-correlations.R
├── input
│   ├── dipg-dmg-cluster6-mirna-immune-target-sig-interactions.txt
│   └── mb-cluster1-mirna-immune-target-sig-interactions.txt
├── plots
│   ├── dmg-cluster6-target-myeloid-cell-dotplot.pdf
│   ├── dmg-cluster6-target-other-term-dotplot.pdf
│   ├── dmg-cluster6-target-t-cell-dotplot.pdf
│   ├── mb-cluster1-target-myeloid-cell-dotplot.pdf
│   ├── mb-cluster1-target-other-term-dotplot.pdf
│   ├── mb-cluster1-target-t-cell-dotplot.pdf
│   └── correlation-plots/
├── README.md
├── results
│   ├── dmg-cluster6-mirna-immune-target-interactions.tsv
│   ├── dmg-cluster6-mirna-immune-target-sig-interactions.tsv
│   ├── dmg-cluster6-mirna-immune-go-target-differential-expression.tsv
│   ├── dmg-cluster6-mirna-pathway-target-mean-pearson-correlations.tsv
│   ├── dmg-cluster6-mirna-pathway-target-mean-tpm.tsv
│   ├── mb-cluster1-mirna-immune-target-interactions.tsv
│   ├── mb-cluster1-mirna-immune-target-sig-interactions.tsv
│   ├── mb-cluster1-mirna-immune-go-target-differential-expression.tsv
│   ├── mb-cluster1-mirna-pathway-target-mean-pearson-correlations.tsv
│   └── mb-cluster1-mirna-pathway-target-mean-tpm.tsv
└── run_module.sh
```
