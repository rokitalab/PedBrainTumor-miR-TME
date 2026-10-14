# Human mouse miRNA orthology analyses 

This module identifies mouse orthologs to human miRNAs that are DE in DMG versus normal brain, and assesses their expression patterns in B7H3 treated versus untreated mice. 

## Usage

`bash run_module.sh`

## Folder contents

1. `01-summarize-orthology-expression.Rmd`: identifies mouse orthologs and assessses DMG cluster-specific expression patterns in B7H3 treated mice
2. `02-plot-cluster6-ortholog-expr.R`: Plots expression of human DMG cluster6 orthologs that are DE in B7H3 treated versus untreated mice. 
3. `03-plot-cluster6-ortholog-target-expr.R`: Identifies high-confidence miRanda-predicted targets (Max_Score >= 145, Max_Energy <= -20) of DE cluster 6 orthologs, plots target expression (mean log2(TPM + 1) per gene) in B7H3 treated versus untreated mice at Day 14 & Day 21 for all targets and for immune-related targets (same immune GO terms as `04-immune-target-expression.R`), and compares treatments within each timepoint using paired Wilcoxon signed-rank tests (paired by gene; BH-adjusted). Also compares per-gene CAR vs. untreated expression differences of predicted targets to those of all other expressed genes (Wilcoxon rank-sum tests; BH-adjusted).
4. `04-immune-target-expression.R`: Identifies high-confidence miRanda-predicted targets of DE cluster 6 orthologs that are annotated to immune-related GO biological processes (immune system process, inflammatory response, cytokine production, cytokine-mediated signaling, defense response, response to cytokine, phagocytosis, and descendant terms), and tests whether these targets are preferentially up- or downregulated in B7H3 treated versus untreated mice at Day 14 & Day 21 (binomial tests among DE targets; Fisher's exact tests versus immune non-targets; BH-adjusted).

## Analysis module directory structure

```
.
├── 01-summarize-orthology-expression.Rmd
├── 02-plot-cluster6-ortholog-expr.R
├── 03-plot-cluster6-ortholog-target-expr.R
├── 04-immune-target-expression.R
├── plots
│   ├── hsa-dmg-clust6-mouse-ortholog-de-tpm.pdf
│   ├── hsa-dmg-clust6-mouse-ortholog-immune-target-de-barplot.pdf
│   ├── hsa-dmg-clust6-mouse-ortholog-immune-target-tpm-violin.pdf
│   ├── hsa-dmg-clust6-mouse-ortholog-predicted-target-tpm-violin.pdf
│   ├── hsa-dmg-de-mirna-mmu-ortho-barplot.pdf
│   └── human-dmg-de-mirna-orthology-mouse-expr-enrichment-plot.pdf
├── README.md
├── results
│   ├── hsa-dmg-clust6-mouse-ortholog-immune-target-direction-tests.tsv
│   ├── hsa-dmg-clust6-mouse-ortholog-immune-target-expr-wilcoxon.tsv
│   ├── hsa-dmg-clust6-mouse-ortholog-immune-targets.tsv
│   ├── hsa-dmg-clust6-mouse-ortholog-predicted-targets.tsv
│   ├── hsa-dmg-clust6-mouse-ortholog-target-expr-wilcoxon.tsv
│   └── hsa-dmg-clust6-mouse-ortholog-target-vs-nontarget-wilcoxon.tsv
└── run_module.sh
```