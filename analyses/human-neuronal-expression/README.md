# Human neuronal expression analyses

## Usage

`bash run_module.sh`

## Folder contents

1. `01-human-neuronal-expression.Rmd`: performs neuronal gene–level comparisons in tumor vs. normal brain samples across DIPG/DMG, MB, EPN, and LGG cohorts.
2. `02-human-neuronal-mirna-target-correlation.Rmd`: performs correlation analysis between DE miRNAs and GSVA scores of neuronal-related GO terms in DIPG/DMG or MB.
3. `03-human-neuronal-heatmap.R`: regenerates the miRNA clustering heatmap with neuronal GO BP correlation coefficients as annotations in DIPG/DMG or MB.
4. `04-human-neuronal-go-enrichment.R`: plots cluster-specific miRNA target GO enrichment for neuronal-related terms in DIPG/DMG or MB.
5. `05-dmg-cluster6-mirna-target-interactions.R` identifies and summarizes predicted target gene interactions for cluster 6 miRNAs in DIPG/DMG.
6. `06-dmg-cluster6-mirna-target-dotplot.R` generates dot plots for selected neuronal target gene categories regulated by cluster 6 miRNAs in DIPG/DMG.
7. `07-mb-cluster1-mirna-target-interactions.R` identifies and summarizes predicted target gene interactions for cluster 1 miRNAs in MB.
8. `08-mb-cluster1-mirna-target-dotplot.R` generates dot plots for selected neuronal target gene categories regulated by cluster 1 miRNAs in MB.

## Analysis module directory structure

```
.
├── 01-human-neuronal-expression.Rmd
├── 01-human-neuronal-expression.html
├── 02-human-neuronal-mirna-target-correlation.Rmd
├── 02-human-neuronal-mirna-target-correlation.html
├── 03-human-neuronal-heatmap.R
├── 03-human-neuronal-heatmap.html
├── 04-human-neuronal-go-enrichment.R
├── 04-human-neuronal-go-enrichment.html
├── 05-dmg-cluster6-mirna-target-interactions.R
├── 05-dmg-cluster6-mirna-target-interactions.html
├── 06-dmg-cluster6-mirna-target-dotplot.R
├── 06-dmg-cluster6-mirna-target-dotplot.html
├── 07-mb-cluster1-mirna-target-interactions.R
├── 07-mb-cluster1-mirna-target-interactions.html
├── 08-mb-cluster1-mirna-target-dotplot.R
├── 08-mb-cluster1-mirna-target-dotplot.html
├── README.md
├── input
│   └── neuronal_go_term_with_id.tsv
├── plots
│   ├── DIPG_or_DMG-cluster6-neuronal-go-term-dotplot.pdf
│   ├── DIPG_or_DMG-cluster6-target-synaptic-dotplot.pdf
│   ├── DIPG_or_DMG_mirna_neuronal_gobp_heatmap.pdf
│   ├── MB-cluster1-neuronal-go-term-dotplot.pdf
│   ├── MB-cluster1-target-neurotransmitter-dotplot.pdf
│   ├── MB-cluster1-target-postsynaptic-dotplot.pdf
│   ├── MB-cluster1-target-synaptic-dotplot.pdf
│   ├── MB_mirna_neuronal_gobp_heatmap.pdf
│   ├── neuronal_DIPG_or_DMG_healthyNormal.pdf
│   ├── neuronal_DIPG_or_DMG_paired.pdf
│   ├── neuronal_EPN_healthyNormal.pdf
│   ├── neuronal_EPN_paired.pdf
│   ├── neuronal_LGG_healthyNormal.pdf
│   ├── neuronal_MB_healthyNormal.pdf
│   └── neuronal_MB_paired.pdf
├── results
│   ├── DIPG_or_DMG-cluster6-mirna-neuronal-target-interactions.tsv
│   ├── DIPG_or_DMG-cluster6-mirna-neuronal-target-sig-interactions.tsv
│   ├── DIPG_or_DMG-cluster6-mirna-target-go-enr-neuronal-terms-full.tsv
│   ├── DIPG_or_DMG-cluster6-mirna-target-go-enr-neuronal-terms-reduced.tsv
│   ├── DIPG_or_DMG_mirna_neuronal_cluster_membership.tsv
│   ├── DIPG_or_DMG_mirna_neuronal_gobp_spearman.tsv
│   ├── DIPG_or_DMG_mirna_neuronal_gobp_spearman_FDR01.tsv
│   ├── MB-cluster1-mirna-immune-target-interactions.tsv
│   ├── MB-cluster1-mirna-immune-target-sig-interactions.tsv
│   ├── MB-cluster1-mirna-target-go-enr-neuronal-terms-full.tsv
│   ├── MB-cluster1-mirna-target-go-enr-neuronal-terms-reduced.tsv
│   ├── MB_mirna_neuronal_cluster_membership.tsv
│   ├── MB_mirna_neuronal_gobp_spearman.tsv
│   ├── MB_mirna_neuronal_gobp_spearman_FDR01.tsv
│   └── neuronal_gobp_terms_from_gsva.tsv
└── run_module.sh
```

