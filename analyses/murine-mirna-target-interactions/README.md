# Murine miRNA–Target Interactions

## Usage

`bash run_module.sh`

## Folder contents

1. `01-murine-de-mirna-clustering.R` integrates miRNA target predictions (miRanda), immune-related GO enrichment results, and RNA-seq differential expression data to identify significant miRNA–target interactions per cluster.
2. `murine-mirna-target-dotplot.R` generates dot plots visualizing miRNA–immune target interactions for a given cluster. This script takes a single argument — the path to the interaction table and automatically detects the cluster ID for figure titles and output filenames.

## Analysis module directory structure

```
.
├── 01-murine-mirna-target-interactions.R
├── README.md
├── murine-mirna-target-dotplot.R
├── plots
│   ├── cluster1-immune-mirna-target-dotplot-filtered.pdf
│   ├── cluster1-immune-mirna-target-dotplot-full.pdf
│   ├── cluster2-immune-mirna-target-dotplot-filtered.pdf
│   ├── cluster2-immune-mirna-target-dotplot-full.pdf
│   ├── cluster5-immune-mirna-target-dotplot-filtered.pdf
│   └── cluster5-immune-mirna-target-dotplot-full.pdf
├── results
│   ├── cluster1-mirna-immune-target-interactions.tsv
│   ├── cluster2-mirna-immune-target-interactions.tsv
│   └── cluster5-mirna-immune-target-interactions.tsv
└── run_module.sh
```