# Human RNA expression analysis

## Usage

`bash run_module.sh`

## Folder contents

1. `rna_seq_pca.Rmd`: investigates whether RNA-seq samples show batch effects or biological clustering based on histology and sample type (tumor, tumor-adjacent normal, and healthy normal brain). The goal is to visualize global expression patterns using PCA.

## Analysis module directory structure

```
.
├── rna_seq_pca.Rmd
├── rna_seq_pca.html
├── README.md
├── plots
│   ├── rna_seq_pca.pdf
│   ├── pca_ATRT_vs_controls.pdf
│   ├── pca_DIPG_vs_controls.pdf
│   ├── pca_Ependymoma_vs_controls.pdf
│   ├── pca_HGG_vs_controls.pdf
│   └── pca_Medulloblastoma_vs_controls.pdf
└── run_module.sh
```