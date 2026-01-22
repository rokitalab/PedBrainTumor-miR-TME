# DE miRNA clustering

## Usage

`bash run_module.sh`

## Folder contents

1. `01-murine-de-mirna-clustering.R`: clusters DE miRNAs from CAR-T versus untreated mice by expression
2. `02-plot-de-mirna.R`: plots DE miRNA expression across time points

## Analysis module directory structure

```
.
├── 01-murine-de-mirna-clustering.R
├── 02-plot-de-mirna.R
├── README.md
├── plots
│   ├── NovelmiRNA-1265_log2TPM_Day14_Day21.pdf
│   ├── NovelmiRNA-1309_log2TPM_Day14_Day21.pdf
│   ├── NovelmiRNA-325_log2TPM_Day14_Day21.pdf
│   ├── de-mirna-heatmap-Day14.pdf
│   ├── de-mirna-heatmap-Day21.pdf
│   ├── de-mirna-heatmap.pdf
│   ├── mmu-let-7a-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-let-7b-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-let-7c-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-let-7d-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-let-7j_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-106b-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-107-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-10a-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-10b-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-1191a_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-1198-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-126a-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-128-2-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-128-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-133b-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-135b-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-137-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-138-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-146a-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-148a-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-151-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-155-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-17-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-17-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-181a-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-1843a-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-185-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-191-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-1943-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-1947-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-1964-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-196b-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-1981-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-204-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-210-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-21a-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-22-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-23a-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-25-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-26a-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-296-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-299a-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-30b-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-30c-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-30e-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-32-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-326-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-330-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-340-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-344-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-361-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-362-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-365-1-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-365-2-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-374b-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-376b-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-376c-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-380-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-384-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-410-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-421-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-423-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-467d-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-487b-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-490-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-500-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-501-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-503-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-5126_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-532-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-543-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-668-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-669c-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-671-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-672-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-6958-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-7043-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-7224-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-760-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-7689-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-770-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-8120_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-9-5p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-92a-1-3p_log2TPM_Day14_Day21.pdf
│   ├── mmu-miR-92a-2-3p_log2TPM_Day14_Day21.pdf
│   └── mmu-miR-93-5p_log2TPM_Day14_Day21.pdf
├── results
│   ├── mouse_sig_DE_miRNA_list.csv
│   ├── mouse-de-mirna-cluster-membership.tsv
│   └── murine-mirna-tpm.rds
└── run_module.sh
```
