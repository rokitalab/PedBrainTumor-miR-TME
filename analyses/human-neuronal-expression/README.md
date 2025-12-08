# Human neuronal expression analyses

## Usage

`bash run_module.sh`

## Folder contents

1. `01-human-neuronal-expression.Rmd`: performs neuronal gene–level comparisons in tumor vs. normal brain samples across DIPG/DMG, MB, EPN, and LGG cohorts. The script checks for significant expression differences using the provided gene sets.

## Analysis module directory structure

```
.
├── 01-human-neuronal-expression.Rmd
├── README.md
├── plots
│   ├── neuronal_DIPG_or_DMG_healthyNormal.pdf
│   ├── neuronal_DIPG_or_DMG_paired.pdf
│   ├── neuronal_EPN_healthyNormal.pdf
│   ├── neuronal_EPN_paired.pdf
│   ├── neuronal_LGG_healthyNormal.pdf
│   ├── neuronal_MB_healthyNormal.pdf
│   └── neuronal_MB_paired.pdf
├── results
└── run_module.sh
```