# Merge miRNA and RNA differential expression results in DIPG/DMG versus normal pons

## Usage

`bash run_module.sh`

## Folder contents

1. `01-merge-mirna-rna-differential-expression.Rmd`; merges DESeq2 results from miRNA and DE miRNA target genes in DIPG/DMG versus normal pons.
2. `02-GO_enrichment.Rmd`; runs gene ontology enrichment analyses on DE miRNA target genes.
3. `03-plot_gsea.Rmd`; plots signficantly enriched GO terms for each DE miRNA target list. 

## Input files

miRTarBase files were pulled from [this database](https://awi.cuhk.edu.cn/~miRTarBase/miRTarBase_2025/), and included predicted miRNA targets separated by strong evidence (SE) versus weak evidence (WE) and experimental evidence type (W = western blot, R = reported assay, Clip = CLIP-seq).

##Analysis module directory structure

```
.
├── 01-merge-mirna-rna-differential-expression.Rmd
├── 02-GO_enrichment.Rmd
├── 03-plot_gsea.Rmd
├── README.md
├── input
│   ├── miRTarBase_SE_R.csv
│   ├── miRTarBase_SE_W.csv
│   ├── miRTarBase_SE_WR.csv
│   ├── miRTarBase_WE_Clip.tsv.gz
│   └── miRTarBase_WE_Other.csv
├── plots
│   ├── hsa-miR-103a-3p-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-23c-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-330-5p-targets-immune-go-term-enrichment-dotplot.pdf
│   └── hsa-miR-361-3p-targets-immune-go-term-enrichment-dotplot.pdf
├── results
│   ├── hsa-miR-103a-3p-target-go-enrichment.tsv
│   ├── hsa-miR-23c-target-go-enrichment.tsv
│   ├── hsa-miR-330-5p-target-go-enrichment.tsv
│   ├── hsa-miR-361-3p-target-go-enrichment.tsv
│   ├── human-mirna-target-predictions-miRTarBase.tsv.gz
│   ├── known-mirna-target-gene-differential-expr-dipg-dmg-versus-normal.tsv
│   └── novel-mirna-target-gene-differential-expr-dipg-dmg-versus-normal.tsv
└── run_module.sh
```