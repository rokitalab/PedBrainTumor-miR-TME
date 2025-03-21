# Differential mRNA expression in DIPG/DMG versus normal pons 

## Usage

`bash run_module.sh`

## Folder contents

1. `01-differential-expression.Rmd`; run differential expression analysis
2. `02-GO_enrichment.Rmd`; run GO enrichment on differentially expressed genes
3. `03-kegg_gsva.Rmd`; run GSVA on expression data and perform pathway differential expression analysis

##Analysis module directory structure

```
.
├── 01-differential-expression.Rmd
├── 01-differential-expression.html
├── 02-GO_enrichment.Rmd
├── 02-GO_enrichment.html
├── 03-kegg_gsva.Rmd
├── README.md
├── plots
│   ├── ACTC1_cts_dipg_dmg_versus_normal.pdf
│   ├── CDKN2B_cts_dipg_dmg_versus_normal.pdf
│   ├── COL14A1_cts_dipg_dmg_versus_normal.pdf
│   ├── COL6A3_cts_dipg_dmg_versus_normal.pdf
│   ├── CRB1_cts_dipg_dmg_versus_normal.pdf
│   ├── HS3ST3A1_cts_dipg_dmg_versus_normal.pdf
│   ├── HS3ST3B1_cts_dipg_dmg_versus_normal.pdf
│   ├── LHX2_cts_dipg_dmg_versus_normal.pdf
│   ├── MYH6_cts_dipg_dmg_versus_normal.pdf
│   ├── MYH7_cts_dipg_dmg_versus_normal.pdf
│   ├── MYL2_cts_dipg_dmg_versus_normal.pdf
│   ├── NPTX2_cts_dipg_dmg_versus_normal.pdf
│   ├── PAX3_cts_dipg_dmg_versus_normal.pdf
│   ├── PAX8-AS1_cts_dipg_dmg_versus_normal.pdf
│   ├── RPL7_cts_dipg_dmg_versus_normal.pdf
│   ├── RPS3A_cts_dipg_dmg_versus_normal.pdf
│   ├── SNORA53_cts_dipg_dmg_versus_normal.pdf
│   ├── THBS1_cts_dipg_dmg_versus_normal.pdf
│   ├── TM4SF1_cts_dipg_dmg_versus_normal.pdf
│   ├── TPT1_cts_dipg_dmg_versus_normal.pdf
│   ├── de-mirna-dipg-dmg-versus-normal-volcano-plot.pdf
│   ├── de-pathway-expr-heatmap.pdf
│   ├── dipg-dmg-upregulated-genes-go-term-enrichment-dotplot.pdf
│   ├── immune_modulation_degs_dotplots.pdf
│   ├── immunosuppressive_degs_dotplots.pdf
│   ├── mir330_degs_dotplots.pdf
│   ├── myeloid_recruitment_degs_dotplots.pdf
│   ├── novelmiR_degs_dotplots.pdf
│   └── rna-pca-plot.pdf
├── results
│   ├── dmg-downregulated-genes.tsv
│   ├── dmg-upregulated-genes.tsv
│   ├── enriched-go-terms-dipg-dmg-upregulated-genes.tsv
│   └── mrna-differential-expression-deseq2-dipg-dmg-versus-normal.tsv
└── run_module.sh
```