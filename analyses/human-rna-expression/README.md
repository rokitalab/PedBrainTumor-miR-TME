# Human miRNA expression analysis

## Usage

`bash run_module.sh`

## Folder contents

2. `02-differential-expression.Rmd`: Run DESeq2 differential expression analyses on tumors versus matched normal and healthy normals. 
3. `03-GO-enrichment.Rmd`: Run GO enrichment using TopGO on tumor DEGs

## Analysis module directory structure

```
.
├── 02-differential-expression.html
├── 02-differential-expression.Rmd
├── 03-GO-enrichment.html
├── 03-GO-enrichment.Rmd
├── plots
│   ├── DIPG or DMG-downregulated-genes-healthyNormal-go-term-enrichment-dotplot.pdf
│   ├── DIPG or DMG-downregulated-genes-paired-go-term-enrichment-dotplot.pdf
│   ├── DIPG or DMG-upregulated-genes-healthyNormal-go-term-enrichment-dotplot.pdf
│   ├── DIPG or DMG-upregulated-genes-paired-go-term-enrichment-dotplot.pdf
│   ├── EPN-downregulated-genes-healthyNormal-go-term-enrichment-dotplot.pdf
│   ├── EPN-downregulated-genes-paired-go-term-enrichment-dotplot.pdf
│   ├── EPN-upregulated-genes-healthyNormal-go-term-enrichment-dotplot.pdf
│   ├── EPN-upregulated-genes-paired-go-term-enrichment-dotplot.pdf
│   ├── LGG-downregulated-genes-healthyNormal-go-term-enrichment-dotplot.pdf
│   ├── LGG-upregulated-genes-healthyNormal-go-term-enrichment-dotplot.pdf
│   ├── MB-downregulated-genes-paired-go-term-enrichment-dotplot.pdf
│   ├── MB-upregulated-genes-healthyNormal-go-term-enrichment-dotplot.pdf
│   ├── MB-upregulated-genes-paired-go-term-enrichment-dotplot.pdf
│   ├── volcano_ATRT_vs_adjNormal.pdf
│   ├── volcano_ATRT_vs_healthyNormal.pdf
│   ├── volcano_DIPG or DMG_paired.pdf
│   ├── volcano_DIPG or DMG_vs_adjNormal.pdf
│   ├── volcano_DIPG or DMG_vs_healthyNormal.pdf
│   ├── volcano_EPN_paired.pdf
│   ├── volcano_EPN_vs_adjNormal.pdf
│   ├── volcano_EPN_vs_healthyNormal.pdf
│   ├── volcano_HGG_paired.pdf
│   ├── volcano_HGG_vs_adjNormal.pdf
│   ├── volcano_HGG_vs_healthyNormal.pdf
│   ├── volcano_LGG_vs_healthyNormal.pdf
│   ├── volcano_MB_paired.pdf
│   ├── volcano_MB_vs_adjNormal.pdf
│   └── volcano_MB_vs_healthyNormal.pdf
├── README.md
├── results
│   ├── DE_summary_counts_paired.csv
│   ├── DE_summary_counts.csv
│   ├── DESeq2_ATRT_vs_adjNormal_sig.csv
│   ├── DESeq2_ATRT_vs_adjNormal.csv
│   ├── DESeq2_ATRT_vs_healthyNormal_sig.csv
│   ├── DESeq2_ATRT_vs_healthyNormal.csv
│   ├── DESeq2_DIPG or DMG_paired_full.csv
│   ├── DESeq2_DIPG or DMG_paired_sig.csv
│   ├── DESeq2_DIPG or DMG_vs_adjNormal_sig.csv
│   ├── DESeq2_DIPG or DMG_vs_adjNormal.csv
│   ├── DESeq2_DIPG or DMG_vs_healthyNormal_sig.csv
│   ├── DESeq2_DIPG or DMG_vs_healthyNormal.csv
│   ├── DESeq2_EPN_paired_full.csv
│   ├── DESeq2_EPN_paired_sig.csv
│   ├── DESeq2_EPN_vs_adjNormal_sig.csv
│   ├── DESeq2_EPN_vs_adjNormal.csv
│   ├── DESeq2_EPN_vs_healthyNormal_sig.csv
│   ├── DESeq2_EPN_vs_healthyNormal.csv
│   ├── DESeq2_HGG_paired_full.csv
│   ├── DESeq2_HGG_paired_sig.csv
│   ├── DESeq2_HGG_vs_adjNormal_sig.csv
│   ├── DESeq2_HGG_vs_adjNormal.csv
│   ├── DESeq2_HGG_vs_healthyNormal_sig.csv
│   ├── DESeq2_HGG_vs_healthyNormal.csv
│   ├── DESeq2_LGG_vs_healthyNormal_sig.csv
│   ├── DESeq2_LGG_vs_healthyNormal.csv
│   ├── DESeq2_MB_paired_full.csv
│   ├── DESeq2_MB_paired_sig.csv
│   ├── DESeq2_MB_vs_adjNormal_sig.csv
│   ├── DESeq2_MB_vs_adjNormal.csv
│   ├── DESeq2_MB_vs_healthyNormal_sig.csv
│   ├── DESeq2_MB_vs_healthyNormal.csv
│   ├── DIPG or DMG-downregulated-genes-healthyNormal-enriched-go-terms.tsv
│   ├── DIPG or DMG-downregulated-genes-paired-enriched-go-terms.tsv
│   ├── DIPG or DMG-upregulated-genes-healthyNormal-enriched-go-terms.tsv
│   ├── DIPG or DMG-upregulated-genes-paired-enriched-go-terms.tsv
│   ├── EPN-downregulated-genes-healthyNormal-enriched-go-terms.tsv
│   ├── EPN-downregulated-genes-paired-enriched-go-terms.tsv
│   ├── EPN-upregulated-genes-healthyNormal-enriched-go-terms.tsv
│   ├── EPN-upregulated-genes-paired-enriched-go-terms.tsv
│   ├── HGG-downregulated-genes-healthyNormal-enriched-go-terms.tsv
│   ├── HGG-upregulated-genes-healthyNormal-enriched-go-terms.tsv
│   ├── HGG-upregulated-genes-paired-enriched-go-terms.tsv
│   ├── LGG-downregulated-genes-healthyNormal-enriched-go-terms.tsv
│   ├── LGG-upregulated-genes-healthyNormal-enriched-go-terms.tsv
│   ├── MB-downregulated-genes-healthyNormal-enriched-go-terms.tsv
│   ├── MB-downregulated-genes-paired-enriched-go-terms.tsv
│   ├── MB-upregulated-genes-healthyNormal-enriched-go-terms.tsv
│   └── MB-upregulated-genes-paired-enriched-go-terms.tsv
└── run_module.sh
```