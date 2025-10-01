# Human miRNA targets GO enrichment analysis

This module identifies and annotates mRNA targets of DE miRNAs across tumor histologies and control contrasts, performs GO enrichment on DE target sets, and visualizes enriched terms.

## Usage

`bash run_module.sh`

## Folder contents

1. `01-merge-mirna-rna-differential-expression.Rmd`; This script identifies and annotates mRNA targets of DE miRNAs across tumor histologies and control contrasts, then saves the combined interaction tables for downstream analysis.
2. `02-GO_enrichment.Rmd`; runs gene ontology enrichment analyses on DE miRNA target genes.
3. `03-plot-enriched-terms.Rmd`; plots significantly enriched GO terms for each DE miRNA target list. 
4. `04-cluster6-go-enrichment.R`; plot enriched immune-related GO terms among DIPG/DMG cluster 6 miRNA targets

## Input files

miRTarBase files were pulled from [this database](https://awi.cuhk.edu.cn/~miRTarBase/miRTarBase_2025/), and included predicted miRNA targets separated by strong evidence (SE) versus weak evidence (WE) and experimental evidence type (W = western blot, R = reported assay, Clip = CLIP-seq).

##Analysis module directory structure

```
.
├── 01-merge-mirna-rna-differential-expression.Rmd
├── 01-merge-mirna-rna-differential-expression.html
├── 02-GO_enrichment.Rmd
├── 02-GO_enrichment.html
├── 03-plot-enriched-terms.Rmd
├── 03-plot-enriched-terms.html
├── 04-cluster6-go-enrichment.R
├── README.md
├── input
│   ├── miRTarBase_SE_R.csv
│   ├── miRTarBase_SE_W.csv
│   ├── miRTarBase_SE_WR.csv
│   ├── miRTarBase_WE_Clip.tsv.gz
│   └── miRTarBase_WE_Other.csv
├── plots
│   ├── DIPG-DMG-cluster-6-immune-go-term-dotplot.pdf
│   ├── hsa-miR-103a-3p-EPN-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-105-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-106b-3p-MB-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-10a-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-10a-5p-EPN-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-10a-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-10b-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-10b-5p-EPN-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-10b-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-1180-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-1197-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-1226-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-1229-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-124-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-1252-5p-EPN-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-129-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-130a-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-130b-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-130b-3p-MB-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-132-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-133a-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-135b-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-138-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-138-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-139-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-143-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-145-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-146b-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-148a-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-148a-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-149-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-152-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-155-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-16-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-17-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-181b-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-182-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-183-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-18a-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-1908-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-193a-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-195-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-196a-5p-EPN-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-199a-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-199a-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-199b-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-19a-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-203a-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-204-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-206-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-206-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-20a-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-21-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-21-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-2113-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-214-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-216a-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-216a-5p-MB-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-216b-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-216b-5p-MB-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-221-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-222-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-224-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-2355-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-2355-3p-MB-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-25-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-28-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-29a-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-301a-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-301b-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-3065-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-3065-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-3074-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-3074-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-31-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-3117-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-3144-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-3167-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-32-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-3200-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-3200-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-320b-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-320c-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-320d-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-320d-MB-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-323b-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-324-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-324-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-328-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-330-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-338-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-342-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-3613-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-3615-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-3615-EPN-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-3616-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-362-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-378d-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-378d-MB-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-378i-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-378i-MB-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-383-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-3943-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-3943-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-410-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-421-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-421-MB-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-424-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-433-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-433-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-449a-EPN-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-449c-5p-EPN-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-450b-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-454-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-455-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-4732-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-4742-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-4787-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-483-3p-EPN-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-490-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-490-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-491-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-491-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-493-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-496-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-5002-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-501-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-5010-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-502-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-503-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-504-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-505-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-511-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-511-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-516a-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-532-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-542-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-548a-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-548ae-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-548ah-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-548av-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-556-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-574-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-576-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-582-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-615-5p-EPN-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-615-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-625-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-625-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-642a-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-6507-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-668-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-671-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-675-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-675-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-7-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-7-5p-EPN-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-7-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-708-5p-EPN-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-744-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-760-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-760-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-766-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-767-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-769-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-769-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-873-5p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-876-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-885-3p-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-892a-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-92b-3p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-93-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-935-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-935-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-940-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   ├── hsa-miR-944-LGG-targets-immune-go-term-enrichment-dotplot.pdf
│   └── hsa-miR-96-5p-DIPG or DMG-targets-immune-go-term-enrichment-dotplot.pdf
├── results
│   ├── DIPG-or-DMG-cluster6-mirna-target-go-enr-immune-terms-full.tsv
│   ├── DIPG-or-DMG-cluster6-mirna-target-go-enr-immune-terms-reduced.tsv
│   ├── hsa-miR-103a-3p-EPN-shared-GO-BP_all.tsv
│   ├── hsa-miR-105-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-10a-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-10a-5p-EPN-shared-GO-BP_all.tsv
│   ├── hsa-miR-10a-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-10b-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-10b-5p-EPN-shared-GO-BP_all.tsv
│   ├── hsa-miR-10b-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-1180-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-1197-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-1226-3p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-1229-3p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-124-3p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-1252-5p-EPN-shared-GO-BP_all.tsv
│   ├── hsa-miR-1269a-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-129-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-130a-3p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-130b-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-132-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-133a-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-135b-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-136-3p-MB-shared-GO-BP_all.tsv
│   ├── hsa-miR-138-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-138-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-139-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-143-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-145-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-146b-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-148a-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-148a-3p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-149-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-152-3p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-155-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-16-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-17-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-181b-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-182-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-183-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-18a-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-1908-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-193a-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-195-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-196a-5p-EPN-shared-GO-BP_all.tsv
│   ├── hsa-miR-199a-3p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-199a-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-199b-3p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-19a-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-203a-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-204-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-206-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-206-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-20a-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-21-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-21-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-2113-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-214-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-216a-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-216a-5p-MB-shared-GO-BP_all.tsv
│   ├── hsa-miR-216b-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-216b-5p-MB-shared-GO-BP_all.tsv
│   ├── hsa-miR-221-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-222-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-224-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-2355-3p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-2355-3p-MB-shared-GO-BP_all.tsv
│   ├── hsa-miR-25-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-28-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-29a-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-301a-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-301b-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-3065-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-3065-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-3074-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-3074-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-31-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-3117-3p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-3144-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-3167-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-32-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-3200-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-3200-3p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-320b-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-320c-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-320d-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-320d-MB-shared-GO-BP_all.tsv
│   ├── hsa-miR-323b-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-324-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-324-3p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-328-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-330-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-338-3p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-342-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-3613-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-3615-EPN-shared-GO-BP_all.tsv
│   ├── hsa-miR-3616-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-362-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-378d-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-378d-MB-shared-GO-BP_all.tsv
│   ├── hsa-miR-378i-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-378i-MB-shared-GO-BP_all.tsv
│   ├── hsa-miR-383-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-3943-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-3943-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-410-3p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-421-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-421-MB-shared-GO-BP_all.tsv
│   ├── hsa-miR-424-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-433-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-433-3p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-449a-EPN-shared-GO-BP_all.tsv
│   ├── hsa-miR-449c-5p-EPN-shared-GO-BP_all.tsv
│   ├── hsa-miR-450b-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-454-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-455-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-4732-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-4742-3p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-483-3p-EPN-shared-GO-BP_all.tsv
│   ├── hsa-miR-490-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-490-3p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-491-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-491-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-493-3p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-496-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-5002-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-500a-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-501-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-5010-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-502-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-503-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-504-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-505-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-511-3p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-511-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-516a-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-532-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-542-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-548a-3p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-548ae-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-548ah-3p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-548av-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-556-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-574-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-576-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-582-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-615-5p-EPN-shared-GO-BP_all.tsv
│   ├── hsa-miR-615-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-625-3p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-625-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-6507-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-660-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-671-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-6733-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-675-3p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-675-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-7-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-7-5p-EPN-shared-GO-BP_all.tsv
│   ├── hsa-miR-7-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-708-5p-EPN-shared-GO-BP_all.tsv
│   ├── hsa-miR-744-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-760-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-760-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-766-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-767-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-769-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-769-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-873-5p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-876-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-885-3p-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-888-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-892a-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-92b-3p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-93-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-935-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── hsa-miR-935-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-940-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-944-LGG-shared-GO-BP_all.tsv
│   ├── hsa-miR-96-5p-DIPG or DMG-shared-GO-BP_all.tsv
│   ├── human-mirna-target-predictions-miRTarBase.tsv.gz
│   ├── known-mirna-target-gene-differential-expr-DIPG or DMG-healthyNormal.tsv
│   ├── known-mirna-target-gene-differential-expr-DIPG or DMG-paired.tsv
│   ├── known-mirna-target-gene-differential-expr-EPN-healthyNormal.tsv
│   ├── known-mirna-target-gene-differential-expr-EPN-paired.tsv
│   ├── known-mirna-target-gene-differential-expr-LGG-healthyNormal.tsv
│   ├── known-mirna-target-gene-differential-expr-MB-healthyNormal.tsv
│   └── known-mirna-target-gene-differential-expr-MB-paired.tsv
└── run_module.sh
```