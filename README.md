# Differential micro-RNA expression in pediatric CNS tumors versus normal brain


## To reproduce the code in this repository:
This repository contains a docker image and code used to conduct analyses

1. Clone the repository
```
git clone git@github.com:childrens-bti/haydar-mirna.git
```

2. Pull the docker container:
```
docker pull pgc-images.sbgenomics.com/rokita-lab/haydar-mirna:v1.0.2
```
NOTE: if running on a Mac with Apple Silicon chip (M1-M4), please add `--platform linux/amd64`; otherwise add `--platform linux/arm64`

3. Start the docker container, from the `haydar-mirna` folder, run:
```
docker run --platform linux/amd64 --name <CONTAINER_NAME> -d -e PASSWORD=ANYTHING -p 8787:8787 -v $PWD:/home/rstudio/haydar-mirna pgc-images.sbgenomics.com/rokita-lab/haydar-mirna:v1.0.2
```
NOTE: if running on a Mac with Apple Silicon chip (M1-M4), please add `platform linux/amd64`

4. To execute shell within the docker image, from the `haydar-mirna` folder, run:
```
docker exec -ti <CONTAINER_NAME> bash
```

5. Run the `download-data.sh` shell script to obtain latest data files: 
```
bash download_data.sh
```

6. Navigate to an analysis module and run the shell script:
```
cd /home/rstudio/haydar-mirna/analyses/module_of_interest
```


### Below is the main directory structure listing the analyses and data files used in this repository

```
.
├── Dockerfile
├── LICENSE
├── README.md
├── analyses
├── data
│   ├── 30-931106737-miRNA-all.fpkm.xls -> v1/30-931106737-miRNA-all.fpkm.xls
│   ├── 30-931106737-miRNA_expression.xls -> v1/30-931106737-miRNA_expression.xls
│   ├── 30-963755216-miRNA-all.fpkm.xls -> v1/30-963755216-miRNA-all.fpkm.xls
│   ├── 30-963755216-miRNA_expression.xls -> v1/30-963755216-miRNA_expression.xls
│   ├── 30-992989426-RNA-TPM_values.csv -> v1/30-992989426-RNA-TPM_values.csv
│   ├── 30-992989426-RNA_raw_counts.csv -> v1/30-992989426-RNA_raw_counts.csv
│   ├── miRNA_Target_anno.xls -> v1/miRNA_Target_anno.xls
│   ├── release-notes.md -> v1/release-notes.md
│   └── v1
├── download_data.sh
├── figures
└── scripts
```

## Code Authors

Ryan Corbett ([@rjcorb](https://github.com/rjcorb)), Bicna Song ([@bicnasong](https://github.com/bicnasong))
