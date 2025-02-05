FROM rocker/tidyverse:4.4.0
LABEL maintainer="Ryan Corbett (rcorbett@childrensnational.org)"
WORKDIR /rocker-build/

RUN apt-get update && apt-get install -y --no-install-recommends apt-utils dialog

# Add curl, bzip2 and some dev libs
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    curl \
    bzip2 \
    zlib1g \
    libbz2-dev \
    liblzma-dev \
    libreadline-dev

# libmagick++-dev is needed for coloblindr to install
RUN apt-get -y --no-install-recommends install \
    libgdal-dev \
    libudunits2-dev \
    libmagick++-dev

# Required for installing pdftools, which is a dependency of gridGraphics
RUN apt-get -y --no-install-recommends install \
    libpoppler-cpp-dev

# Install java
RUN apt-get update && apt-get -y --no-install-recommends install \
   default-jdk \
   libxt6

# Set the Bioconductor repository as the primary repository
RUN R -e "options(repos = BiocManager::repositories())"

# Install BiocManager and the desired version of Bioconductor
RUN R -e "install.packages('BiocManager', dependencies=TRUE)"
RUN R -e "BiocManager::install(version = '3.19', ask = FALSE)"

# Install packages
RUN R -e 'BiocManager::install(c( \
  "AnnotationDbi", \
  "Biobase", \
  "broom", \
  "circlize", \
  "clusterProfiler", \
  "ComplexHeatmap", \
  "corrplot", \
  "cowplot", \
  "DESeq2", \
  "edgeR", \
  "EnhancedVolcano", \
  "fgsea", \
  "ggpubr", \
  "ggstatsplot", \
  "ggthemes", \
  "gridExtra", \
  "GO.db", \
  "GSVA", \
  "limma", \
  "msigdbr", \
  "optparse", \
  "org.Hs.eg.db", \
  "pheatmap", \
  "rtracklayer", \
  "R.utils", \
  "sva", \
  "UpSetR" \
))'


## install GitHub packages
RUN R -e "remotes::install_github('clauswilke/colorblindr', ref = '1ac3d4d62dad047b68bb66c06cee927a4517d678', dependencies = TRUE)"
RUN R -e "remotes::install_github('thomasp85/patchwork', ref = '1cb732b129ed6a65774796dc1f618558c7498b66', dependencies = TRUE)"

WORKDIR /rocker-build/

ADD Dockerfile .
