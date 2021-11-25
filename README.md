
<!-- README.md is generated from README.Rmd. Please edit that file -->

# GEGVICshine

<!-- badges: start -->
<!-- badges: end -->

A shiny app created with the intention to provide a user-friendly access
to the `GEGVIC` [package](https://github.com/oriolarques/GEGVIC).
Together they create a workflow to analyse **G**ene **E**xpression,
**G**enetic **V**ariations and **I**mmune cell **C**omposition of tumour
samples using Next Generation Sequencing data.

Presently GEGVICshine can only be used locally. For that use the
following command on your terminal:

    git clone https://github.com/oriolarques/GEGVICshine

Or manually download the repository (**app.R** file, **R/** folder and
**wwww/** folder).

In R, since the `GEGVIC` package requires many dependencies, it is
recommended to execute the following code before the first usage to
prepare the environment correctly.

``` r
# CRAN packages 
if(!require(shiny)) install.packages("shiny")
if(!require(dplyr)) install.packages("dplyr")
if(!require(tibble)) install.packages("tibble")
if(!require(tidyr)) install.packages("tidyr")
if(!require(ggplot2)) install.packages("ggplot2")
if(!require(ggrepel)) install.packages("ggrepel")
if(!require(rlang)) install.packages("rlang")
if(!require(ggplotify)) install.packages("ggplotify")
if(!require(ggpubr)) install.packages("ggpubr")
if(!require(patchwork)) install.packages("patchwork")
if(!require(pheatmap)) install.packages("pheatmap")
if(!require(devtools)) install.packages("devtools")
if(!require(remotes)) install.packages("remotes")
if(!require(DT)) install.packages("DT")
if(!require(shinyFiles)) install.packages("shinyFiles")
if(!require(shinythemes)) install.packages("shinythemes")
if(!require(tm)) install.packages("tm")

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if(!require(DESeq2)) BiocManager::install("DESeq2")
if(!require(apeglm)) BiocManager::install("apeglm")
if(!require(maftools)) BiocManager::install("maftools")
if(!require(clusterProfiler)) BiocManager::install("clusterProfiler")
if(!require(GSEAmining)) BiocManager::install("GSEAmining")
if(!require(BSgenome.Hsapiens.UCSC.hg19)) BiocManager::install("BSgenome")
if(!require(BSgenome.Hsapiens.UCSC.hg19)) BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
if(!require(BSgenome.Hsapiens.UCSC.hg19)) BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
if(!require(BSgenome.Hsapiens.UCSC.hg19)) BiocManager::install("BSgenome.Mmusculus.UCSC.mm10")
if(!require(BSgenome.Hsapiens.UCSC.hg19)) BiocManager::install("BSgenome.Mmusculus.UCSC.mm39")

# Github packages
remotes::install_github("icbi-lab/immunedeconv")
devtools::install_github('raerose01/deconstructSigs')
devtools::install_github("oriolarques/GEGVIC")
```

Then run the `GEGVICshine` app using the following code:

    shiny::runApp('path to the app.R file')
