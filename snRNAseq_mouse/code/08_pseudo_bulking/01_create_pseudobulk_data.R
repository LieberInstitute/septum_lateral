library(SingleCellExperiment)
library(sessioninfo)
library(spatialLIBD)
library(scater)
library(here)

############################### Load sce object ###############################

load(here("snRNAseq_mouse","processed_data","SCE","sce_updated_LS.rda"),
    verbose = TRUE)
# Loading objects:
#     sce.ls
#     annotationTab.ls
#     cell_colors.ls

sce.ls
# class: SingleCellExperiment
# dim: 32285 22860
# metadata(1): Samples
# assays(3): counts binomial_pearson_residuals logcounts
# rownames(32285): ENSMUSG00000051951 ENSMUSG00000089699 ...
# ENSMUSG00000095019 ENSMUSG00000095041
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(22860): 1_AAACCCAAGGTACATA-1 1_AAACCCACATCCGAGC-1 ...
# 4_TTTGTTGCATACAGCT-1 4_TTTGTTGGTCAAACGG-1
# colData names(17): Sample Barcode ... cellType.final cellType.broad
# reducedDimNames(4): GLMPCA_approx UMAP TSNE GLMPCA_50
# mainExpName: NULL
# altExpNames(0):

###############################################################################



###################### Function to perform pseudobulking ######################

do_pseudobulk <- function(sce, cell_cluster) {
    sce_pseudo <-
        registration_pseudobulk(sce,
            var_registration = cell_cluster,
            var_sample_id = "Sample",
            min_ncells = 10
        )

    ## Compute PCs
    pca <- prcomp(t(assays(sce_pseudo)$logcounts))
    metadata(sce_pseudo) <- list("PCA_var_explained" = (summary(pca))$importance[2, 1:20])
    metadata(sce_pseudo)
    pca_pseudo <- pca$x[, seq_len(20)]
    colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
    reducedDims(sce_pseudo) <- list(PCA = pca_pseudo)

    ## Compute some reduced dims
    set.seed(20230509)
    sce_pseudo <- scater::runMDS(sce_pseudo, ncomponents = 20)
    sce_pseudo <- scater::runPCA(sce_pseudo, name = "runPCA")

    return(sce_pseudo)
}

###############################################################################



############################## Filter sce object ##############################

## sce.ls.filter has all cell types except the ones that needed to be dropped
## sce.ls.LS has only LS clusters

#Drop the "drop" and "mixed" clusters.
all_names <- levels(colData(sce.ls)$cellType.broad)[-grep("drop|mixed",
                                                          levels(colData(sce.ls)$cellType.final))]
all_names <- all_names[-grep("mixed", all_names)]

#Pull LS clusters.
LS_names <- levels(colData(sce.ls)$cellType.final)[grep("LS",
                                                        levels(colData(sce.ls)$cellType.final))]
#LS + septal clusters.
LSept_names <- levels(colData(sce.ls)$cellType.final)[grep("LS|Sept",
                                                           levels(colData(sce.ls)$cellType.final))]
#All neuronal clusters.
neuronal_names <- levels(colData(sce.ls)$cellType.final)[grep("LS|Sept|MS|TT|TNoS",
                                                              levels(colData(sce.ls)$cellType.final))]
neuronalbr_names <- c("Chol", "LS", "Sept", "Str", "MS", "TNoS", "TT.IG.SH", "Thal", "IoC")

#Filter all of the SCE objects.
sce.ls.filter         <- sce.ls[,sce.ls$cellType.broad %in% all_names]
dim(sce.ls.filter)
# [1] 32285 21884

sce.ls.LS             <- sce.ls[,sce.ls$cellType.final %in% LS_names]
dim(sce.ls.LS)
# [1] 32285  1841

sce.ls.LS.Sept        <- sce.ls[,sce.ls$cellType.final %in% LSept_names]
dim(sce.ls.LS.Sept)
# [1] 32285  3728

sce.ls.neuronal       <- sce.ls[,sce.ls$cellType.final %in% neuronal_names]
dim(sce.ls.neuronal)
# [1] 32285  5919

sce.ls.neuronal.broad <- sce.ls[,sce.ls$cellType.broad %in% neuronalbr_names]
dim(sce.ls.neuronal.broad)
# [1] 32285 12245
###############################################################################



############################# pseudobulking for LS ############################

## pseudobulking for broad cell type clusters
sce_pseudo_all <- do_pseudobulk(sce.ls.filter, "cellType.broad")
# 2023-11-06 16:29:27.06635 make pseudobulk object
# 2023-11-06 16:29:30.347203 dropping 5 pseudo-bulked samples that are below 'min_ncells'.
# 2023-11-06 16:29:30.386218 drop lowly expressed genes
# 2023-11-06 16:29:30.603524 normalize expression
# Warning in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth = TRUE,  :
#                           You're computing too large a percentage of total singular values, use a standard svd instead.
#

## pseudobulking across LS clusters
sce_pseudo_LS <- do_pseudobulk(sce.ls.LS, "cellType.final")
# 2023-11-06 16:30:31.905901 make pseudobulk object
# 2023-11-06 16:30:32.685797 dropping 4 pseudo-bulked samples that are below 'min_ncells'.
# 2023-11-06 16:30:32.71675 drop lowly expressed genes
# 2023-11-06 16:30:32.814695 normalize expression
# Warning in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth = TRUE,  :
#                           You're computing too large a percentage of total singular values, use a standard svd instead.
# Warning message:
# In check_numbers(k = k, nu = nu, nv = nv, limit = min(dim(x)) -  :
#   more singular values/vectors requested than available


## pseudobulking across LS and Sept clusters
sce_pseudo_LS.Sept <- do_pseudobulk(sce.ls.LS.Sept, "cellType.final")
# 2023-11-06 16:30:51.801836 make pseudobulk object
# 2023-11-06 16:30:52.676572 dropping 6 pseudo-bulked samples that are below 'min_ncells'.
# 2023-11-06 16:30:52.705673 drop lowly expressed genes
# 2023-11-06 16:30:52.797302 normalize expression
# Warning in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth = TRUE,  :
#                           You're computing too large a percentage of total singular values, use a standard svd instead.
# Warning message:
# In check_numbers(k = k, nu = nu, nv = nv, limit = min(dim(x)) -  :
#   more singular values/vectors requested than available

## pseudobulking for neuronal clusters
sce_pseudo_neuronal <- do_pseudobulk(sce.ls.neuronal, "cellType.final")
# 2023-11-06 16:31:29.561479 make pseudobulk object
# 2023-11-06 16:31:31.001732 dropping 8 pseudo-bulked samples that are below 'min_ncells'.
# 2023-11-06 16:31:31.047799 drop lowly expressed genes
# 2023-11-06 16:31:31.146155 normalize expression
# Warning in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth = TRUE,  :
#                           You're computing too large a percentage of total singular values, use a standard svd instead.

## pseudobulking for neuronal broad clusters
sce_pseudo_neuronal.broad <- do_pseudobulk(sce.ls.neuronal.broad, "cellType.broad")
# 2023-11-06 16:31:47.902001 make pseudobulk object
# 2023-11-06 16:31:49.967362 dropping 1 pseudo-bulked samples that are below 'min_ncells'.
# 2023-11-06 16:31:49.99616 drop lowly expressed genes
# 2023-11-06 16:31:50.085289 normalize expression
# Warning in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth = TRUE,  :
#                           You're computing too large a percentage of total singular values, use a standard svd instead.
# Warning message:
# In check_numbers(k = k, nu = nu, nv = nv, limit = min(dim(x)) -  :
#   more singular values/vectors requested than available

###############################################################################



#################### Save both pseudobulking results to rda ###################

save(sce_pseudo_all,
     sce_pseudo_LS,
     sce_pseudo_LS.Sept,
     sce_pseudo_neuronal,
     sce_pseudo_neuronal.broad,
     file = here("snRNAseq_mouse","processed_data","SCE",
                 "sce_pseudobulking_LS_and_broad_110623.rda"))

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
# [1] "Reproducibility information:"
# [1] "2023-11-06 16:34:24 EST"
# user  system elapsed
# 75.550   5.648 615.975
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.1 Patched (2023-07-19 r84711)
# os       Rocky Linux 9.2 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2023-11-06
# pandoc   3.1.3 @ /jhpce/shared/community/core/conda_R/4.3/bin/pandoc
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────
# package                * version   date (UTC) lib source
# abind                    1.4-5     2016-07-21 [2] CRAN (R 4.3.1)
# AnnotationDbi            1.62.2    2023-07-02 [2] Bioconductor
# AnnotationHub            3.8.0     2023-04-25 [2] Bioconductor
# attempt                  0.3.1     2020-05-03 [2] CRAN (R 4.3.1)
# beachmat                 2.16.0    2023-04-25 [2] Bioconductor
# beeswarm                 0.4.0     2021-06-01 [2] CRAN (R 4.3.1)
# benchmarkme              1.0.8     2022-06-12 [2] CRAN (R 4.3.1)
# benchmarkmeData          1.0.4     2020-04-23 [2] CRAN (R 4.3.1)
# Biobase                * 2.60.0    2023-04-25 [2] Bioconductor
# BiocFileCache            2.8.0     2023-04-25 [2] Bioconductor
# BiocGenerics           * 0.46.0    2023-04-25 [2] Bioconductor
# BiocIO                   1.10.0    2023-04-25 [2] Bioconductor
# BiocManager              1.30.22   2023-08-08 [2] CRAN (R 4.3.1)
# BiocNeighbors            1.18.0    2023-04-25 [2] Bioconductor
# BiocParallel             1.34.2    2023-05-22 [2] Bioconductor
# BiocSingular             1.16.0    2023-04-25 [2] Bioconductor
# BiocVersion              3.17.1    2022-11-04 [2] Bioconductor
# Biostrings               2.68.1    2023-05-16 [2] Bioconductor
# bit                      4.0.5     2022-11-15 [2] CRAN (R 4.3.1)
# bit64                    4.0.5     2020-08-30 [2] CRAN (R 4.3.1)
# bitops                   1.0-7     2021-04-24 [2] CRAN (R 4.3.1)
# blob                     1.2.4     2023-03-17 [2] CRAN (R 4.3.1)
# bslib                    0.5.1     2023-08-11 [2] CRAN (R 4.3.1)
# cachem                   1.0.8     2023-05-01 [2] CRAN (R 4.3.1)
# cli                      3.6.1     2023-03-23 [2] CRAN (R 4.3.1)
# codetools                0.2-19    2023-02-01 [3] CRAN (R 4.3.1)
# colorout               * 1.2-2     2023-09-22 [1] Github (jalvesaq/colorout@79931fd)
# colorspace               2.1-0     2023-01-23 [2] CRAN (R 4.3.1)
# config                   0.3.2     2023-08-30 [2] CRAN (R 4.3.1)
# cowplot                  1.1.1     2020-12-30 [1] CRAN (R 4.3.1)
# crayon                   1.5.2     2022-09-29 [2] CRAN (R 4.3.1)
# curl                     5.0.2     2023-08-14 [2] CRAN (R 4.3.1)
# data.table               1.14.8    2023-02-17 [2] CRAN (R 4.3.1)
# DBI                      1.1.3     2022-06-18 [2] CRAN (R 4.3.1)
# dbplyr                   2.3.3     2023-07-07 [2] CRAN (R 4.3.1)
# DelayedArray             0.26.7    2023-07-28 [2] Bioconductor
# DelayedMatrixStats       1.22.6    2023-08-28 [2] Bioconductor
# digest                   0.6.33    2023-07-07 [2] CRAN (R 4.3.1)
# doParallel               1.0.17    2022-02-07 [2] CRAN (R 4.3.1)
# dotCall64                1.0-2     2022-10-03 [2] CRAN (R 4.3.1)
# dplyr                    1.1.3     2023-09-03 [2] CRAN (R 4.3.1)
# dqrng                    0.3.1     2023-08-30 [2] CRAN (R 4.3.1)
# DropletUtils             1.20.0    2023-04-25 [2] Bioconductor
# DT                       0.29      2023-08-29 [2] CRAN (R 4.3.1)
# edgeR                    3.42.4    2023-05-31 [2] Bioconductor
# ellipsis                 0.3.2     2021-04-29 [2] CRAN (R 4.3.1)
# ExperimentHub            2.8.1     2023-07-12 [2] Bioconductor
# fansi                    1.0.4     2023-01-22 [2] CRAN (R 4.3.1)
# fastmap                  1.1.1     2023-02-24 [2] CRAN (R 4.3.1)
# fields                   15.2      2023-08-17 [2] CRAN (R 4.3.1)
# filelock                 1.0.2     2018-10-05 [2] CRAN (R 4.3.1)
# foreach                  1.5.2     2022-02-02 [2] CRAN (R 4.3.1)
# generics                 0.1.3     2022-07-05 [2] CRAN (R 4.3.1)
# GenomeInfoDb           * 1.36.3    2023-09-07 [2] Bioconductor
# GenomeInfoDbData         1.2.10    2023-07-20 [2] Bioconductor
# GenomicAlignments        1.36.0    2023-04-25 [2] Bioconductor
# GenomicRanges          * 1.52.0    2023-04-25 [2] Bioconductor
# ggbeeswarm               0.7.2     2023-04-29 [2] CRAN (R 4.3.1)
# ggplot2                * 3.4.3     2023-08-14 [2] CRAN (R 4.3.1)
# ggrepel                  0.9.3     2023-02-03 [2] CRAN (R 4.3.1)
# glue                     1.6.2     2022-02-24 [2] CRAN (R 4.3.1)
# golem                    0.4.1     2023-06-05 [2] CRAN (R 4.3.1)
# gridExtra                2.3       2017-09-09 [2] CRAN (R 4.3.1)
# gtable                   0.3.4     2023-08-21 [2] CRAN (R 4.3.1)
# HDF5Array                1.28.1    2023-05-01 [2] Bioconductor
# here                   * 1.0.1     2020-12-13 [2] CRAN (R 4.3.1)
# htmltools                0.5.6     2023-08-10 [2] CRAN (R 4.3.1)
# htmlwidgets              1.6.2     2023-03-17 [2] CRAN (R 4.3.1)
# httpuv                   1.6.11    2023-05-11 [2] CRAN (R 4.3.1)
# httr                     1.4.7     2023-08-15 [2] CRAN (R 4.3.1)
# interactiveDisplayBase   1.38.0    2023-04-25 [2] Bioconductor
# IRanges                * 2.34.1    2023-06-22 [2] Bioconductor
# irlba                    2.3.5.1   2022-10-03 [2] CRAN (R 4.3.1)
# iterators                1.0.14    2022-02-05 [2] CRAN (R 4.3.1)
# jquerylib                0.1.4     2021-04-26 [2] CRAN (R 4.3.1)
# jsonlite                 1.8.7     2023-06-29 [2] CRAN (R 4.3.1)
# KEGGREST                 1.40.0    2023-04-25 [2] Bioconductor
# later                    1.3.1     2023-05-02 [2] CRAN (R 4.3.1)
# lattice                  0.21-8    2023-04-05 [3] CRAN (R 4.3.1)
# lazyeval                 0.2.2     2019-03-15 [2] CRAN (R 4.3.1)
# lifecycle                1.0.3     2022-10-07 [2] CRAN (R 4.3.1)
# limma                    3.56.2    2023-06-04 [2] Bioconductor
# locfit                   1.5-9.8   2023-06-11 [2] CRAN (R 4.3.1)
# magick                   2.7.5     2023-08-07 [2] CRAN (R 4.3.1)
# magrittr                 2.0.3     2022-03-30 [2] CRAN (R 4.3.1)
# maps                     3.4.1     2022-10-30 [2] CRAN (R 4.3.1)
# Matrix                   1.6-1.1   2023-09-18 [3] CRAN (R 4.3.1)
# MatrixGenerics         * 1.12.3    2023-07-30 [2] Bioconductor
# matrixStats            * 1.0.0     2023-06-02 [2] CRAN (R 4.3.1)
# memoise                  2.0.1     2021-11-26 [2] CRAN (R 4.3.1)
# mime                     0.12      2021-09-28 [2] CRAN (R 4.3.1)
# munsell                  0.5.0     2018-06-12 [2] CRAN (R 4.3.1)
# paletteer                1.5.0     2022-10-19 [2] CRAN (R 4.3.1)
# pillar                   1.9.0     2023-03-22 [2] CRAN (R 4.3.1)
# pkgconfig                2.0.3     2019-09-22 [2] CRAN (R 4.3.1)
# plotly                   4.10.2    2023-06-03 [2] CRAN (R 4.3.1)
# png                      0.1-8     2022-11-29 [2] CRAN (R 4.3.1)
# promises                 1.2.1     2023-08-10 [2] CRAN (R 4.3.1)
# purrr                    1.0.2     2023-08-10 [2] CRAN (R 4.3.1)
# R.methodsS3              1.8.2     2022-06-13 [2] CRAN (R 4.3.1)
# R.oo                     1.25.0    2022-06-12 [2] CRAN (R 4.3.1)
# R.utils                  2.12.2    2022-11-11 [2] CRAN (R 4.3.1)
# R6                       2.5.1     2021-08-19 [2] CRAN (R 4.3.1)
# rappdirs                 0.3.3     2021-01-31 [2] CRAN (R 4.3.1)
# RColorBrewer             1.1-3     2022-04-03 [2] CRAN (R 4.3.1)
# Rcpp                     1.0.11    2023-07-06 [2] CRAN (R 4.3.1)
# RCurl                    1.98-1.12 2023-03-27 [2] CRAN (R 4.3.1)
# rematch2                 2.1.2     2020-05-01 [2] CRAN (R 4.3.1)
# restfulr                 0.0.15    2022-06-16 [2] CRAN (R 4.3.1)
# rhdf5                    2.44.0    2023-04-25 [2] Bioconductor
# rhdf5filters             1.12.1    2023-04-30 [2] Bioconductor
# Rhdf5lib                 1.22.1    2023-09-10 [2] Bioconductor
# rjson                    0.2.21    2022-01-09 [2] CRAN (R 4.3.1)
# rlang                    1.1.1     2023-04-28 [2] CRAN (R 4.3.1)
# rprojroot                2.0.3     2022-04-02 [2] CRAN (R 4.3.1)
# Rsamtools                2.16.0    2023-04-25 [2] Bioconductor
# RSQLite                  2.3.1     2023-04-03 [2] CRAN (R 4.3.1)
# rsvd                     1.0.5     2021-04-16 [2] CRAN (R 4.3.1)
# rtracklayer              1.60.1    2023-08-15 [2] Bioconductor
# S4Arrays                 1.0.6     2023-08-30 [2] Bioconductor
# S4Vectors              * 0.38.1    2023-05-02 [2] Bioconductor
# sass                     0.4.7     2023-07-15 [2] CRAN (R 4.3.1)
# ScaledMatrix             1.8.1     2023-05-03 [2] Bioconductor
# scales                   1.2.1     2022-08-20 [2] CRAN (R 4.3.1)
# scater                 * 1.28.0    2023-04-25 [2] Bioconductor
# scuttle                * 1.10.2    2023-08-03 [2] Bioconductor
# sessioninfo            * 1.2.2     2021-12-06 [2] CRAN (R 4.3.1)
# shiny                    1.7.5     2023-08-12 [2] CRAN (R 4.3.1)
# shinyWidgets             0.8.0     2023-08-30 [2] CRAN (R 4.3.1)
# SingleCellExperiment   * 1.22.0    2023-04-25 [2] Bioconductor
# spam                     2.9-1     2022-08-07 [2] CRAN (R 4.3.1)
# sparseMatrixStats        1.12.2    2023-07-02 [2] Bioconductor
# SpatialExperiment      * 1.10.0    2023-04-25 [2] Bioconductor
# spatialLIBD            * 1.12.0    2023-04-27 [2] Bioconductor
# statmod                  1.5.0     2023-01-06 [2] CRAN (R 4.3.1)
# SummarizedExperiment   * 1.30.2    2023-06-06 [2] Bioconductor
# tibble                   3.2.1     2023-03-20 [2] CRAN (R 4.3.1)
# tidyr                    1.3.0     2023-01-24 [2] CRAN (R 4.3.1)
# tidyselect               1.2.0     2022-10-10 [2] CRAN (R 4.3.1)
# utf8                     1.2.3     2023-01-31 [2] CRAN (R 4.3.1)
# vctrs                    0.6.3     2023-06-14 [2] CRAN (R 4.3.1)
# vipor                    0.4.5     2017-03-22 [2] CRAN (R 4.3.1)
# viridis                  0.6.4     2023-07-22 [2] CRAN (R 4.3.1)
# viridisLite              0.4.2     2023-05-02 [2] CRAN (R 4.3.1)
# withr                    2.5.0     2022-03-03 [2] CRAN (R 4.3.1)
# XML                      3.99-0.14 2023-03-19 [2] CRAN (R 4.3.1)
# xtable                   1.8-4     2019-04-21 [2] CRAN (R 4.3.1)
# XVector                  0.40.0    2023-04-25 [2] Bioconductor
# yaml                     2.3.7     2023-01-23 [2] CRAN (R 4.3.1)
# zlibbioc                 1.46.0    2023-04-25 [2] Bioconductor
#
# [1] /users/rphillip/R/4.3
# [2] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────
