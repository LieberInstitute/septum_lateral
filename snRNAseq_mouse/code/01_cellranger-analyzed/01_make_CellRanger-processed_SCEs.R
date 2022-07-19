### LS PRELIM snRNA-seq analysis (mouse)
### Build SCE from Cell Ranger-processed count matrices
### qrsh -l bluejay,mf=20G,h_vmem=22G
### Initiated MNT 19Jan2022

# test.edit

library(SingleCellExperiment)
library(DropletUtils)
library(rtracklayer)
library(BiocParallel)
library(scater)
library(scuttle)
library(jaffelab)
library(here)
library(sessioninfo)

here()
# [1] "/dcs04/lieber/marmaypag/pilotLS_LIBD1070"

# hm this is kinda weird bc we are in:
getwd()
# [1] "/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/code/01_cellranger-analyzed"


## These operations can definitely be built into a function to streamline, but let's
# first read in the count data from the processed (nuclei-already-called) first sample,
# '1M_C_LS':
Sys.time()
# [1] "2022-01-19 13:16:11 EST"
sce.m1 <- read10xCounts(
    samples = "/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/processed_data/cellranger/1M_C_LS/outs/filtered_feature_bc_matrix/",
    sample.names = "1M_C_LS",
    type = "sparse",
    col.names = TRUE
)
Sys.time()
# [1] "2022-01-19 13:16:35 EST"

# A summary of this object:
sce.m1
# class: SingleCellExperiment
# dim: 32285 8304
# metadata(1): Samples
# assays(1): counts
# rownames(32285): ENSMUSG00000051951 ENSMUSG00000089699 ...
#   ENSMUSG00000095019 ENSMUSG00000095041
# rowData names(3): ID Symbol Type
# colnames(8304): AAACCCAAGGTACATA-1 AAACCCACATCCGAGC-1 ...
#   TTTGTTGTCACGGAGA-1 TTTGTTGTCCAACCGG-1
# colData names(2): Sample Barcode
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):


## Currently there are just Ensembl ID, Symbol, and 'Type' in the rowData
head(rowData(sce.m1))

# Note that the rownames are the Ensembl gene names - let's switch to the more familiar gene symbols:
# ...but also notice that length(unique(rowData(sce.m1)$Symbol)) != nrow(sce.m1)
#   - That's because some gene symbols are used multiple times.

# -> To avoid getting into sticky situations, let's 'uniquify' these names:
rowData(sce.m1)$Symbol.uniq <- scuttle::uniquifyFeatureNames(ID = rowData(sce.m1)$ID, names = rowData(sce.m1)$Symbol)
rownames(sce.m1) <- rowData(sce.m1)$Symbol.uniq


## Read in Cell Ranger's t-SNE and add as an entry to the reducedDims slot, calling it 'TSNE_CR'
TSNE.coords <- read.csv(
    file = "/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/processed_data/cellranger/1M_C_LS/outs/analysis/tsne/2_components/projection.csv",
    header = TRUE, row.names = 1
)

head(TSNE.coords)
#                       TSNE.1     TSNE.2
# AAACCCAAGGTACATA-1 -15.54532 -33.666412
# AAACCCACATCCGAGC-1 -10.86300  11.466032
# AAACCCAGTAGCTCGC-1 -33.65646  14.447342
# AAACCCAGTTACGGAG-1  15.50977   1.180676
# AAACCCATCACGGACC-1 -16.26428  18.028536
# AAACCCATCCTCTGCA-1  15.79586  28.272553

# Check it's the same order as the colnames of the SCE
table(rownames(TSNE.coords) == colnames(sce.m1)) # all TRUE, good.
reducedDim(sce.m1, "TSNE_CR") <- TSNE.coords


## Do the same for the graph-based clustering - in this case as a column vector (class 'factor') in the colData:
graph.clust <- read.csv(
    file = "/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/processed_data/cellranger/1M_C_LS/outs/analysis/clustering/graphclust/clusters.csv",
    header = TRUE, row.names = 1
)

table(rownames(graph.clust) == colnames(sce.m1)) # all TRUE, good.
colData(sce.m1)$cluster.graph <- factor(graph.clust$Cluster)

# What's the distribution?
table(sce.m1$cluster.graph)
#    1    2    3    4    5    6    7    8    9   10
# 2034 1767 1197 1111  555  453  416  341  290  140


## Plot TSNE, coloring by graph-based clusters - this is the same as what's in the web summary :)
# plotReducedDim(sce.m1, dimred="TSNE_CR", colour_by="cluster.graph", text_by="cluster.graph")

# Can also print some violin plots of your favorite genes - however we first want to log-transform-normalize
#   the counts, due to such variable total counts/cell type(/nuclei)
sce.m1 <- logNormCounts(sce.m1, assay.type = "counts", log = TRUE, pseudo.count = 1)
# This creates & stores another assays entry, called "logcounts"
sce.m1

pdf("./graphics/1M_C_LS_CR-analysis_TSNE_and_violinPlots-std_graphClusters.pdf")
# TSNE:
plotReducedDim(sce.m1, dimred = "TSNE_CR", colour_by = "cluster.graph", text_by = "cluster.graph") +
    ggtitle("Cell-Ranger-automated analysis: TSNE (sample '1M_C_LS')")
# Violin plots of whatever genes you want to print
plotExpression(sce.m1,
    x = "cluster.graph", colour_by = "cluster.graph", exprs_values = "logcounts",
    features = c("Snap25", "Mbp", "Gfap", "Drd3", "Oxtr", "Nts"),
    ncol = 2, scales = "free_y"
)
dev.off()


# Using a neater 'plotExpressionCustom' Louise helped me create, for prettier graphics:
#   (italicizes the gene names and plots the median)
source("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/plotExpressionCustom.R")

pdf("./graphics/1M_C_LS_CR-analysis_violinPlots-custom_broadMarkers_graphClusters.pdf", width = 7, height = 10)
plotExpressionCustom(sce.m1,
    anno_name = "cluster.graph", exprs_values = "logcounts",
    # Broad cell types
    features = c(
        "Snap25", "Syt1", # neuronal
        "Slc17a7", "Slc7a6", # excit
        "Gad1", "Gad2", # inhib
        "Mbp", "Plp1", # oligo
        "Gfap", "Slc1a2", # astro
        "Pdgfra", "Vcan", # OPC
        "Cd74", "C3", # micro
        "Cldn5", "Flt1", # endothelial
        "Col1a2", "Tbx18"
    ), # mural
    features_name = "custom-selected",
    ncol = 2, scales = "free_y"
)
dev.off()

# From Keri/various RNAscope probes recently shown (human) data on
pdf("./graphics/1M_C_LS_CR-analysis_violinPlots-custom_specificMarkers_graphClusters.pdf")
plotExpressionCustom(sce.m1,
    anno_name = "cluster.graph", exprs_values = "logcounts",
    # Broad cell types
    features = c(
        "Drd3", "Oxtr", "Nts",
        "Gad1", "Ppp1r1b", "Slc17a6", "Slc32a1",
        "Calb1", "Pvalb", "Chat"
    ),
    features_name = "custom-selected",
    ncol = 3, scales = "free_y"
)
dev.off()


# For reference - what I usually look at, in human:
# markers.mathys.tran = list(
#   'neurons' = c('SYT1', 'SNAP25', 'GRIN1'),
#   'excitatory_neuron' = c('CAMK2A', 'NRGN','SLC17A7', 'SLC17A6', 'SLC17A8'),
#   'inhibitory_neuron' = c('GAD1', 'GAD2', 'SLC32A1'),
#   'oligodendrocyte' = c('MBP', 'MOBP', 'PLP1'),
#   'oligodendrocyte_precursor' = c('PDGFRA', 'VCAN', 'CSPG4'),
#   'microglia' = c('CD74', 'CSF1R', 'C3'),
#   'astrocytes' = c('GFAP', 'TNC', 'AQP4', 'SLC1A2'),
#   'endothelial' = c('CLDN5', 'FLT1', 'VTN'),
#   # MSN markers
#   'MSNs.pan' = c("PPP1R1B","BCL11B"),# "CTIP2")
#   'MSNs.D1' = c("DRD1", "PDYN", "TAC1"),
#   'MSNs.D2' = c("DRD2", "PENK"),
#   # Post-hoc from Tran-Maynard, et al. Neuron 2021
#   'differn_committed_OPC' = c("SOX4", "BCAN", "GPR17", "TNS3"),
#   'Tcell' = c('SKAP1', 'ITK', 'CD247'),
#   'Mural' = c('COL1A2', 'TBX18', 'RBPMS'),
#   'Macro' = c('CD163', 'SIGLEC1', 'F13A1')
# )







## Reproducibility information ====
print("Reproducibility information:")
Sys.time()
# [1] "2022-01-19 17:04:56 EST"
proc.time()
#     user    system   elapsed
#  277.934     8.779 14192.693
options(width = 120)
session_info()
# ─ Session info ───────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.1.2 Patched (2021-11-04 r81138)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2022-01-19
# pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
#
# ─ Packages ───────────────────────────────────────────────────────────────────
# package              * version  date (UTC) lib source
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
# beachmat               2.10.0   2021-10-26 [2] Bioconductor
# beeswarm               0.4.0    2021-06-01 [2] CRAN (R 4.1.2)
# Biobase              * 2.54.0   2021-10-26 [2] Bioconductor
# BiocGenerics         * 0.40.0   2021-10-26 [2] Bioconductor
# BiocIO                 1.4.0    2021-10-26 [2] Bioconductor
# BiocNeighbors          1.12.0   2021-10-26 [2] Bioconductor
# BiocParallel         * 1.28.3   2021-12-09 [2] Bioconductor
# BiocSingular           1.10.0   2021-10-26 [2] Bioconductor
# Biostrings             2.62.0   2021-10-26 [2] Bioconductor
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
# cli                    3.1.0    2021-10-27 [2] CRAN (R 4.1.2)
# colorspace             2.0-2    2021-06-24 [2] CRAN (R 4.1.0)
# cowplot                1.1.1    2020-12-30 [2] CRAN (R 4.1.2)
# crayon                 1.4.2    2021-10-29 [2] CRAN (R 4.1.2)
# DBI                    1.1.2    2021-12-20 [2] CRAN (R 4.1.2)
# DelayedArray           0.20.0   2021-10-26 [2] Bioconductor
# DelayedMatrixStats     1.16.0   2021-10-26 [2] Bioconductor
# digest                 0.6.29   2021-12-01 [2] CRAN (R 4.1.2)
# dplyr                  1.0.7    2021-06-18 [2] CRAN (R 4.1.0)
# dqrng                  0.3.0    2021-05-01 [2] CRAN (R 4.1.2)
# DropletUtils         * 1.14.2   2022-01-09 [2] Bioconductor
# edgeR                  3.36.0   2021-10-26 [2] Bioconductor
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
# fansi                  1.0.2    2022-01-14 [2] CRAN (R 4.1.2)
# farver                 2.1.0    2021-02-28 [2] CRAN (R 4.1.0)
# fs                     1.5.2    2021-12-08 [2] CRAN (R 4.1.2)
# gargle                 1.2.0    2021-07-02 [2] CRAN (R 4.1.0)
# generics               0.1.1    2021-10-25 [2] CRAN (R 4.1.2)
# GenomeInfoDb         * 1.30.0   2021-10-26 [2] Bioconductor
# GenomeInfoDbData       1.2.7    2021-11-01 [2] Bioconductor
# GenomicAlignments      1.30.0   2021-10-26 [2] Bioconductor
# GenomicRanges        * 1.46.1   2021-11-18 [2] Bioconductor
# ggbeeswarm             0.6.0    2017-08-07 [2] CRAN (R 4.1.2)
# ggplot2              * 3.3.5    2021-06-25 [2] CRAN (R 4.1.0)
# ggrepel                0.9.1    2021-01-15 [2] CRAN (R 4.1.0)
# glue                   1.6.0    2021-12-17 [2] CRAN (R 4.1.2)
# googledrive            2.0.0    2021-07-08 [2] CRAN (R 4.1.0)
# gridExtra              2.3      2017-09-09 [2] CRAN (R 4.1.0)
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
# HDF5Array              1.22.1   2021-11-14 [2] Bioconductor
# here                 * 1.0.1    2020-12-13 [2] CRAN (R 4.1.2)
# IRanges              * 2.28.0   2021-10-26 [2] Bioconductor
# irlba                  2.3.5    2021-12-06 [2] CRAN (R 4.1.2)
# jaffelab             * 0.99.31  2021-12-13 [1] Github (LieberInstitute/jaffelab@2cbd55a)
# labeling               0.4.2    2020-10-20 [2] CRAN (R 4.1.0)
# lattice                0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
# lifecycle              1.0.1    2021-09-24 [2] CRAN (R 4.1.2)
# limma                  3.50.0   2021-10-26 [2] Bioconductor
# locfit                 1.5-9.4  2020-03-25 [2] CRAN (R 4.1.0)
# magrittr               2.0.1    2020-11-17 [2] CRAN (R 4.1.0)
# Matrix                 1.4-0    2021-12-08 [3] CRAN (R 4.1.2)
# MatrixGenerics       * 1.6.0    2021-10-26 [2] Bioconductor
# matrixStats          * 0.61.0   2021-09-17 [2] CRAN (R 4.1.2)
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
# pillar                 1.6.4    2021-10-18 [2] CRAN (R 4.1.2)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# purrr                  0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# R.methodsS3            1.8.1    2020-08-26 [2] CRAN (R 4.1.0)
# R.oo                   1.24.0   2020-08-26 [2] CRAN (R 4.1.0)
# R.utils                2.11.0   2021-09-26 [2] CRAN (R 4.1.2)
# R6                     2.5.1    2021-08-19 [2] CRAN (R 4.1.2)
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.1.2)
# RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 4.1.0)
# Rcpp                   1.0.8    2022-01-13 [2] CRAN (R 4.1.2)
# RCurl                  1.98-1.5 2021-09-17 [2] CRAN (R 4.1.2)
# restfulr               0.0.13   2017-08-06 [2] CRAN (R 4.1.0)
# rhdf5                  2.38.0   2021-10-26 [2] Bioconductor
# rhdf5filters           1.6.0    2021-10-26 [2] Bioconductor
# Rhdf5lib               1.16.0   2021-10-26 [2] Bioconductor
# rjson                  0.2.21   2022-01-09 [2] CRAN (R 4.1.2)
# rlang                  0.4.12   2021-10-18 [2] CRAN (R 4.1.2)
# rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
# Rsamtools              2.10.0   2021-10-26 [2] Bioconductor
# rsvd                   1.0.5    2021-04-16 [2] CRAN (R 4.1.2)
# rtracklayer          * 1.54.0   2021-10-26 [2] Bioconductor
# S4Vectors            * 0.32.3   2021-11-21 [2] Bioconductor
# ScaledMatrix           1.2.0    2021-10-26 [2] Bioconductor
# scales                 1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
# scater               * 1.22.0   2021-10-26 [2] Bioconductor
# scuttle              * 1.4.0    2021-10-26 [2] Bioconductor
# segmented              1.3-4    2021-04-22 [1] CRAN (R 4.1.2)
# sessioninfo          * 1.2.2    2021-12-06 [2] CRAN (R 4.1.2)
# SingleCellExperiment * 1.16.0   2021-10-26 [2] Bioconductor
# sparseMatrixStats      1.6.0    2021-10-26 [2] Bioconductor
# SummarizedExperiment * 1.24.0   2021-10-26 [2] Bioconductor
# tibble                 3.1.6    2021-11-07 [2] CRAN (R 4.1.2)
# tidyselect             1.1.1    2021-04-30 [2] CRAN (R 4.1.0)
# utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
# vctrs                  0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
# vipor                  0.4.5    2017-03-22 [2] CRAN (R 4.1.2)
# viridis                0.6.2    2021-10-13 [2] CRAN (R 4.1.2)
# viridisLite            0.4.0    2021-04-13 [2] CRAN (R 4.1.0)
# withr                  2.4.3    2021-11-30 [2] CRAN (R 4.1.2)
# XML                    3.99-0.8 2021-09-17 [2] CRAN (R 4.1.2)
# XVector                0.34.0   2021-10-26 [2] Bioconductor
# yaml                   2.2.1    2020-02-01 [2] CRAN (R 4.1.0)
# zlibbioc               1.40.0   2021-10-26 [2] Bioconductor
#
# [1] /users/ntranngu/R/4.1.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
#
# ──────────────────────────────────────────────────────────────────────────────
