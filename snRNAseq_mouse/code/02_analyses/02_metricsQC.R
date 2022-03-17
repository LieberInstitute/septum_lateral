### LS snRNA-seq analysis
### Mito rate QC & doublet assessment
### Initiated: MNT/LAR/ST 25Mar2022

library(SingleCellExperiment)
library(DropletUtils)
library(scater)
library(scran)
library(rtracklayer)
library(BiocParallel)
library(jaffelab)
library(gridExtra)
library(here)
library(sessioninfo)

here()



### Mito rate QC ============

# Load nuclei-called SCE
load(here("snRNAseq_mouse","processed_data","SCE", "sce_working_LS.rda"), verbose=T)
    # sce.ls

sce.ls
    # class: SingleCellExperiment
    # dim: 32285 25527
    # metadata(1): Samples
    # assays(1): counts
    # rownames(32285): ENSMUSG00000051951 ENSMUSG00000089699 ...
    #   ENSMUSG00000095019 ENSMUSG00000095041
    # rowData names(6): source type ... gene_name gene_type
    # colnames(25527): 1_AAACCCAAGGTACATA-1 1_AAACCCACATCCGAGC-1 ...
    #   4_TTTGTTGCATACAGCT-1 4_TTTGTTGGTCAAACGG-1
    # colData names(4): Sample Barcode mouseNums Sex
    # reducedDimNames(0):
    # mainExpName: NULL
    # altExpNames(0):

head(colData(sce.ls))

# Compute mito rates at the sample-level
sample.idx <- splitit(sce.ls$Sample)

stats <- list()
for(i in names(sample.idx)){
  stats[[i]] <- perCellQCMetrics(sce.ls[ ,sample.idx[[i]]],
                                 subsets=list(Mito = which(seqnames(sce.ls) == "chrM")))
  }
names(stats) <- names(sample.idx)

save(sce.ls, stats,
     file=here("snRNAseq_mouse","processed_data","SCE", "sce_working_LS.rda"))


## 17Mar2022 work ===
load(here("snRNAseq_mouse","processed_data","SCE", "sce_working_LS.rda"), verbose=T)
    # sce.ls, stats




### Trick: Add a pseudo-count==1 for a 'MT transcript' ===
# Note: This was implemented because we realized samples with mito rate distributions that
#       were 'clean' and tightly distributed about 0 would yield a 3x MAD = 0, thus over-penalizing
#       nuclei even if they had a single MT transcript (throwing out upwards of 50% of the sample)

# First check computation of mito percent:
table(stats[[1]]$subsets_Mito_percent == (stats[[1]]$subsets_Mito_sum/stats[[1]]$sum)*100)
    # All TRUE

test.stats <- stats

for(i in 1:length(test.stats)){
  test.stats[[i]]$pseudo_subsets_Mito_sum <- test.stats[[i]]$subsets_Mito_sum + 1
  test.stats[[i]]$pseudo_subsets_Mito_percent <- test.stats[[i]]$pseudo_subsets_Mito_sum / (test.stats[[i]]$sum+1) * 100
}

## Lapply: MAD approach for mito rate thresholding
pseudo.high.mito <- lapply(test.stats, function(x) {
  isOutlier(x$pseudo_subsets_Mito_percent,
            nmads=3,
            type="higher")
  })
pseudo.high.mito.table <- lapply(pseudo.high.mito, table)

# Percent dropped
sapply(pseudo.high.mito.table, function(x) round(x[2]/sum(x), 3))
    #1M.TRUE 2F.TRUE 3M.TRUE 4F.TRUE 
    #0.096   0.073   0.128   0.107

# Thresholds
sapply(pseudo.high.mito, function(x){round(attributes(x)[["thresholds"]]["higher"], 4)})
    #1M.higher 2F.higher 3M.higher 4F.higher 
    #0.3476    0.2545    0.1644    0.1693


## Bind [true] stats to colData
# (we'll just keep the 'pseudo' result since this was made/meant to work with a range of data)
table(rownames(do.call("rbind", stats)) == colnames(sce.ls))
    # all 25527 TRUE

colData(sce.ls) <- cbind(colData(sce.ls),
                         do.call("rbind", stats),
                         do.call("c", pseudo.high.mito))
#colnames(colData(sce.ls))[which(colnames(colData(sce.ls)) == "do.call(\"c\", pseudo.high.mito)")] <- "high.mito"
colnames(colData(sce.ls))[11] <- "high.mito"

# $sum == $total
sce.ls$total <- NULL

# Store original for comparison/plotting
sce.ls.unfiltered <- sce.ls
sce.ls <- sce.ls[ ,!sce.ls$high.mito]




## Plot some metrics
mitoCutoffs <- sapply(pseudo.high.mito, function(x){round(attributes(x)[["thresholds"]]["higher"], 3)})
names(mitoCutoffs) <- gsub(".higher","", names(mitoCutoffs))

pdf(here("snRNAseq_mouse","plots","LS-n3_QCmetrics_high-mitoColored.pdf"), height=4)
for(i in names(sample.idx)){
  grid.arrange(
    plotColData(sce.ls.unfiltered[ ,sample.idx[[i]]], y="sum", colour_by="high.mito", point_alpha=0.4) +
      scale_y_log10() + ggtitle(paste0("Total count: ", i)),
    plotColData(sce.ls.unfiltered[ ,sample.idx[[i]]], y="detected", colour_by="high.mito", point_alpha=0.4) +
      scale_y_log10() + ggtitle("Detected features"),
    plotColData(sce.ls.unfiltered[ ,sample.idx[[i]]], y="subsets_Mito_percent",
                colour_by="high.mito", point_alpha=0.4) +
      ggtitle(paste0("Mito % (cutoff = ", mitoCutoffs[i],")")),
    ncol=3
  )
  # Mito rate vs n detected features
  print(
    plotColData(sce.ls.unfiltered[ ,sample.idx[[i]]], x="detected", y="subsets_Mito_percent",
                colour_by="high.mito", point_size=2.5, point_alpha=0.5) +
      ggtitle(paste0("Sample: ", i,
                     ";     pre-QC nNuclei: ", ncol(sce.ls.unfiltered[ ,sce.ls.unfiltered$Sample==i]),";      ",
                     "nNuclei kept: ", ncol(sce.ls[ ,sce.ls$Sample==i])," (",
                     round(ncol(sce.ls[ ,sce.ls$Sample==i]) /
                             ncol(sce.ls.unfiltered[ ,sce.ls.unfiltered$Sample==i]) * 100, 2), "%)"
      ))
  )
  # Detected features vs total count
  print(
    plotColData(sce.ls.unfiltered[ ,sample.idx[[i]]], x="sum", y="detected",
                colour_by="high.mito", point_size=2.5, point_alpha=0.5) +
      ggtitle(paste0("Sample: ", i,
                     ";     pre-QC nNuclei: ", ncol(sce.ls.unfiltered[ ,sce.ls.unfiltered$Sample==i]),";      ",
                     "nNuclei kept: ", ncol(sce.ls[ ,sce.ls$Sample==i])," (",
                     round(ncol(sce.ls[ ,sce.ls$Sample==i]) /
                             ncol(sce.ls.unfiltered[ ,sce.ls.unfiltered$Sample==i]) * 100, 2), "%)"
      ))
  )
}
dev.off()



## NEXT TIME =================
##
## compute doublet scores and append to colData 
##
                          ## ==================


# Save
save(sce.ls, sce.ls.unfiltered, 
     file=here("snRNAseq_mouse", "processed_data","SCE", "sce_working_LS.rda"))



## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
    #[1] "2022-03-17 14:34:24 EDT"
proc.time()
    #    user   system  elapsed 
    # 252.845   23.536 3534.600 
options(width = 120)
session_info()
# ─ Session info ────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.1.2 Patched (2021-11-04 r81138)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2022-03-17
# pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
# 
# ─ Packages ────────────────────────────────────────────────────────────────────
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
# bluster                1.4.0    2021-10-26 [2] Bioconductor
# cli                    3.2.0    2022-02-14 [2] CRAN (R 4.1.2)
# cluster                2.1.2    2021-04-17 [3] CRAN (R 4.1.2)
# colorspace             2.0-3    2022-02-21 [2] CRAN (R 4.1.2)
# cowplot                1.1.1    2020-12-30 [2] CRAN (R 4.1.2)
# crayon                 1.5.0    2022-02-14 [2] CRAN (R 4.1.2)
# DBI                    1.1.2    2021-12-20 [2] CRAN (R 4.1.2)
# DelayedArray           0.20.0   2021-10-26 [2] Bioconductor
# DelayedMatrixStats     1.16.0   2021-10-26 [2] Bioconductor
# digest                 0.6.29   2021-12-01 [2] CRAN (R 4.1.2)
# dplyr                  1.0.8    2022-02-08 [2] CRAN (R 4.1.2)
# dqrng                  0.3.0    2021-05-01 [2] CRAN (R 4.1.2)
# DropletUtils         * 1.14.2   2022-01-09 [2] Bioconductor
# edgeR                  3.36.0   2021-10-26 [2] Bioconductor
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
# fansi                  1.0.2    2022-01-14 [2] CRAN (R 4.1.2)
# farver                 2.1.0    2021-02-28 [2] CRAN (R 4.1.0)
# fs                     1.5.2    2021-12-08 [2] CRAN (R 4.1.2)
# gargle                 1.2.0    2021-07-02 [2] CRAN (R 4.1.0)
# generics               0.1.2    2022-01-31 [2] CRAN (R 4.1.2)
# GenomeInfoDb         * 1.30.1   2022-01-30 [2] Bioconductor
# GenomeInfoDbData       1.2.7    2021-11-01 [2] Bioconductor
# GenomicAlignments      1.30.0   2021-10-26 [2] Bioconductor
# GenomicRanges        * 1.46.1   2021-11-18 [2] Bioconductor
# ggbeeswarm             0.6.0    2017-08-07 [2] CRAN (R 4.1.2)
# ggplot2              * 3.3.5    2021-06-25 [2] CRAN (R 4.1.0)
# ggrepel                0.9.1    2021-01-15 [2] CRAN (R 4.1.0)
# glue                   1.6.2    2022-02-24 [2] CRAN (R 4.1.2)
# googledrive            2.0.0    2021-07-08 [2] CRAN (R 4.1.0)
# gridExtra            * 2.3      2017-09-09 [2] CRAN (R 4.1.0)
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
# HDF5Array              1.22.1   2021-11-14 [2] Bioconductor
# here                 * 1.0.1    2020-12-13 [2] CRAN (R 4.1.2)
# igraph                 1.2.11   2022-01-04 [2] CRAN (R 4.1.2)
# IRanges              * 2.28.0   2021-10-26 [2] Bioconductor
# irlba                  2.3.5    2021-12-06 [2] CRAN (R 4.1.2)
# jaffelab             * 0.99.31  2021-12-13 [1] Github (LieberInstitute/jaffelab@2cbd55a)
# labeling               0.4.2    2020-10-20 [2] CRAN (R 4.1.0)
# lattice                0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
# lifecycle              1.0.1    2021-09-24 [2] CRAN (R 4.1.2)
# limma                  3.50.1   2022-02-17 [2] Bioconductor
# locfit                 1.5-9.5  2022-03-03 [2] CRAN (R 4.1.2)
# magrittr               2.0.2    2022-01-26 [2] CRAN (R 4.1.2)
# Matrix                 1.4-0    2021-12-08 [3] CRAN (R 4.1.2)
# MatrixGenerics       * 1.6.0    2021-10-26 [2] Bioconductor
# matrixStats          * 0.61.0   2021-09-17 [2] CRAN (R 4.1.2)
# metapod                1.2.0    2021-10-26 [2] Bioconductor
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
# pillar                 1.7.0    2022-02-01 [2] CRAN (R 4.1.2)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# purrr                  0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# R.methodsS3            1.8.1    2020-08-26 [2] CRAN (R 4.1.0)
# R.oo                   1.24.0   2020-08-26 [2] CRAN (R 4.1.0)
# R.utils                2.11.0   2021-09-26 [2] CRAN (R 4.1.2)
# R6                     2.5.1    2021-08-19 [2] CRAN (R 4.1.2)
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.1.2)
# RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 4.1.0)
# Rcpp                   1.0.8.3  2022-03-17 [2] CRAN (R 4.1.2)
# RCurl                  1.98-1.6 2022-02-08 [2] CRAN (R 4.1.2)
# restfulr               0.0.13   2017-08-06 [2] CRAN (R 4.1.0)
# rhdf5                  2.38.1   2022-03-10 [2] Bioconductor
# rhdf5filters           1.6.0    2021-10-26 [2] Bioconductor
# Rhdf5lib               1.16.0   2021-10-26 [2] Bioconductor
# rjson                  0.2.21   2022-01-09 [2] CRAN (R 4.1.2)
# rlang                  1.0.2    2022-03-04 [2] CRAN (R 4.1.2)
# rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
# Rsamtools              2.10.0   2021-10-26 [2] Bioconductor
# rsvd                   1.0.5    2021-04-16 [2] CRAN (R 4.1.2)
# rtracklayer          * 1.54.0   2021-10-26 [2] Bioconductor
# S4Vectors            * 0.32.3   2021-11-21 [2] Bioconductor
# ScaledMatrix           1.2.0    2021-10-26 [2] Bioconductor
# scales                 1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
# scater               * 1.22.0   2021-10-26 [2] Bioconductor
# scran                * 1.22.1   2021-11-14 [2] Bioconductor
# scuttle              * 1.4.0    2021-10-26 [2] Bioconductor
# segmented              1.3-4    2021-04-22 [1] CRAN (R 4.1.2)
# sessioninfo          * 1.2.2    2021-12-06 [2] CRAN (R 4.1.2)
# SingleCellExperiment * 1.16.0   2021-10-26 [2] Bioconductor
# sparseMatrixStats      1.6.0    2021-10-26 [2] Bioconductor
# statmod                1.4.36   2021-05-10 [2] CRAN (R 4.1.0)
# SummarizedExperiment * 1.24.0   2021-10-26 [2] Bioconductor
# tibble                 3.1.6    2021-11-07 [2] CRAN (R 4.1.2)
# tidyselect             1.1.2    2022-02-21 [2] CRAN (R 4.1.2)
# utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
# vctrs                  0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
# vipor                  0.4.5    2017-03-22 [2] CRAN (R 4.1.2)
# viridis                0.6.2    2021-10-13 [2] CRAN (R 4.1.2)
# viridisLite            0.4.0    2021-04-13 [2] CRAN (R 4.1.0)
# withr                  2.5.0    2022-03-03 [2] CRAN (R 4.1.2)
# XML                    3.99-0.9 2022-02-24 [2] CRAN (R 4.1.2)
# XVector                0.34.0   2021-10-26 [2] Bioconductor
# yaml                   2.3.5    2022-02-21 [2] CRAN (R 4.1.2)
# zlibbioc               1.40.0   2021-10-26 [2] Bioconductor
# 
# [1] /users/ntranngu/R/4.1.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
# 
# ──────────────────────────────────────────────────────────────────

