### LC snRNA-seq analysis
### Build SCE from raw count matrices; perform nuclei calling
### qrsh -l bluejay,mf=52G,h_vmem=56G
### Adapted from https://github.com/lmweber/locus-c/blob/main/code/analyses_sn/01_buildSCE_callNuclei.R
### Initiated MNT 08Mar2022

library(SingleCellExperiment)
library(DropletUtils)
library(rtracklayer)
library(BiocParallel)
library(jaffelab)
library(here)
library(sessioninfo)

here()
# [1] "/dcs04/lieber/marmaypag/pilotLS_LIBD1070"

### Create raw SCE ====
# Basic sample info
sample.info <- read.table(here(
    "snRNAseq_mouse", "raw_data", "sample_info",
    "sample_info_simplified.tsv"
))

sample.info$path <- file.path(
    here("snRNAseq_mouse", "processed_data", "cellranger"),
    sample.info$V1,
    "outs",
    "raw_feature_bc_matrix"
)
stopifnot(all(file.exists(sample.info$path)))

## Build basic SCE (will add more subject-level colData later, once obtained)
Sys.time()
# [1] "2022-03-10 07:53:33 EST"
sce <- read10xCounts(
    samples = sample.info$path,
    sample.names = ss(sample.info$V1, "_", 1),
    type = "sparse",
    col.names = TRUE
)
Sys.time()
# [1] "2022-03-10 07:56:43 EST"


## Read in the gene information from the annotation GTF file
# (following Leo's method in https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/build_basic_sce.R)
gtf <-
    rtracklayer::import(
        "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-mm10-2020-A/genes/genes.gtf"
    )
gtf <- gtf[gtf$type == "gene"]
names(gtf) <- gtf$gene_id

## Match the genes
match_genes <- match(rownames(sce), gtf$gene_id)
stopifnot(all(!is.na(match_genes)))

## Keep only some columns from the gtf
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_version", "gene_name", "gene_type")]

## Add the gene info to our SCE object
rowRanges(sce) <- gtf[match_genes]

## Inspect object
sce
# class: SingleCellExperiment
# dim: 32285 9521311
# metadata(1): Samples
# assays(1): counts
# rownames(32285): ENSMUSG00000051951 ENSMUSG00000089699 ...
# ENSMUSG00000095019 ENSMUSG00000095041
# rowData names(6): source type ... gene_name gene_type
# colnames(9521311): 1_AAACCCAAGAAACACT-1 1_AAACCCAAGAAACCAT-1 ...
# 4_TTTGTTGTCTTTGATC-1 4_TTTGTTGTCTTTGCGC-1
# colData names(2): Sample Barcode
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):

table(sce$Sample)
#     1M      2F      3M      4F
# 2310153 2555720 2382021 2273417

# Add mouse IDs too and sex
sce$mouseNums <- sample.info$V3[match(
    sce$Sample,
    ss(sample.info$V1, "_", 1)
)]
sce$Sex <- substr(sce$Sample, 2, 2)

sce.ls <- sce
sample.idx <- splitit(sce.ls$Sample)

## First use barcodeRanks to ID the '2nd knee' and make some plots to diagnose ===
# dir.create(here("snRNAseq_mouse","plots"))
bcRanks.ls <- list()

Sys.time()
# [1] "2022-03-10 08:09:53 EST"
for (i in names(sample.idx)) {

    # Re-define the `lower=` param by identifying the 'second knee point'
    bcRanks.ls[[i]] <- barcodeRanks(counts(sce.ls[, sample.idx[[i]]]),
        fit.bounds = c(10, 1e3)
    )
    cat(paste0(
        "'Second knee point' for ", i, " is: ",
        metadata(bcRanks.ls[[i]])$knee, "\n"
    ))

    # Plot the barcode rank plot with +100 and +300
    png(here(
        "snRNAseq_mouse", "plots",
        paste0("barcodeRankPlot_", i, "_w-fitbounds_k2plus100-300UMIs.png")
    ))
    plot(bcRanks.ls[[i]]$rank, bcRanks.ls[[i]]$total,
        log = "xy", xlab = "Barcode Rank", ylab = "Total UMIs",
        main = paste0("Barcode rank plot for: ", i, "\n( `fit.bounds=c(10,1e3)` )"),
        cex.axis = 0.8, cex.lab = 1.2, las = 1
    )
    o <- order(bcRanks.ls[[i]]$rank)
    lines(bcRanks.ls[[i]]$rank[o], bcRanks.ls[[i]]$fitted[o], col = "red")
    # k_2 from above
    abline(h = metadata(bcRanks.ls[[i]])$knee, col = "darkblue", lty = 2)
    # k_2 + 100:
    abline(h = metadata(bcRanks.ls[[i]])$knee + 100, col = "dodgerblue", lty = 2)
    # k_2 + 300:
    abline(h = metadata(bcRanks.ls[[i]])$knee + 300, col = "darkgreen", lty = 2)
    legend("bottomleft",
        lty = 2, col = c("darkblue", "dodgerblue", "darkgreen"),
        legend = c("knee_2", "knee_2+100", "knee_2+300")
    )
    dev.off()
}
Sys.time()
# 'Second knee point' for 1M is: 430
# 'Second knee point' for 2F is: 560
# 'Second knee point' for 3M is: 583
# 'Second knee point' for 4F is: 399

# Diagnosis: +300 UMIs to the '2nd knee point' will be better than +100 that sufficed
#            for the LC project.  However, it would be better to calculate the 2nd deriv
#            minimum to adaptively identify where this '2nd plateau' starts




## NOW run emptyDrops ===
e.out.custom <- list()

Sys.time()
# [1] "2022-03-10 08:46:33 EST"
for (i in names(sample.idx)) {
    cat(paste0(
        "'Second knee point' for ", i, " is: ",
        metadata(bcRanks.ls[[i]])$knee, "\n"
    ))
    cat(paste0(
        "\tSetting `lower` for `empyDrops()` to = ",
        metadata(bcRanks.ls[[i]])$knee + 300, "\n\n"
    ))

    # Set `lower = bcRanks.ls[[i]] + 300` (to capture the 'plateau' mode of low UMI totals)
    cat(paste0("\tSimulating empty drops for: ", i, "... \n"))
    set.seed(109)
    e.out.custom[[i]] <- emptyDrops(counts(sce.ls[, sample.idx[[i]]]),
        niters = 20000,
        lower = metadata(bcRanks.ls[[i]])$knee + 300
    )
    cat(paste0("\n\t...Simulations complete. \n\t", Sys.time(), "\n\n\n"))
}
#


## emptyDrops() results ===
for (i in 1:length(e.out.custom)) {
    print(names(e.out.custom)[[i]])
    print(table(
        Signif = e.out.custom[[i]]$FDR <= 0.001,
        Limited = e.out.custom[[i]]$Limited
    ))
}

# [1] "1M"
#       Limited
# Signif  FALSE TRUE
#   FALSE  7918    0
#   TRUE     67 6661

# [1] "2F"
#       Limited
# Signif  FALSE  TRUE
#   FALSE 23312     0
#   TRUE   1033 18588   - this is still too high...

# [1] "3M"
#       Limited
# Signif  FALSE TRUE
#   FALSE  6172    0
#   TRUE     54 6907

# [1] "4F"
#       Limited
# Signif  FALSE TRUE
#   FALSE  3128    0
#   TRUE     70 7786

# For reference:
sample.info
#        V1             V2            V3      V4
# 1 1M_C_LS Lateral_Septum ms5478_ms5480 1M-C-LS
# 2 2F_C_LS Lateral_Septum ms5483_ms5484 2F-C-LS
# 3 3M_C_LS Lateral_Septum ms5479_ms5481 3M-C-LS
# 4 4F_C_LS Lateral_Septum ms5485_ms5487 4F-C-LS


## Save this with the raw SCE for interactive downstream analyses ===
# dir.create(here("snRNAseq_mouse","processed_data","SCE"))
README.custom <- "Object 'e.out.custom' is the output of emptyDrops() with lower= set to the quantified
  'second knee point' (+300) to better model empty/debris-containing droplets"
save(sce.ls,
    e.out.custom, README.custom,
    file = here("snRNAseq_mouse", "processed_data", "SCE", "sce_raw_LS.rda")
)


## Btw without a 'lower=' param for barcodeRanks, it'll most likely pick the
#    '2nd' inflection point
bcRanks.2F <- barcodeRanks(counts(sce.ls[, sample.idx[["2F"]]]),
    fit.bounds = c(1e3, 1e5),
    lower = metadata(bcRanks.ls[[i]])$knee + 300
)

metadata(bcRanks.2F)
#      knee inflection
#     10524       6491

# Re-plot this with that inflection point
i <- "2F"
png(here(
    "snRNAseq_mouse", "plots",
    paste0("barcodeRankPlot_2F_w-fitbounds_k2plus100-300UMIs.png")
))
plot(bcRanks.ls[[i]]$rank, bcRanks.ls[[i]]$total,
    log = "xy", xlab = "Barcode Rank", ylab = "Total UMIs",
    main = paste0("Barcode rank plot for: ", i, "\n( `fit.bounds=c(10,1e3)` )"),
    cex.axis = 0.8, cex.lab = 1.2, las = 1
)
o <- order(bcRanks.ls[[i]]$rank)
lines(bcRanks.ls[[i]]$rank[o], bcRanks.ls[[i]]$fitted[o], col = "red")
# k_2 from above
abline(h = metadata(bcRanks.ls[[i]])$knee, col = "darkblue", lty = 2)
# k_2 + 100:
abline(h = metadata(bcRanks.ls[[i]])$knee + 100, col = "dodgerblue", lty = 2)
# k_2 + 300:
abline(h = metadata(bcRanks.ls[[i]])$knee + 300, col = "darkgreen", lty = 2)
# Re-computed inflection:
abline(h = metadata(bcRanks.2F)$inflection, col = "darkred", lty = 2)
legend("bottomleft",
    lty = 2, col = c("darkblue", "dodgerblue", "darkgreen", "darkred"),
    legend = c("knee_2", "knee_2+100", "knee_2+300", "inflection")
)
dev.off()


# What about restricting to this and having a high scores emptyDrops?
table(rownames(bcRanks.2F) == rownames(e.out.custom[["2F"]])) # all good

table(bcRanks.2F$total >= metadata(bcRanks.2F)$inflection &
    e.out.custom[["2F"]]$FDR <= 0.001)
#  FALSE    TRUE
# 2551738    3982




## Filter based on this revised nuclei-calling approach and save into another .rda ===
BCs.2F <- rownames(e.out.custom[["2F"]])[bcRanks.2F$total >= metadata(bcRanks.2F)$inflection &
    e.out.custom[["2F"]]$FDR <= 0.001]
e.out.rest <- e.out.custom
e.out.rest[["2F"]] <- NULL

BCs.rest <- unlist(
    lapply(e.out.rest, function(x) {
        rownames(x)[which(x$FDR <= 0.001)]
    })
)

BCs.keep <- c(BCs.2F, BCs.rest)
sce.ls <- sce.ls[, colnames(sce.ls) %in% BCs.keep]

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

save(sce.ls, file = here("snRNAseq_mouse", "processed_data", "SCE", "sce_working_LS.rda"))


## Reproducibility information ====
print("Reproducibility information:")
Sys.time()
# [1] "2022-03-14 17:52:23 EDT"
proc.time()
#    user   system  elapsed
# 756.542   48.333 6516.563
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
# date     2022-03-14
# pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
#
# ─ Packages ────────────────────────────────────────────────────────────────────
# package              * version  date (UTC) lib source
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
# beachmat               2.10.0   2021-10-26 [2] Bioconductor
# Biobase              * 2.54.0   2021-10-26 [2] Bioconductor
# BiocGenerics         * 0.40.0   2021-10-26 [2] Bioconductor
# BiocIO                 1.4.0    2021-10-26 [2] Bioconductor
# BiocParallel         * 1.28.3   2021-12-09 [2] Bioconductor
# Biostrings             2.62.0   2021-10-26 [2] Bioconductor
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
# cli                    3.2.0    2022-02-14 [2] CRAN (R 4.1.2)
# crayon                 1.5.0    2022-02-14 [2] CRAN (R 4.1.2)
# DBI                    1.1.2    2021-12-20 [2] CRAN (R 4.1.2)
# DelayedArray           0.20.0   2021-10-26 [2] Bioconductor
# DelayedMatrixStats     1.16.0   2021-10-26 [2] Bioconductor
# dplyr                  1.0.8    2022-02-08 [2] CRAN (R 4.1.2)
# dqrng                  0.3.0    2021-05-01 [2] CRAN (R 4.1.2)
# DropletUtils         * 1.14.2   2022-01-09 [2] Bioconductor
# edgeR                  3.36.0   2021-10-26 [2] Bioconductor
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
# fansi                  1.0.2    2022-01-14 [2] CRAN (R 4.1.2)
# fs                     1.5.2    2021-12-08 [2] CRAN (R 4.1.2)
# gargle                 1.2.0    2021-07-02 [2] CRAN (R 4.1.0)
# generics               0.1.2    2022-01-31 [2] CRAN (R 4.1.2)
# GenomeInfoDb         * 1.30.1   2022-01-30 [2] Bioconductor
# GenomeInfoDbData       1.2.7    2021-11-01 [2] Bioconductor
# GenomicAlignments      1.30.0   2021-10-26 [2] Bioconductor
# GenomicRanges        * 1.46.1   2021-11-18 [2] Bioconductor
# glue                   1.6.2    2022-02-24 [2] CRAN (R 4.1.2)
# googledrive            2.0.0    2021-07-08 [2] CRAN (R 4.1.0)
# HDF5Array              1.22.1   2021-11-14 [2] Bioconductor
# here                 * 1.0.1    2020-12-13 [2] CRAN (R 4.1.2)
# IRanges              * 2.28.0   2021-10-26 [2] Bioconductor
# jaffelab             * 0.99.31  2021-12-13 [1] Github (LieberInstitute/jaffelab@2cbd55a)
# lattice                0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
# lifecycle              1.0.1    2021-09-24 [2] CRAN (R 4.1.2)
# limma                  3.50.1   2022-02-17 [2] Bioconductor
# locfit                 1.5-9.5  2022-03-03 [2] CRAN (R 4.1.2)
# magrittr               2.0.2    2022-01-26 [2] CRAN (R 4.1.2)
# Matrix                 1.4-0    2021-12-08 [3] CRAN (R 4.1.2)
# MatrixGenerics       * 1.6.0    2021-10-26 [2] Bioconductor
# matrixStats          * 0.61.0   2021-09-17 [2] CRAN (R 4.1.2)
# pillar                 1.7.0    2022-02-01 [2] CRAN (R 4.1.2)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# purrr                  0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# R.methodsS3            1.8.1    2020-08-26 [2] CRAN (R 4.1.0)
# R.oo                   1.24.0   2020-08-26 [2] CRAN (R 4.1.0)
# R.utils                2.11.0   2021-09-26 [2] CRAN (R 4.1.2)
# R6                     2.5.1    2021-08-19 [2] CRAN (R 4.1.2)
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.1.2)
# RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 4.1.0)
# Rcpp                   1.0.8.2  2022-03-11 [2] CRAN (R 4.1.2)
# RCurl                  1.98-1.6 2022-02-08 [2] CRAN (R 4.1.2)
# restfulr               0.0.13   2017-08-06 [2] CRAN (R 4.1.0)
# rhdf5                  2.38.1   2022-03-10 [2] Bioconductor
# rhdf5filters           1.6.0    2021-10-26 [2] Bioconductor
# Rhdf5lib               1.16.0   2021-10-26 [2] Bioconductor
# rjson                  0.2.21   2022-01-09 [2] CRAN (R 4.1.2)
# rlang                  1.0.2    2022-03-04 [2] CRAN (R 4.1.2)
# rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
# Rsamtools              2.10.0   2021-10-26 [2] Bioconductor
# rtracklayer          * 1.54.0   2021-10-26 [2] Bioconductor
# S4Vectors            * 0.32.3   2021-11-21 [2] Bioconductor
# scuttle                1.4.0    2021-10-26 [2] Bioconductor
# segmented              1.3-4    2021-04-22 [1] CRAN (R 4.1.2)
# sessioninfo          * 1.2.2    2021-12-06 [2] CRAN (R 4.1.2)
# SingleCellExperiment * 1.16.0   2021-10-26 [2] Bioconductor
# sparseMatrixStats      1.6.0    2021-10-26 [2] Bioconductor
# SummarizedExperiment * 1.24.0   2021-10-26 [2] Bioconductor
# tibble                 3.1.6    2021-11-07 [2] CRAN (R 4.1.2)
# tidyselect             1.1.2    2022-02-21 [2] CRAN (R 4.1.2)
# utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
# vctrs                  0.3.8    2021-04-29 [2] CRAN (R 4.1.0)
# XML                    3.99-0.9 2022-02-24 [2] CRAN (R 4.1.2)
# XVector                0.34.0   2021-10-26 [2] Bioconductor
# yaml                   2.3.5    2022-02-21 [2] CRAN (R 4.1.2)
# zlibbioc               1.40.0   2021-10-26 [2] Bioconductor
#
# [1] /users/ntranngu/R/4.1.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
#
# ───────────────────────────────────────────────────────────────────────────────
