library("SingleCellExperiment")
library("here")
library("lobstr")
library("scuttle")
library("sessioninfo")

### For iSEE shiny app =============
load(
    here(
        "snRNAseq_mouse",
        "processed_data",
        "SCE",
        "sce_updated_LS.rda"
    ),
    verbose = TRUE
)

# As-is size:
lobstr::obj_size(sce.ls)
# 7.70 GB

# Remove the clusters we're not reporting
sce.ls <- sce.ls[, -grep("drop.", sce.ls$cellType.final)]
sce.ls <- sce.ls[, -grep("Neuron.mixed", sce.ls$cellType.final)]
sce.ls$cellType.final <- droplevels(sce.ls$cellType.final)
sce.ls$cellType.broad <- droplevels(sce.ls$cellType.broad)

# Corresponding hexadecimal colors
cell_colors.ls <-
    cell_colors.ls[-grep("drop", names(cell_colors.ls))]
cell_colors.ls <-
    cell_colors.ls[-grep("Neuron.mixed", names(cell_colors.ls))]

sce.hold <- sce.ls

# Remove 'counts' & 'binomial_pearson_residuals' (just keep 'logcounts')
assay(sce.ls, "counts") <- NULL
assay(sce.ls, "binomial_pearson_residuals") <- NULL
lobstr::obj_size(sce.ls)
# 878.46 MB

# 'GLMPCA_50' is just the top 50 PCs from the 100 computed in 'GLMPCA_approx'
reducedDim(sce.ls, "GLMPCA_approx") <- NULL
lobstr::obj_size(sce.ls)
# 860.78 MB

# Because rownames are currently in Ensembl ID, re-assign & 'uniquify' to gene symbols
rownames(sce.ls) <-
    uniquifyFeatureNames(rowData(sce.ls)$gene_id, rowData(sce.ls)$gene_name)

lobstr::obj_size(sce.ls)
# 860.78 MB

## Some add'l metrics for rowData ===
#    As created/implemented in https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/master/shiny_apps/00_clean_functions.R
rowData(sce.ls)$propNucleiExprs <- apply(
    assay(sce.hold, "counts"),
    1,
    function(x) {
        round(mean(x != 0), 3)
    }
)

# The above, by cell type ===
cellType.idx <- rafalib::splitit(sce.ls$cellType.final)
rowdat.sce <- rowData(sce.ls)
for (i in names(cellType.idx)) {
    message(Sys.time(), " computing propNucleiExprs for ", i)
    rowdat.sce[, paste0("propExprsIn.", i)] <- apply(
        assay(sce.hold, "counts")[, cellType.idx[[i]]],
        1,
        function(x) {
            round(mean(x != 0), 3)
        }
    )
}
rowData(sce.ls) <- rowdat.sce

## Clean up unnecessary info
keepCols <- setdiff(
    colnames(colData(sce.ls)),
    c(
        "high.mito",
        "sizeFactor",
        "clusters.glmpca",
        "cellType",
        "cellType.exp"
    )
)
colData(sce.ls) <- colData(sce.ls)[, keepCols]

lobstr::obj_size(sce.ls)
# 867.90 MB

# Re-name these
sce.ls.small <- sce.ls
cell_cols.clean <- cell_colors.ls

## Save this
save(
    sce.ls.small,
    cell_cols.clean,
    file = here(
        "snRNAseq_mouse",
        "code",
        "03_iSEE_apps",
        "snRNA-seq_app",
        "sce_for_iSEE_LS.rda"
    )
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

# ─ Session info ──────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.3.1 Patched (2023-06-20 r84596)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2023-06-22
#  pandoc   3.1.1 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.3/bin/pandoc
#
# ─ Packages ──────────────────────────────────────────────────────────────────────
#  package              * version   date (UTC) lib source
#  beachmat               2.16.0    2023-04-25 [2] Bioconductor
#  Biobase              * 2.60.0    2023-04-25 [2] Bioconductor
#  BiocGenerics         * 0.46.0    2023-04-25 [2] Bioconductor
#  BiocParallel           1.34.2    2023-05-22 [2] Bioconductor
#  bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.0)
#  cli                    3.6.1     2023-03-23 [2] CRAN (R 4.3.0)
#  codetools              0.2-19    2023-02-01 [3] CRAN (R 4.3.1)
#  colorout               1.2-2     2023-05-06 [1] Github (jalvesaq/colorout@79931fd)
#  colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.3.0)
#  crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.0)
#  DelayedArray           0.26.3    2023-05-22 [2] Bioconductor
#  DelayedMatrixStats     1.22.1    2023-06-09 [2] Bioconductor
#  digest                 0.6.31    2022-12-11 [2] CRAN (R 4.3.0)
#  dplyr                  1.1.2     2023-04-20 [2] CRAN (R 4.3.0)
#  fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.3.0)
#  fastmap                1.1.1     2023-02-24 [2] CRAN (R 4.3.0)
#  generics               0.1.3     2022-07-05 [2] CRAN (R 4.3.0)
#  GenomeInfoDb         * 1.36.0    2023-04-25 [2] Bioconductor
#  GenomeInfoDbData       1.2.10    2023-04-11 [2] Bioconductor
#  GenomicRanges        * 1.52.0    2023-04-25 [2] Bioconductor
#  ggplot2                3.4.2     2023-04-03 [2] CRAN (R 4.3.0)
#  glue                   1.6.2     2022-02-24 [2] CRAN (R 4.3.0)
#  gtable                 0.3.3     2023-03-21 [2] CRAN (R 4.3.0)
#  here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.3.0)
#  htmltools              0.5.5     2023-03-23 [2] CRAN (R 4.3.0)
#  htmlwidgets            1.6.2     2023-03-17 [2] CRAN (R 4.3.0)
#  httpuv                 1.6.11    2023-05-11 [2] CRAN (R 4.3.0)
#  IRanges              * 2.34.0    2023-04-25 [2] Bioconductor
#  jsonlite               1.8.5     2023-06-05 [2] CRAN (R 4.3.0)
#  later                  1.3.1     2023-05-02 [2] CRAN (R 4.3.0)
#  lattice                0.21-8    2023-04-05 [3] CRAN (R 4.3.1)
#  lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.3.0)
#  lobstr               * 1.1.2     2022-06-22 [2] CRAN (R 4.3.0)
#  magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.3.0)
#  Matrix                 1.5-4.1   2023-05-18 [3] CRAN (R 4.3.1)
#  MatrixGenerics       * 1.12.2    2023-06-09 [2] Bioconductor
#  matrixStats          * 1.0.0     2023-06-02 [2] CRAN (R 4.3.0)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 4.3.0)
#  pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.3.0)
#  pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.3.0)
#  png                    0.1-8     2022-11-29 [2] CRAN (R 4.3.0)
#  prettyunits            1.1.1     2020-01-24 [2] CRAN (R 4.3.0)
#  promises               1.2.0.1   2021-02-11 [2] CRAN (R 4.3.0)
#  R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.0)
#  rafalib                1.0.0     2015-08-09 [1] CRAN (R 4.3.0)
#  RColorBrewer           1.1-3     2022-04-03 [2] CRAN (R 4.3.0)
#  Rcpp                   1.0.10    2023-01-22 [2] CRAN (R 4.3.0)
#  RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.3.0)
#  rlang                  1.1.1     2023-04-28 [2] CRAN (R 4.3.0)
#  rmote                  0.3.4     2023-05-06 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.3.0)
#  S4Arrays               1.0.4     2023-05-14 [2] Bioconductor
#  S4Vectors            * 0.38.1    2023-05-02 [2] Bioconductor
#  scales                 1.2.1     2022-08-20 [2] CRAN (R 4.3.0)
#  scuttle              * 1.10.1    2023-05-02 [2] Bioconductor
#  servr                  0.27      2023-05-02 [1] CRAN (R 4.3.0)
#  sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.3.0)
#  SingleCellExperiment * 1.22.0    2023-04-25 [2] Bioconductor
#  sparseMatrixStats      1.12.1    2023-06-20 [2] Bioconductor
#  SummarizedExperiment * 1.30.2    2023-06-06 [2] Bioconductor
#  tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.3.0)
#  tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.3.0)
#  utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.3.0)
#  vctrs                  0.6.3     2023-06-14 [2] CRAN (R 4.3.1)
#  xfun                   0.39      2023-04-20 [2] CRAN (R 4.3.0)
#  XVector                0.40.0    2023-04-25 [2] Bioconductor
#  zlibbioc               1.46.0    2023-04-25 [2] Bioconductor
#
#  [1] /users/lcollado/R/4.3
#  [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.3/R/4.3/lib64/R/site-library
#  [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.3/R/4.3/lib64/R/library
#
# ─────────────────────────────────────────────────────────────────────────────────

