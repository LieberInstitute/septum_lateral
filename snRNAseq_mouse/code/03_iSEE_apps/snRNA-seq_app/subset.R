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
# 860.78 MB

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

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 4.3.2 (2023-10-31)
#  os       macOS Ventura 13.6
#  system   x86_64, darwin20
#  ui       RStudio
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       America/New_York
#  date     2023-11-17
#  rstudio  2023.09.1+494 Desert Sunflower (desktop)
#  pandoc   NA
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date (UTC) lib source
#  abind                  1.4-5     2016-07-21 [1] CRAN (R 4.3.0)
#  beachmat               2.16.0    2023-05-11 [1] Bioconductor
#  Biobase              * 2.60.0    2023-05-11 [1] Bioconductor
#  BiocGenerics         * 0.46.0    2023-05-11 [1] Bioconductor
#  BiocParallel           1.34.2    2023-05-22 [1] Bioconductor
#  bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.0)
#  cli                    3.6.1     2023-03-23 [1] CRAN (R 4.3.0)
#  codetools              0.2-19    2023-02-01 [1] CRAN (R 4.3.2)
#  crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
#  DelayedArray           0.26.7    2023-07-28 [1] Bioconductor
#  DelayedMatrixStats     1.22.6    2023-08-28 [1] Bioconductor
#  fansi                  1.0.5     2023-10-08 [1] CRAN (R 4.3.0)
#  GenomeInfoDb         * 1.36.4    2023-10-02 [1] Bioconductor
#  GenomeInfoDbData       1.2.10    2023-05-22 [1] Bioconductor
#  GenomicRanges        * 1.52.1    2023-10-08 [1] Bioconductor
#  glue                   1.6.2     2022-02-24 [1] CRAN (R 4.3.0)
#  here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.3.0)
#  IRanges              * 2.34.1    2023-06-22 [1] Bioconductor
#  lattice                0.21-9    2023-10-01 [1] CRAN (R 4.3.2)
#  lifecycle              1.0.3     2022-10-07 [1] CRAN (R 4.3.0)
#  lobstr               * 1.1.2     2022-06-22 [1] CRAN (R 4.3.0)
#  magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
#  Matrix                 1.6-1.1   2023-09-18 [1] CRAN (R 4.3.2)
#  MatrixGenerics       * 1.12.3    2023-07-31 [1] Bioconductor
#  matrixStats          * 1.0.0     2023-06-02 [1] CRAN (R 4.3.0)
#  paletteer              1.5.0     2022-10-19 [1] CRAN (R 4.3.0)
#  pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
#  pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
#  prettyunits            1.2.0     2023-09-24 [1] CRAN (R 4.3.0)
#  rafalib                1.0.0     2015-08-09 [1] CRAN (R 4.3.0)
#  RColorBrewer           1.1-3     2022-04-03 [1] CRAN (R 4.3.0)
#  Rcpp                   1.0.11    2023-07-06 [1] CRAN (R 4.3.0)
#  RCurl                  1.98-1.12 2023-03-27 [1] CRAN (R 4.3.0)
#  rematch2               2.1.2     2020-05-01 [1] CRAN (R 4.3.0)
#  rlang                  1.1.1     2023-04-28 [1] CRAN (R 4.3.0)
#  rprojroot              2.0.3     2022-04-02 [1] CRAN (R 4.3.0)
#  rstudioapi             0.15.0    2023-07-07 [1] CRAN (R 4.3.0)
#  S4Arrays               1.0.6     2023-08-30 [1] Bioconductor
#  S4Vectors            * 0.38.2    2023-09-22 [1] Bioconductor
#  scuttle              * 1.10.3    2023-10-10 [1] Bioconductor
#  sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.0)
#  SingleCellExperiment * 1.22.0    2023-05-11 [1] Bioconductor
#  sparseMatrixStats      1.12.2    2023-07-02 [1] Bioconductor
#  SummarizedExperiment * 1.30.2    2023-06-06 [1] Bioconductor
#  tibble                 3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
#  utf8                   1.2.3     2023-01-31 [1] CRAN (R 4.3.0)
#  vctrs                  0.6.4     2023-10-12 [1] CRAN (R 4.3.0)
#  XVector                0.40.0    2023-05-11 [1] Bioconductor
#  zlibbioc               1.46.0    2023-05-11 [1] Bioconductor
#
#  [1] /Library/Frameworks/R.framework/Versions/4.3-x86_64/Resources/library
#
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
