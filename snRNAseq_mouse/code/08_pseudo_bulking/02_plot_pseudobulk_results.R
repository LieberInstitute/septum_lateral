library(SingleCellExperiment)
library(sessioninfo)
library(scater)
library(here)

################### Load sce object and pseudobulk results ####################

load(here("snRNAseq_mouse",
          "processed_data",
          "SCE",
          "sce_updated_LS.rda"),
    verbose = TRUE)
# Loading objects:
#     sce.ls
#     annotationTab.ls
#     cell_colors.ls

load(here("snRNAseq_mouse",
          "processed_data",
          "SCE",
          "sce_pseudobulking_LS_and_broad_110623.rda"),
    verbose = TRUE)
# Loading objects:
#     sce_pseudo_all
#     sce_pseudo_LS
#     sce_pseudo_LS.Sept
#     sce_pseudo_neuronal
#     sce_pseudo_neuronal.broad

###############################################################################



############################## Setting up colors ##############################
cell_colors_all <- c(1:length(levels(colData(sce_pseudo_all)$cellType.broad)))
names(cell_colors_all) <- levels(colData(sce_pseudo_all)$cellType.broad)
neuronal_names <- c("Chol", "LS", "Sept", "Str", "MS", "TNoS", "TT.IG.SH", "Thal", "IoC")
cell_colors_all[names(cell_colors_all) %in% neuronal_names] <- "#e68f00"
cell_colors_all[!names(cell_colors_all) %in% neuronal_names] <- "#006164"

cell_colors_LS <- cell_colors.ls[grep("LS",names(cell_colors.ls))]
cell_colors_LS[1:length(cell_colors_LS)] <- c("#D62728","#FF9896","#C5B0D5",
                                              "#8C564B", "#C49C94", "#E377C2",
                                              "#F7B6D2", "#7F7F7F")

cell_colors_LS.Sept <- cell_colors.ls[grep("LS|Sept",names(cell_colors.ls))]
cell_colors_LS.Sept[1:length(cell_colors_LS.Sept)] <- c("#D62728","#FF9896","#C5B0D5",
                                                        "#8C564B", "#C49C94", "#E377C2",
                                                        "#F7B6D2", "#7F7F7F", "#ED665D",
                                                        "#AD8BC9")

cell_colors_neur <- cell_colors.ls[grep("LS|Sept|MS|TT|TNoS",names(cell_colors.ls))]
cell_colors_neur[1:length(cell_colors_neur)] <- c("#D62728","#FF9896","#C5B0D5",
                                                  "#8C564B", "#C49C94", "#E377C2",
                                                  "#F7B6D2", "#7F7F7F", "#BCBD22",
                                                  "#DBDB8D", "#ED665D","#AD8BC9",
                                                  "#6DCCDA", "#999999","#E69F00",
                                                  "#56B4E9")

###############################################################################



############################ Plot LS and broad PCs ############################

colData(sce_pseudo_all)$cell.group <- "non neuronal"
colData(sce_pseudo_all)$cell.group[colData(sce_pseudo_all)$cellType.broad %in% neuronal_names] <- "neuronal"

## Broad clusters
p_all <- plotPCA(
    sce_pseudo_all,
    colour_by = "cellType.broad",
    shape_by = "cell.group",
    ncomponents = 4,
    point_size = 2,
    label_format = c("%s %02i", " (%i%%)"),
    percentVar = (metadata(sce_pseudo_all)$PCA_var_explained) * 100
    ) + scale_color_manual("cellType.broad", values = cell_colors_all)

pdf(file = here("snRNAseq_mouse/plots/08_pseudo_bulking/sce_pseudo_broad_PCs_110623.pdf"), width = 8, height = 8)
print(p_all)
dev.off()

## Broad clusters colored by cell.group
p_all <- plotPCA(
    sce_pseudo_all,
    colour_by = "cell.group",
    ncomponents = 4,
    point_size = 2,
    label_format = c("%s %02i", " (%i%%)"),
    percentVar = (metadata(sce_pseudo_all)$PCA_var_explained) * 100
    ) + scale_color_manual("cellType.broad", values = cell_colors_all)

pdf(file = here("snRNAseq_mouse/plots/08_pseudo_bulking/sce_pseudo_broad_PCs_v2_110623.pdf"), width = 8, height = 8)
print(p_all)
dev.off()

## LS clusters
p_LS <- plotPCA(
    sce_pseudo_LS,
    colour_by = "cellType.final",
    ncomponents = 4,
    point_size = 2,
    label_format = c("%s %02i", " (%i%%)"),
    percentVar = (metadata(sce_pseudo_LS)$PCA_var_explained) * 100
    ) + scale_color_manual("cellType.final",values = cell_colors_LS)

pdf(file = here("snRNAseq_mouse/plots/08_pseudo_bulking/sce_pseudo_LS_PCs_110623.pdf"), width = 8, height = 8)
print(p_LS)
dev.off()

## LS and Sept clusters
p_LS_Sept <- plotPCA(
    sce_pseudo_LS.Sept,
    colour_by = "cellType.final",
    ncomponents = 5,
    point_size = 2,
    label_format = c("%s %02i", " (%i%%)"),
    percentVar = (metadata(sce_pseudo_LS.Sept)$PCA_var_explained) * 100
    ) + scale_color_manual("cellType.final",values = cell_colors_LS.Sept)

pdf(file = here("snRNAseq_mouse/plots/08_pseudo_bulking/sce_pseudo_LS_Sept_PCs_110623.pdf"), width = 10, height = 10)
print(p_LS_Sept)
dev.off()

p_neur <- plotPCA(
    sce_pseudo_neuronal,
    colour_by = "cellType.final",
    ncomponents = 5,
    point_size = 2,
    label_format = c("%s %02i", " (%i%%)"),
    percentVar = (metadata(sce_pseudo_neuronal)$PCA_var_explained) * 100
    ) + scale_color_manual("cellType.final", values = cell_colors_neur)

pdf(file = here("snRNAseq_mouse/plots/08_pseudo_bulking/sce_pseudo_neuronal_PCs_110623.pdf"), width = 10, height = 10)
print(p_neur)
dev.off()

p_neur <- plotPCA(
    sce_pseudo_neuronal.broad,
    colour_by = "cellType.broad",
    ncomponents = 12,
    point_size = 1,
    label_format = c("%s %02i", " (%i%%)"),
    percentVar = (metadata(sce_pseudo_neuronal.broad)$PCA_var_explained) * 100
    )
#+ scale_color_manual("cellType.broad", values = cell_colors_all)

pdf(file = here("snRNAseq_mouse/plots/08_pseudo_bulking/sce_pseudo_neuronal_broad_PCs_110623.pdf"), width = 14, height = 14)
print(p_neur)
dev.off()


###############################################################################
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
# [1] "Reproducibility information:"
# [1] "2023-11-06 17:10:37 EST"
# user  system elapsed
# 46.450   2.166 914.126
# ─ Session info ─────────────────────────────
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
# ─ Packages ─────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [2] CRAN (R 4.3.1)
# beachmat               2.16.0    2023-04-25 [2] Bioconductor
# beeswarm               0.4.0     2021-06-01 [2] CRAN (R 4.3.1)
# Biobase              * 2.60.0    2023-04-25 [2] Bioconductor
# BiocGenerics         * 0.46.0    2023-04-25 [2] Bioconductor
# BiocNeighbors          1.18.0    2023-04-25 [2] Bioconductor
# BiocParallel           1.34.2    2023-05-22 [2] Bioconductor
# BiocSingular           1.16.0    2023-04-25 [2] Bioconductor
# bitops                 1.0-7     2021-04-24 [2] CRAN (R 4.3.1)
# cli                    3.6.1     2023-03-23 [2] CRAN (R 4.3.1)
# codetools              0.2-19    2023-02-01 [3] CRAN (R 4.3.1)
# colorout             * 1.2-2     2023-09-22 [1] Github (jalvesaq/colorout@79931fd)
# colorspace             2.1-0     2023-01-23 [2] CRAN (R 4.3.1)
# crayon                 1.5.2     2022-09-29 [2] CRAN (R 4.3.1)
# DelayedArray           0.26.7    2023-07-28 [2] Bioconductor
# DelayedMatrixStats     1.22.6    2023-08-28 [2] Bioconductor
# dplyr                  1.1.3     2023-09-03 [2] CRAN (R 4.3.1)
# fansi                  1.0.4     2023-01-22 [2] CRAN (R 4.3.1)
# farver                 2.1.1     2022-07-06 [2] CRAN (R 4.3.1)
# generics               0.1.3     2022-07-05 [2] CRAN (R 4.3.1)
# GenomeInfoDb         * 1.36.3    2023-09-07 [2] Bioconductor
# GenomeInfoDbData       1.2.10    2023-07-20 [2] Bioconductor
# GenomicRanges        * 1.52.0    2023-04-25 [2] Bioconductor
# ggbeeswarm             0.7.2     2023-04-29 [2] CRAN (R 4.3.1)
# ggplot2              * 3.4.3     2023-08-14 [2] CRAN (R 4.3.1)
# ggrepel                0.9.3     2023-02-03 [2] CRAN (R 4.3.1)
# glue                   1.6.2     2022-02-24 [2] CRAN (R 4.3.1)
# gridExtra              2.3       2017-09-09 [2] CRAN (R 4.3.1)
# gtable                 0.3.4     2023-08-21 [2] CRAN (R 4.3.1)
# here                 * 1.0.1     2020-12-13 [2] CRAN (R 4.3.1)
# IRanges              * 2.34.1    2023-06-22 [2] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [2] CRAN (R 4.3.1)
# labeling               0.4.3     2023-08-29 [2] CRAN (R 4.3.1)
# lattice                0.21-8    2023-04-05 [3] CRAN (R 4.3.1)
# lifecycle              1.0.3     2022-10-07 [2] CRAN (R 4.3.1)
# magrittr               2.0.3     2022-03-30 [2] CRAN (R 4.3.1)
# Matrix                 1.6-1.1   2023-09-18 [3] CRAN (R 4.3.1)
# MatrixGenerics       * 1.12.3    2023-07-30 [2] Bioconductor
# matrixStats          * 1.0.0     2023-06-02 [2] CRAN (R 4.3.1)
# munsell                0.5.0     2018-06-12 [2] CRAN (R 4.3.1)
# pillar                 1.9.0     2023-03-22 [2] CRAN (R 4.3.1)
# pkgconfig              2.0.3     2019-09-22 [2] CRAN (R 4.3.1)
# R6                     2.5.1     2021-08-19 [2] CRAN (R 4.3.1)
# Rcpp                   1.0.11    2023-07-06 [2] CRAN (R 4.3.1)
# RCurl                  1.98-1.12 2023-03-27 [2] CRAN (R 4.3.1)
# rlang                  1.1.1     2023-04-28 [2] CRAN (R 4.3.1)
# rprojroot              2.0.3     2022-04-02 [2] CRAN (R 4.3.1)
# rsvd                   1.0.5     2021-04-16 [2] CRAN (R 4.3.1)
# S4Arrays               1.0.6     2023-08-30 [2] Bioconductor
# S4Vectors            * 0.38.1    2023-05-02 [2] Bioconductor
# ScaledMatrix           1.8.1     2023-05-03 [2] Bioconductor
# scales                 1.2.1     2022-08-20 [2] CRAN (R 4.3.1)
# scater               * 1.28.0    2023-04-25 [2] Bioconductor
# scuttle              * 1.10.2    2023-08-03 [2] Bioconductor
# sessioninfo          * 1.2.2     2021-12-06 [2] CRAN (R 4.3.1)
# SingleCellExperiment * 1.22.0    2023-04-25 [2] Bioconductor
# sparseMatrixStats      1.12.2    2023-07-02 [2] Bioconductor
# SummarizedExperiment * 1.30.2    2023-06-06 [2] Bioconductor
# tibble                 3.2.1     2023-03-20 [2] CRAN (R 4.3.1)
# tidyselect             1.2.0     2022-10-10 [2] CRAN (R 4.3.1)
# utf8                   1.2.3     2023-01-31 [2] CRAN (R 4.3.1)
# vctrs                  0.6.3     2023-06-14 [2] CRAN (R 4.3.1)
# vipor                  0.4.5     2017-03-22 [2] CRAN (R 4.3.1)
# viridis                0.6.4     2023-07-22 [2] CRAN (R 4.3.1)
# viridisLite            0.4.2     2023-05-02 [2] CRAN (R 4.3.1)
# withr                  2.5.0     2022-03-03 [2] CRAN (R 4.3.1)
# XVector                0.40.0    2023-04-25 [2] Bioconductor
# zlibbioc               1.46.0    2023-04-25 [2] Bioconductor
#
# [1] /users/rphillip/R/4.3
# [2] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.3/R/lib64/R/library
#
# ────────────────────────────────────────────
#
