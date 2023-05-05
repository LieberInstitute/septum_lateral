library("dplyr")
library("here")
library("SingleCellExperiment")
library("ggplot2")
library("patchwork")
library("scater")
library("CATALYST")
library("sessioninfo")

load(
    here(
        "snRNAseq_mouse",
        "processed_data",
        "SCE",
        "sce_updated_LS.rda"
    ),
    verbose = TRUE
)

# Loading objects:
#   sce.ls
#   annotationTab.ls
#   cell_colors.ls

sce.ls

