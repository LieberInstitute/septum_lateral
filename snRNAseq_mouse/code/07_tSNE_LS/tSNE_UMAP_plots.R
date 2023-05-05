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

keepnames <- levels(colData(sce.ls)$cellType.final)[-c(4:6, 24)]
LSnames <- levels(colData(sce.ls)$cellType.final)[10:18]

sce.ls.filter <- filterSCE(sce.ls, cellType.final %in% keepnames)
sce.ls.LS <- filterSCE(sce.ls, cellType.final %in% LSnames)

cell_colors.ls[names(cell_colors.ls)[-c(10:18)]] <- "#595E60"

