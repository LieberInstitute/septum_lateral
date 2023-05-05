library("dplyr")
library("here")
library("SingleCellExperiment")
library("ggplot2")
library("patchwork")
library("scater")
library("CATALYST")
library("sessioninfo")



############################### Load sce object ###############################

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

###############################################################################



############################## Filter sce object ##############################

## sce.ls.filter has all cell types except the ones that needed to be dropped
## sce.ls.LS has only LS clusters

keepnames <- levels(colData(sce.ls)$cellType.final)[-c(4:6, 24)]
LSnames <- levels(colData(sce.ls)$cellType.final)[10:18]

sce.ls.filter <- filterSCE(sce.ls, cellType.final %in% keepnames)
sce.ls.LS <- filterSCE(sce.ls, cellType.final %in% LSnames)

###############################################################################



############################### Plot aesthetics ###############################

my_theme <- theme_bw() +
    theme(
        text = element_text(size = 15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    )

cell_colors.ls[names(cell_colors.ls)[-c(10:18)]] <- "#595E60"

###############################################################################



################################# tSNE plots ##################################

tSNE_cellTypes_no_legend <- ggcells(sce.ls.filter, mapping = aes(x = TSNE.1, y = TSNE.2, colour = cellType.final)) +
    geom_point(size = 0.2, alpha = 0.3) +
    scale_color_manual(values = cell_colors.ls) +
    my_theme +
    coord_equal() +
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2") +
    theme(legend.position = "None")


tSNE_LS_no_legend <- ggcells(sce.ls.LS, mapping = aes(x = TSNE.1, y = TSNE.2, colour = cellType.final)) +
    geom_point(size = 0.2, alpha = 0.3) +
    scale_color_manual(values = cell_colors.ls) +
    my_theme +
    coord_equal() +
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2") +
    theme(legend.position = "None")

tSNE_LS_no_legend +
    facet_wrap(~cellType.final)

###############################################################################



################################# UMAP plots ##################################

UMAP_cellTypes_no_legend <- ggcells(sce.ls.filter, mapping = aes(x = UMAP.1, y = UMAP.2, colour = cellType.final)) +
    geom_point(size = 0.2, alpha = 0.3) +
    scale_color_manual(values = cell_colors.ls) +
    my_theme +
    coord_equal() +
    labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2") +
    theme(legend.position = "None")

UMAP_LS_no_legend <- ggcells(sce.ls.LS, mapping = aes(x = UMAP.1, y = UMAP.2, colour = cellType.final)) +
    geom_point(size = 0.2, alpha = 0.3) +
    scale_color_manual(values = cell_colors.ls) +
    my_theme +
    coord_equal() +
    labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2") +
    theme(legend.position = "None")

UMAP_LS_no_legend +
    facet_wrap(~cellType.final)

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
