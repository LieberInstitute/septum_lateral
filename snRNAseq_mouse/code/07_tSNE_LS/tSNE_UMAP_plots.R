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

all_names <- levels(colData(sce.ls)$cellType.final)[-grep("drop|mixed", levels(colData(sce.ls)$cellType.final))]
LS_names <- levels(colData(sce.ls)$cellType.final)[grep("LS", levels(colData(sce.ls)$cellType.final))]

sce.ls.filter <- filterSCE(sce.ls, cellType.final %in% all_names)
sce.ls.LS <- filterSCE(sce.ls, cellType.final %in% LS_names)

###############################################################################



############################### Plot aesthetics ###############################

my_theme <- theme_bw() +
    theme(
        text = element_text(size = 15),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()
    )

cell_colors.ls[-grep("LS", names(cell_colors.ls))] <- "#edede9"

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
    geom_point(size = 0.8, alpha = 0.3) +
    scale_color_manual(values = cell_colors.ls) +
    my_theme +
    coord_equal() +
    labs(x = "TSNE Dimension 1", y = "TSNE Dimension 2") +
    theme(legend.position = "None")

tSNE_LS_facet <- tSNE_LS_no_legend +
    facet_wrap(~cellType.final)

ggsave(
    tSNE_cellTypes_no_legend +
        tSNE_LS_facet +
        theme(axis.title.y = element_blank()),
    filename = here("snRNAseq_mouse","plots", "07_tSNE_LS", "tSNE_cellType_full_facet.pdf"),
    width = 13
)

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
    geom_point(size = 0.8, alpha = 0.3) +
    scale_color_manual(values = cell_colors.ls) +
    my_theme +
    coord_equal() +
    labs(x = "UMAP Dimension 1", y = "UMAP Dimension 2") +
    theme(legend.position = "None")

UMAP_LS_facet <- UMAP_LS_no_legend +
    facet_wrap(~cellType.final)

ggsave(
    UMAP_cellTypes_no_legend +
        UMAP_LS_facet +
        theme(axis.title.y = element_blank()),
    filename = here("snRNAseq_mouse","plots", "07_tSNE_LS", "UMAP_cellType_full_facet.pdf"),
    width = 13
)

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
