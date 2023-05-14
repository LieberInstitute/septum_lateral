library("here")
library("SingleCellExperiment")
library("scater")
library("sessioninfo")

################### Load sce object and pseudobulk results ####################

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

load(
    here(
        "snRNAseq_mouse",
        "processed_data",
        "SCE",
        "sce_pseudobulking_LS_and_broad.rda"
    ),
    verbose = TRUE
)
# Loading objects:
#   sce_pseudo_LS
#   sce_pseudo_all

###############################################################################



############################## Correct LS colors ##############################

cell_colors_LS <- cell_colors.ls[grep("LS",names(cell_colors.ls))]

cell_colors_LS[1:9] <- c("#58B4E4","#169F74", "#0673B4", "#D56128", "#CC79A8", "#1878B6", "#AEC7E6","#F57E20", "#EFE642")

###############################################################################



############################ Plot LS and broad PCs ############################

## LS clusters
p_LS <- plotPCA(
    sce_pseudo_LS,
    colour_by = "cellType.final",
    ncomponents = 4,
    point_size = 2,
    label_format = c("%s %02i", " (%i%%)")
    #percentVar = metadata(sce_pseudo)$PCA_var_explained
    ) + scale_color_manual("cellType.final",values = cell_colors_LS)

pdf(file = here("snRNAseq_mouse/plots/sce_pseudo_LS_PCs.pdf"), width = 8, height = 8)
print(p_LS)
dev.off()

## Broad clusters
p_all <- plotPCA(
    sce_pseudo_all,
    colour_by = "cellType.broad",
    ncomponents = 12,
    point_size = 1,
    label_format = c("%s %02i", " (%i%%)"),
    #percentVar = metadata(sce_pseudo)$PCA_var_explained
    )

# This needs to be change
#p_all <- p_all + scale_color_manual("cellType.final", values = cell_colors.ls)

pdf(file = here("snRNAseq_mouse/plots/sce_pseudo_broad_PCs.pdf"), width = 14, height = 14)
print(p_all)
dev.off()

###############################################################################

