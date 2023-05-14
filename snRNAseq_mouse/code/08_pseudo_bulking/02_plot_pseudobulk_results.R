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

