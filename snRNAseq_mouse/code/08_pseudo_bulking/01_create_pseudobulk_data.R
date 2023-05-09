library("here")
library("spatialLIBD")
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

