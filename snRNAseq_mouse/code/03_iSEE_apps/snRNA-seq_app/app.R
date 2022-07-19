
library("SingleCellExperiment")
library("iSEE")
library("shiny")
# library("RColorBrewer")

load("sce_for_iSEE_LS.rda", verbose = TRUE)

## Don't run this on app.R since we don't want to run this every single time
# lobstr::obj_size(sce.ls.small)
# 876.33 MB

source("initial.R", print.eval = TRUE)

## From https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/810b47364af4c8afe426bd2a6b559bd6a9f1cc98/shiny_apps/tran2021_AMY/app.R#L10-L14
## Related to https://github.com/iSEE/iSEE/issues/568
# colData(sce.ls.small) <- cbind(
#   colData(sce.ls.small)[, !colnames(colData(sce.ls.small)) %in% c("donor", "cellType.final")],
#   colData(sce.ls.small)[, c("cellType.final", "donor")]
# )


iSEE(
    sce.ls.small,
    appTitle = "mm_LS_2022",
    initial = initial,
    # colormap = ExperimentColorMap(colData = list(
    #     donor = function(n) {
    #         cols <- RColorBrewer::brewer.pal(8, "Dark2")
    #         names(cols) <- paste0("donor", seq_len(8))
    #         return(cols)
    #     },
    #     cellType.final = function(n) {
    #         cell_colors.ls[!grepl("drop", names(cell_colors))]
    #     }
    # ))
)
