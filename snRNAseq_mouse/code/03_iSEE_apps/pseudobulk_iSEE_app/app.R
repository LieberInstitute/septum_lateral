library("SingleCellExperiment")
library("iSEE")
library("shiny")
library("paletteer")
library("scuttle")

load("rse_gene_TrkB_KO_LS_n8_wm.Rdata", verbose = TRUE)

#topifnot(all(unique(rse_gene$SampleID)))

## Don't run this on app.R since we don't want to run this every single time
# lobstr::obj_size(rse_gene)
# 876.34 MB


## Make unique gene names
rownames(rse_gene) <-
    scuttle::uniquifyFeatureNames(rowData(rse_gene)$gencodeID, rowData(rse_gene)$Symbol)


source("initial.R", print.eval = TRUE)

## From https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/810b47364af4c8afe426bd2a6b559bd6a9f1cc98/shiny_apps/tran2021_AMY/app.R#L10-L14
## Related to https://github.com/iSEE/iSEE/issues/568
colData(rse_gene) <- cbind(
    colData(rse_gene)[, !colnames(colData(rse_gene)) %in% c("Condition", "SampleID")],
    colData(rse_gene)[, c("SampleID", "Condition")]
)

rse_gene$Condition <- as.factor(rse_gene$Condition)

#rse_gene <- registerAppOptions(rse_gene, color.maxlevels = length(cell_cols.clean))
iSEE(
    rse_gene,
    appTitle = "bulkRNA-seq_lateral_septum",
    initial = initial,
    colormap = ExperimentColorMap(colData = list(
        Sample = function(n) {
            cols <- paletteer::paletteer_d(
                palette = "RColorBrewer::Dark2",
                n = length(unique(rse_gene$Condition))
            )
            cols <- as.vector(cols)
            names(cols) <- levels(rse_gene$Condition)
            return(cols)
        }#,
        # SampleID = function(n) {
        #     return(cell_cols.clean)
        # }
    ))
)
