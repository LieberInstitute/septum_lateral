### For iSEE shiny app =============
load(
    here(
        "snRNAseq_mouse",
        "processed_data",
        "SCE",
        "sce_updated_LS.rda"
    ),
    verbose = TRUE
)

# As-is size:
lobstr::obj_size(sce.ls) / 1024^3 # 7.17 GB

# Remove the clusters we're not reporting
sce.ls <- sce.ls[, -grep("drop.", sce.ls$cellType.final)]
sce.ls <- sce.ls[, -grep("Neuron.mixed", sce.ls$cellType.final)]
sce.ls$cellType.final <- droplevels(sce.ls$cellType.final)

# Corresponding hexadecimal colors
cell_colors.ls <-
    cell_colors.ls[-grep("drop", names(cell_colors.ls))]
cell_colors.ls <-
    cell_colors.ls[-grep("Neuron.mixed", names(cell_colors.ls))]

sce.hold <- sce.ls

# Remove 'counts' & 'binomial_pearson_residuals' (just keep 'logcounts')
assay(sce.ls, "counts") <- NULL
assay(sce.ls, "binomial_pearson_residuals") <- NULL
lobstr::obj_size(sce.ls) / 1024^3 # 0.82 GB

# 'GLMPCA_50' is just the top 50 PCs from the 100 computed in 'GLMPCA_approx'
reducedDim(sce.ls, "GLMPCA_approx") <- NULL
lobstr::obj_size(sce.ls) / 1024^3 # 0.80 GB

# Because rownames are currently in Ensembl ID, re-assign & 'uniquify' to gene symbols
rownames(sce.ls) <-
    uniquifyFeatureNames(rowData(sce.ls)$gene_id, rowData(sce.ls)$gene_name)

lobstr::obj_size(sce.ls) / 1024^3

## Some add'l metrics for rowData ===
#    As created/implemented in https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/master/shiny_apps/00_clean_functions.R
rowData(sce.ls)$propNucleiExprs <- apply(
    assay(sce.hold, "counts"),
    1,
    function(x) {
        round(mean(x != 0), 3)
    }
)

# The above, by cell type ===
cellType.idx <- splitit(sce.ls$cellType.final)
rowdat.sce <- rowData(sce.ls)
for (i in names(cellType.idx)) {
    message(Sys.time(), " computing propNucleiExprs for ", i)
    rowdat.sce[, paste0("propExprsIn.", i)] <- apply(
        assay(sce.hold, "counts")[, cellType.idx[[i]]],
        1,
        function(x) {
            round(mean(x != 0), 3)
        }
    )
}
rowData(sce.ls) <- rowdat.sce

## Clean up unnecessary info
keepCols <- setdiff(
    colnames(colData(sce.ls)),
    c(
        "high.mito",
        "sizeFactor",
        "clusters.glmpca",
        "cellType",
        "cellType.exp"
    )
)
colData(sce.ls) <- colData(sce.ls)[, keepCols]

lobstr::obj_size(sce.ls) / 1024^3 # 0.82 GB

# Re-name these
sce.ls.small <- sce.ls
cell_cols.clean <- cell_colors.ls

## Save this
save(
    sce.ls.small,
    cell_cols.clean,
    file = here(
        "snRNAseq_mouse",
        "processed_data",
        "SCE",
        "sce_for_iSEE_LS.rda"
    )
)
