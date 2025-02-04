library("here")
library("sessioninfo")
library("biomaRt")
library("dplyr")
library("purrr")
library("stringr")
library("data.table")

source(
    here(
        "snRNAseq_mouse",
        "code",
        "04_clinical_set_enrichment",
        "reshape_modeling_results.R"
    )
)



####################### Load data for modeling_results ########################

## Load markers by Tran et al
load(
    here(
        "snRNAseq_mouse",
        "processed_data",
        "SCE",
        "markers-stats_LS-n4_findMarkers_33cellTypes.rda"
    ),
    verbose = TRUE
)
# Loading objects:
#   markers.ls.t.pw
#   markers.ls.t.1vAll
#   medianNon0.ls

lobstr::obj_size(markers.ls.t.pw)
# 1.22 GB
lobstr::obj_size(markers.ls.t.1vAll)
# 134.29 MB

class(markers.ls.t.1vAll[[1]])
# [1] "SimpleList"
# attr(,"package")
# [1] "S4Vectors"
names(markers.ls.t.1vAll$LS_In.R)
# [1] "0" "1"
class(markers.ls.t.pw)
# [1] "SimpleList"
# attr(,"package")
# [1] "S4Vectors"

###############################################################################



##### Edit markers.ls.t.1vAll to match markers.ls.t.1vAll.broad structure #####

# Do some reorganizing
markers.ls.t.1vAll <- lapply(markers.ls.t.1vAll, function(x) {
    # Basically take the 'stats.[1 or 0]' since is redundant with the 'summary'-level stats
    lapply(x, function(y) {
        y[, 4]
    })
})

# Re-name std.lfc column and the entries; add non-0-median info
for (i in names(markers.ls.t.1vAll)) {
    colnames(markers.ls.t.1vAll[[i]][["0"]])[1] <- "std.logFC"
    colnames(markers.ls.t.1vAll[[i]][["1"]])[1] <- "std.logFC"
    # Add non0median Boolean - might be informative for both sets of stats
    markers.ls.t.1vAll[[i]][["0"]] <- cbind(
        markers.ls.t.1vAll[[i]][["0"]],
        medianNon0.ls[[i]][match(
            rownames(markers.ls.t.1vAll[[i]][["0"]]),
            names(medianNon0.ls[[i]])
        )]
    )
    colnames(markers.ls.t.1vAll[[i]][["0"]])[4] <- "non0median"

    # "1" aka 'enriched'
    markers.ls.t.1vAll[[i]][["1"]] <- cbind(
        markers.ls.t.1vAll[[i]][["1"]],
        medianNon0.ls[[i]][match(
            rownames(markers.ls.t.1vAll[[i]][["1"]]),
            names(medianNon0.ls[[i]])
        )]
    )
    colnames(markers.ls.t.1vAll[[i]][["1"]])[4] <- "non0median"

    # Then re-name the entries to more interpretable, because we'll keeping both contrasts
    names(markers.ls.t.1vAll[[i]]) <- paste0(i, c("_depleted", "_enriched"))
}

###############################################################################



#################### Reshape markers.ls.t.1vAll (enriched) ####################

markers.ls.t.1vAll$LS_In.C$LS_In.C_enriched

markers.ls.t.1vAll_subset <- markers.ls.t.1vAll[grep("LS|Sept|MS|TNoS|TT\\.IG\\.SH|Chol", names(markers.ls.t.1vAll))]

modeling_result_1vsAll <- reshape_1vsAll(OnevsAll = markers.ls.t.1vAll_subset)

colSums(modeling_result_1vsAll[, grep("fdr_", colnames(modeling_result_1vsAll))] < 0.05)
#     fdr_Chol_Ex.D       fdr_LS_In.C       fdr_LS_In.D       fdr_LS_In.M
#               711              1235              1720              1114
#       fdr_LS_In.N       fdr_LS_In.O       fdr_LS_In.P       fdr_LS_In.Q
#               474               623              1631               594
#       fdr_LS_In.R       fdr_MS_In.J       fdr_MS_In.K     fdr_Sept_In.G
#              1275              1244              1592               788
#     fdr_Sept_In.I     fdr_TNoS_Ex.A fdr_TT.IG.SH_Ex.C fdr_TT.IG.SH_Ex.E
#              1574               848              1325               594
# fdr_TT.IG.SH_Ex.F
#              1729

colSums(modeling_result_1vsAll[, grep("fdr_", colnames(modeling_result_1vsAll))] < 0.1)
#     fdr_Chol_Ex.D       fdr_LS_In.C       fdr_LS_In.D       fdr_LS_In.M
#               780              1357              1761              1193
#       fdr_LS_In.N       fdr_LS_In.O       fdr_LS_In.P       fdr_LS_In.Q
#               534               694              1772               692
#       fdr_LS_In.R       fdr_MS_In.J       fdr_MS_In.K     fdr_Sept_In.G
#              1400              1416              1677               831
#     fdr_Sept_In.I     fdr_TNoS_Ex.A fdr_TT.IG.SH_Ex.C fdr_TT.IG.SH_Ex.E
#              1607               886              1431               656
# fdr_TT.IG.SH_Ex.F
#              1803

###############################################################################



########## Reshape markers.ls.t.pw for LS and Sept cell types (1vs1) ##########

markers.ls.t.pw

markers.ls.t.1vs1_subset <- markers.ls.t.pw[grep("LS|Sept|MS|TNoS|TT\\.IG\\.SH|Chol", names(markers.ls.t.pw))]

## Select genes with non0median == TRUE
non0med_genes <- lapply(markers.ls.t.1vs1_subset, function(x) {
    rownames(x[x$non0median == TRUE, ])
})

non0med_genes <- unique(unlist(non0med_genes))
non0med_genes <- non0med_genes[order(non0med_genes)]

## Select only LS and Sept stats for each LS and Sept
OnevsOne_modified <- lapply(markers.ls.t.1vs1_subset, function(celltype) {
    enriched <- celltype[non0med_genes, grep("LS|Sept|MS|TNoS|TT\\.IG\\.SH|median", names(celltype))]
    return(enriched)
})

modeling_result_1vs1 <- reshape_1vs1(OnevsOne_modified)

###############################################################################



######################## Create modeling_results object #######################

modeling_results <- list(
    "enrichment" = as.data.frame(modeling_result_1vsAll),
    "pairwise" = as.data.frame(modeling_result_1vs1)
)

###############################################################################



######################## Save modeling_results to rda #########################

save(modeling_results, file = here(
    "snRNAseq_mouse",
    "processed_data",
    "SCE",
    "gene_set_enrichment_objects.rda"
))

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
