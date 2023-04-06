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



#################################### BROAD ####################################

########################################
#### Load data for modeling_results ####
########################################

## Load markers by Tran et al
load(
    here(
        "snRNAseq_mouse",
        "processed_data",
        "SCE",
        "markers-stats_LS-n4_findMarkers_33cellTypes.broad.rda"
    ),
    verbose = TRUE
)
# Loading objects:
#   markers.ls.t.pw.broad
#   markers.ls.t.1vAll.broad
#   medianNon0.ls.broad

lobstr::obj_size(markers.ls.t.pw.broad)
# 423.75 MB
lobstr::obj_size(markers.ls.t.1vAll.broad)
# 63.08 MB

class(markers.ls.t.1vAll.broad$LS$LS_enriched)
# [1] "DFrame"
# attr(,"package")
# [1] "S4Vectors"
names(markers.ls.t.1vAll.broad$LS$LS_enriched)
# [1] "std.logFC"   "log.p.value" "log.FDR"     "non0median"
class(markers.ls.t.pw.broad$LS)
# [1] "DFrame"
# attr(,"package")
# [1] "S4Vectors"
names(markers.ls.t.pw.broad$LS)
names(markers.ls.t.pw.broad$LS$stats.Astro)
# [1] "logFC"       "log.p.value" "log.FDR"


#####################################################
#### Reshape markers.ls.t.1vAll.broad (enriched) ####
#####################################################

markers.ls.t.1vAll.broad$LS$LS_enriched

modeling_result_broad_1vsAll <- reshape_1vsAll(OnevsAll = markers.ls.t.1vAll.broad)

colSums(modeling_result_broad_1vsAll[, grep("fdr_", colnames(modeling_result_broad_1vsAll))] < 0.05)
#      fdr_Astro       fdr_Chol        fdr_ChP       fdr_Endo  fdr_Ependymal
#            488            700           1392            890           1439
#        fdr_IoC         fdr_LS      fdr_Micro         fdr_MS      fdr_Mural
#            607           1925            453           1757            253
# fdr_Neuroblast      fdr_Oligo        fdr_OPC       fdr_Sept        fdr_Str
#            602            388            633           1503           1592
#       fdr_Thal       fdr_TNoS   fdr_TT.IG.SH
#           1900            841           1790

colSums(modeling_result_broad_1vsAll[, grep("fdr_", colnames(modeling_result_broad_1vsAll))] < 0.1)
#      fdr_Astro       fdr_Chol        fdr_ChP       fdr_Endo  fdr_Ependymal
#            491            765           1496            908           1463
#        fdr_IoC         fdr_LS      fdr_Micro         fdr_MS      fdr_Mural
#            630           1945            466           1833            269
# fdr_Neuroblast      fdr_Oligo        fdr_OPC       fdr_Sept        fdr_Str
#            613            389            653           1547           1609
#       fdr_Thal       fdr_TNoS   fdr_TT.IG.SH
#           1978            872           1843


#####################################################
#### Reshape markers.ls.t.pw.broad$LS (pairwise) ####
#####################################################

markers.ls.t.pw.broad$LS
markers.ls.t.pw.broad$LS$stats.Astro

## Convert S4Vector object to DataFrame, select genes with non0median == TRUE
OnevsOne_broad_modified <- as.data.frame(markers.ls.t.pw.broad$LS) %>%
    dplyr::filter(non0median == TRUE) %>%
    dplyr::select(matches("stats\\..+"))

## Unlog pvalues and FDRs
OnevsOne_broad_modified <- OnevsOne_broad_modified %>% mutate_at(vars(contains('FDR')), exp)
OnevsOne_broad_modified <- OnevsOne_broad_modified %>% mutate_at(vars(contains('value')), exp)

## Change column names
names(OnevsOne_broad_modified) <- gsub(names(OnevsOne_broad_modified), pattern = "stats\\.", replacement = "LS-")
names(OnevsOne_broad_modified) <- sapply(
    lapply(strsplit(names(OnevsOne_broad_modified), "\\.log"),
        rev),
    paste, collapse = "_"
    )
names(OnevsOne_broad_modified) <- gsub(names(OnevsOne_broad_modified), pattern = "FC", replacement = "t_stat")
names(OnevsOne_broad_modified) <- gsub(names(OnevsOne_broad_modified), pattern = "\\.p\\.value", replacement = "p_value")
names(OnevsOne_broad_modified) <- gsub(names(OnevsOne_broad_modified), pattern = "\\.FDR", replacement = "fdr")
OnevsOne_broad_modified$ensembl <- rownames(OnevsOne_broad_modified)
rownames(OnevsOne_broad_modified) <- NULL

modeling_result_broad_1vs1 <- add_gene_names(OnevsOne_broad_modified)


##############################################
#### Create modeling_results_broad object ####
##############################################

modeling_results_broad <- list(
    "enrichment" = as.data.frame(modeling_result_broad_1vsAll),
    "pairwise" = as.data.frame(modeling_result_broad_1vs1)
)

###############################################################################



################################# LS and Sept #################################

########################################
#### Load data for modeling_results ####
########################################

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


#############################################################################
#### Edit markers.ls.t.1vAll to match markers.ls.t.1vAll.broad structure ####
#############################################################################

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


###############################################
#### Reshape markers.ls.t.1vAll (enriched) ####
###############################################

markers.ls.t.1vAll$LS_In.C$LS_In.C_enriched

markers.ls.t.1vAll_subset <- markers.ls.t.1vAll[grep("LS|Sept", names(markers.ls.t.1vAll))]

modeling_result_1vsAll <- reshape_1vsAll(OnevsAll = markers.ls.t.1vAll_subset)

colSums(modeling_result_1vsAll[, grep("fdr_", colnames(modeling_result_1vsAll))] < 0.05)
# fdr_LS_In.C   fdr_LS_In.D   fdr_LS_In.M   fdr_LS_In.N   fdr_LS_In.O
#        1235          1720          1114           474           623
# fdr_LS_In.P   fdr_LS_In.Q   fdr_LS_In.R fdr_Sept_In.G fdr_Sept_In.I
#        1631           594          1275           788          1574

colSums(modeling_result_1vsAll[, grep("fdr_", colnames(modeling_result_1vsAll))] < 0.1)
# fdr_LS_In.C   fdr_LS_In.D   fdr_LS_In.M   fdr_LS_In.N   fdr_LS_In.O
#        1357          1761          1193           534           694
# fdr_LS_In.P   fdr_LS_In.Q   fdr_LS_In.R fdr_Sept_In.G fdr_Sept_In.I
#        1772           692          1400           831          1607


###################################################################
#### Reshape markers.ls.t.pw for LS and Sept cell types (1vs1) ####
###################################################################

markers.ls.t.pw

markers.ls.t.1vs1_subset <- markers.ls.t.pw[grep("LS|Sept", names(markers.ls.t.pw))]

## Select genes with non0median == TRUE
non0med_genes <- lapply(markers.ls.t.1vs1_subset, function(x) {
    rownames(x[x$non0median == TRUE, ])
})

non0med_genes <- unique(unlist(non0med_genes))
non0med_genes <- non0med_genes[order(non0med_genes)]

## Select only LS and Sept stats for each LS and Sept
OnevsOne_modified <- lapply(markers.ls.t.1vs1_subset, function(celltype) {
    enriched <- celltype[non0med_genes, grep("LS|Sept|median", names(celltype))]
    return(enriched)
})

## Change pvalues. fdrs and t-stats for genes non0median == FALSE
OnevsOne_modified <- lapply(OnevsOne_modified, function(enriched) {
    enriched <- as.data.frame(enriched)
    enriched[enriched$non0median == FALSE, grep("FC", names(enriched))] <- 0
    enriched[enriched$non0median == FALSE, grep("value", names(enriched))] <- log(1)
    enriched[enriched$non0median == FALSE, grep("FDR", names(enriched))] <- log(1)
    enriched <- enriched %>% select(-non0median)
    return(enriched)
})

## Un log pvalues and FDR
OnevsOne_modified <- lapply(OnevsOne_modified, function(enriched) {
    enriched <- enriched %>% mutate_at(vars(contains('FDR')), exp)
    enriched <- enriched %>% mutate_at(vars(contains('value')), exp)
    return(enriched)
})

## Change column names
OnevsOne_modified <- lapply(OnevsOne_modified, function(enriched) {
    names(enriched) <- gsub(names(enriched), pattern = "stats\\.", replacement = "__")
    return(enriched)
})

## Convert to data frame
OnevsOne_modified <- as.data.frame(OnevsOne_modified)

## Change column names
names(OnevsOne_modified) <- gsub(names(OnevsOne_modified), pattern = "\\.\\_\\_", replacement = "-")
names(OnevsOne_modified) <- sapply(
    lapply(strsplit(names(OnevsOne_modified), "\\.log"),
        rev),
    paste, collapse = "_"
    )
names(OnevsOne_modified) <- gsub(names(OnevsOne_modified), pattern = "FC", replacement = "t_stat")
names(OnevsOne_modified) <- gsub(names(OnevsOne_modified), pattern = "\\.p\\.value", replacement = "p_value")
names(OnevsOne_modified) <- gsub(names(OnevsOne_modified), pattern = "\\.FDR", replacement = "fdr")
OnevsOne_modified$ensembl <- rownames(OnevsOne_modified)
rownames(OnevsOne_modified) <- NULL

## Add names from mgi data base
modeling_result_1vs1 <- add_gene_names(OnevsOne_modified)


########################################
#### Create modeling_results object ####
########################################

modeling_results <- list(
    "enrichment" = as.data.frame(modeling_result_1vsAll),
    "pairwise" = as.data.frame(modeling_result_1vs1)
)

###############################################################################



################################# Gene lists ##################################

############################
#### Load gene set data ####
############################

sigGenes <- fread(file = "/dcl01/lieber/ajaffe/Keri/TrkBKO/tables/sigGenes_DE_FDR05.csv", header = TRUE, sep = ",", data.table = FALSE, stringsAsFactors = FALSE)


##############################
#### Make geneList object ####
##############################

gene_list_FDR05 <- list(
    all = sigGenes %>% dplyr::select(ensemblID) %>% unlist() %>% as.vector(),
    positive = sigGenes %>% filter(logFC > 0) %>% dplyr::select(ensemblID) %>% unlist() %>% as.vector(),
    negative = sigGenes %>% filter(logFC < 0) %>% dplyr::select(ensemblID) %>% unlist() %>% as.vector()
)

gene_list_FDR01 <- list(
    all = sigGenes %>% filter(adj.P.Val < 0.01) %>% dplyr::select(ensemblID) %>% unlist() %>% as.vector(),
    positive = sigGenes %>% filter(logFC > 0, adj.P.Val < 0.01) %>% dplyr::select(ensemblID) %>% unlist() %>% as.vector(),
    negative = sigGenes %>% filter(logFC < 0, adj.P.Val < 0.01) %>% dplyr::select(ensemblID) %>% unlist() %>% as.vector()
)

###############################################################################



####################################################
#### Save modeling_results and gene_list to rda ####
####################################################

save(modeling_results_broad, modeling_results, gene_list_FDR05, gene_list_FDR01, file = here(
    "snRNAseq_mouse",
    "processed_data",
    "SCE",
    "gene_set_enrichment_objects.rda"
))


#####################################
#### Reproducibility information ####
#####################################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
