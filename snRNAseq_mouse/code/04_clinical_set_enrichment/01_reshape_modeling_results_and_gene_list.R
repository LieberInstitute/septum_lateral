library("here")
library("sessioninfo")
library("biomaRt")
library("dplyr")
library("purrr")
library("stringr")
library("data.table")

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

## Select genes with non0median == TRUE
non0med_genes <- lapply(markers.ls.t.1vAll.broad, function(x) {
    rownames(x[[2]][x[[2]]$non0median == TRUE, ])
})

non0med_genes <- unique(unlist(non0med_genes))
non0med_genes <- non0med_genes[order(non0med_genes)]

## Change pvalues. fdrs and t-stats for genes non0median == FALSE
markers.ls.t.1vAll.broad_modified <- lapply(markers.ls.t.1vAll.broad, function(celltype) {
    enriched <- celltype[[2]][non0med_genes, ]
    enriched$std.logFC[!enriched$non0median] <- 0
    enriched$log.p.value[!enriched$non0median] <- log(1)
    enriched$log.FDR[!enriched$non0median] <- log(1)
    return(enriched)
})

## Unlog p-values and FDRs
markers.ls.t.1vAll.broad_modified <- lapply(markers.ls.t.1vAll.broad_modified, function(enriched) {
    res <- enriched
    res$FDR <- exp(enriched$log.FDR)
    res$p.value <- exp(enriched$log.p.value)
    res <- res[, c("std.logFC", "p.value", "FDR")]
    return(res)
})

## Change column names
markers.ls.t.1vAll.broad_modified <- mapply(function(res, names_ct) {
    colnames(res) <- paste0(c("t_stat_", "p_value_", "fdr_"), names_ct)
    res$ensembl <- rownames(res)
    return(res)
}, markers.ls.t.1vAll.broad_modified, names(markers.ls.t.1vAll.broad_modified))

## Convert to data.frame
modeling_result_enrichment <-
    as.data.frame(Reduce(
        function(...) {
            merge(..., by = "ensembl")
        },
        markers.ls.t.1vAll.broad_modified
    ))

colSums(modeling_result_enrichment[, grep("fdr_", colnames(modeling_result_enrichment))] < 0.05)
#      fdr_Astro       fdr_Chol        fdr_ChP       fdr_Endo  fdr_Ependymal
#            488            700           1392            890           1439
#        fdr_IoC         fdr_LS      fdr_Micro         fdr_MS      fdr_Mural
#            607           1925            453           1757            253
# fdr_Neuroblast      fdr_Oligo        fdr_OPC       fdr_Sept        fdr_Str
#            602            388            633           1503           1592
#       fdr_Thal       fdr_TNoS   fdr_TT.IG.SH
#           1900            841           1790

colSums(modeling_result_enrichment[, grep("fdr_", colnames(modeling_result_enrichment))] < 0.1)
#      fdr_Astro       fdr_Chol        fdr_ChP       fdr_Endo  fdr_Ependymal
#            491            765           1496            908           1463
#        fdr_IoC         fdr_LS      fdr_Micro         fdr_MS      fdr_Mural
#            630           1945            466           1833            269
# fdr_Neuroblast      fdr_Oligo        fdr_OPC       fdr_Sept        fdr_Str
#            613            389            653           1547           1609
#       fdr_Thal       fdr_TNoS   fdr_TT.IG.SH
#           1978            872           1843

## Add names from mgi data base
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- modeling_result_enrichment$ensembl
gene_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "mgi_symbol"), values = genes, mart = mart)
modeling_result_enrichment <- merge(modeling_result_enrichment, gene_list, by.x = "ensembl", by.y = "ensembl_gene_id")
l_names <- length(colnames(modeling_result_enrichment))
modeling_result_enrichment <- modeling_result_enrichment[, c(2:(l_names - 1), 1, l_names)]
modeling_result_enrichment <- dplyr::rename(modeling_result_enrichment, gene = mgi_symbol)


#####################################################
#### Reshape markers.ls.t.pw.broad$LS (pairwise) ####
#####################################################

markers.ls.t.pw.broad$LS
markers.ls.t.pw.broad$LS$stats.Astro

## Creat list of name with comparations eg. _LS-Astro, _LS-
stats_tiss <- stringr::str_match(
    string = names(markers.ls.t.pw.broad$LS),
    pattern = "stats\\..+"
) %>%
    na.exclude() %>%
    as.vector() %>%
    str_remove("stats\\.")
stats_tiss <- paste("_", "LS-", stats_tiss, sep = "")

## Convert S4Vector object to DataFrame, select genes with non0median == TRUE
modeling_result_pairwise <- as.data.frame(markers.ls.t.pw.broad$LS) %>%
    dplyr::filter(non0median == TRUE) %>%
    dplyr::select(matches("stats\\..+"))

## Unlog pvalues and FDRs
colchang <- c(1:dim(modeling_result_pairwise)[2])[rep(c(FALSE, TRUE, TRUE))]

for (i in colchang) {
    modeling_result_pairwise[, i] <- exp(modeling_result_pairwise[, i])
}

## Change column names
new_names <- paste(c("t_stat", "p_value", "fdr"), rep(stats_tiss, each = 3), sep = "")
names(modeling_result_pairwise) <- new_names
modeling_result_pairwise$ensembl <- rownames(modeling_result_pairwise)

## Add names from mgi data base
genes <- modeling_result_pairwise$ensembl
gene_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "mgi_symbol"), values = genes, mart = mart)
modeling_result_pairwise <- merge(modeling_result_pairwise, gene_list, by.x = "ensembl", by.y = "ensembl_gene_id")
l_names <- length(colnames(modeling_result_pairwise))
modeling_result_pairwise <- modeling_result_pairwise[, c(2:(l_names - 1), 1, l_names)]
modeling_result_pairwise <- dplyr::rename(modeling_result_pairwise, gene = mgi_symbol)


########################################
#### Create modeling_results object ####
########################################

modeling_results <- list(
    "enrichment" = as.data.frame(modeling_result_enrichment),
    "pairwise" = as.data.frame(modeling_result_pairwise)
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


####################################################
#### Save modeling_results and gene_list to rda ####
####################################################

save(modeling_results, gene_list_FDR05, gene_list_FDR01, file = here(
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
