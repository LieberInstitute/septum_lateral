## Call libraries
library("here")
library("sessioninfo")
library("biomaRt")
library("dplyr")
library("purrr")
library("stringr")
library("data.table")


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
#     markers.ls.t.pw.broad
#     markers.ls.t.1vAll.broad
#     medianNon0.ls.broad

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

## Select only significant genes with media != 0 and FDR < 0.05
LS_enriched_names <- lapply(markers.ls.t.1vAll.broad, function(x) {
    rownames(x[[2]][x[[2]]$log.FDR < log(0.05) & x[[2]]$non0median == TRUE, ])
})

LS_enriched_names <- unique(unlist(LS_enriched_names))
LS_enriched_names <- LS_enriched_names[order(LS_enriched_names)]

## Convert S4Vector object to DataFrame and selecting needed columns
LS_enriched <- as.data.frame(markers.ls.t.1vAll.broad)
LS_enriched <- LS_enriched[LS_enriched_names, ]
LS_enriched <- LS_enriched %>%
    dplyr::select(contains("enriched"))


vecs <- seq(1, length(names(LS_enriched)), 4)
for (i in vecs){
    test_ls <- LS_enriched[,i:(i+3)]
    LS_enriched[test_ls[,4]==FALSE, i:(i+3)] <- NA
}

LS_enriched <- LS_enriched %>%
    dplyr::select(!contains("non0median"))

## Unlog pvalues and FDRs
colchang <- c(1:dim(LS_enriched)[2])[rep(c(FALSE, TRUE, TRUE))]

for (i in colchang) {
    LS_enriched[, i] <- exp(LS_enriched[, i])
}

## Create new column names
cells <- names(markers.ls.t.1vAll.broad)
new_names <- paste(rep(c("t_stat", "p_value", "fdr"), length(cells)), rep(cells, each = 3), sep = "_")

## Change column names
names(LS_enriched) <- new_names
LS_enriched$ensembl <- rownames(LS_enriched)

## Add names from mgi data base
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- LS_enriched$ensembl
gene_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "mgi_symbol"), values = genes, mart = mart)
LS_enriched <- merge(LS_enriched, gene_list, by.x = "ensembl", by.y = "ensembl_gene_id")
l_names <- length(colnames(LS_enriched))
LS_enriched <- LS_enriched[, c(2:(l_names - 1), 1, l_names)]
LS_enriched <- dplyr::rename(LS_enriched, gene = mgi_symbol)


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
LS_pairw <- as.data.frame(markers.ls.t.pw.broad$LS) %>%
    dplyr::filter(non0median == TRUE) %>%
    dplyr::select(matches("stats\\..+"))

## Unlog pvalues and FDRs
colchang <- c(1:dim(LS_pairw)[2])[rep(c(FALSE, TRUE, TRUE))]

for (i in colchang) {
    LS_pairw[, i] <- exp(LS_pairw[, i])
}

## Change column names
new_names <- paste(c("t_stat", "p_value", "fdr"), rep(stats_tiss, each = 3), sep = "")
names(LS_pairw) <- new_names
LS_pairw$ensembl <- rownames(LS_pairw)

## Add names from mgi data base
genes <- LS_pairw$ensembl
gene_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "mgi_symbol"), values = genes, mart = mart)
LS_pairw <- merge(LS_pairw, gene_list, by.x = "ensembl", by.y = "ensembl_gene_id")
l_names <- length(colnames(LS_pairw))
LS_pairw <- LS_pairw[, c(2:(l_names - 1), 1, l_names)]
LS_pairw <- dplyr::rename(LS_pairw, gene = mgi_symbol)


########################################
#### Create modeling_results object ####
########################################

modeling_results <- list(
    "enrichment" = as.data.frame(LS_enriched),
    "pairwise" = as.data.frame(LS_pairw)
)


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
