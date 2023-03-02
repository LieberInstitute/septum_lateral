## Call libraries
library("here")
library("sessioninfo")
library("biomaRt")
library("dplyr")
library("purrr")
library("stringr")

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


#######################################################
#### Reshape markers.ls.t.1vAll.broad (enriched) ####
#######################################################

markers.ls.t.1vAll.broad$LS$LS_enriched

## Change column names
colnames(markers.ls.t.1vAll.broad$LS$LS_enriched)[1:3] <- c("t_stat", "p_value", "fdr")
markers.ls.t.1vAll.broad$LS$LS_enriched$ensembl <- rownames(markers.ls.t.1vAll.broad$LS$LS_enriched)

## Add names from mgi data base
mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- markers.ls.t.1vAll.broad$LS$LS_enriched$ensembl
gene_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "mgi_symbol"), values = genes, mart = mart)
LS_enriched <- merge(markers.ls.t.1vAll.broad$LS$LS_enriched, gene_list, by.x = "ensembl", by.y = "ensembl_gene_id")
LS_enriched <- LS_enriched[, c(2:4, 1, 6)]
colnames(LS_enriched)[5] <- "gene"


#####################################################
#### Reshape markers.ls.t.pw.broad$LS (pairwise) ####
#####################################################

markers.ls.t.pw.broad$LS
markers.ls.t.pw.broad$LS$Astro

## Creat list of name with comparations eg. _LS-Astro, _LS-
stats_tiss <- stringr::str_match(
    string = names(markers.ls.t.pw.broad$LS),
    pattern = "stats\\..+"
) %>%
    na.exclude() %>%
    as.vector() %>%
    str_remove("stats\\.")
stats_tiss <- paste("_", "LS-", stats_tiss, sep = "")

## Convert S4Vector object to DataFrame and selecting needed columns
LS_pairw<-as.data.frame(markers.ls.t.pw.broad$LS) %>% select(matches("stats\\..+"))

## Change column names
new_names<-paste(c("t_stat", "p_value", "fdr"), rep(stats_tiss,each = 3), sep = "")
names(LS_pairw)<-new_names
LS_pairw$ensembl<-rownames(LS_pairw)

## Add names from mgi data base
genes <- LS_pairw$ensembl
gene_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "mgi_symbol"), values = genes, mart = mart)
LS_pairw <- merge(LS_pairw, gene_list, by.x = "ensembl", by.y = "ensembl_gene_id")
l_names<-length(colnames(LS_pairw))
LS_pairw<-LS_pairw[, c(2:(l_names-1),1,l_names)]
LS_pairw<-rename(LS_pairw, gene = mgi_symbol)


########################################
#### Create modeling_results object ####
########################################

modeling_results <- list(
    "enrichment" = LS_enriched
    "pairwise" = LS_pairw
)



## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
