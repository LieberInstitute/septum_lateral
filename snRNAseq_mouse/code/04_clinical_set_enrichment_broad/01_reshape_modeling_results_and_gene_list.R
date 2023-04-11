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
        "04_clinical_set_enrichment_broad",
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

###############################################################################



################# Reshape markers.ls.t.1vAll.broad (enriched) #################

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

###############################################################################



################# Reshape markers.ls.t.pw.broad$LS (pairwise) #################

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

###############################################################################



##################### Create modeling_results_broad object ####################

modeling_results_broad <- list(
    "enrichment" = as.data.frame(modeling_result_broad_1vsAll),
    "pairwise" = as.data.frame(modeling_result_broad_1vs1)
)

###############################################################################



###################### Load gene set data for gene lists ######################

sigGenes <- fread(file = "/dcl01/lieber/ajaffe/Keri/TrkBKO/tables/sigGenes_DE_FDR05.csv", header = TRUE, sep = ",", data.table = FALSE, stringsAsFactors = FALSE)

###############################################################################



############################# Make geneList object ############################

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



################# Save modeling_results and gene_list to rda ##################

save(modeling_results_broad, gene_list_FDR05, gene_list_FDR01, file = here(
    "snRNAseq_mouse",
    "processed_data",
    "SCE",
    "gene_set_enrichment_objects_broad.rda"
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
