## Call libraries
library("here")
library("sessioninfo")
library("dplyr")
library("purrr")
library("spatialLIBD")


##############################################
#### Load objects for gene_set_enrichment ####
##############################################

## Load markers by Tran et al
load(
    here(
        "snRNAseq_mouse",
        "processed_data",
        "SCE",
        "gene_set_enrichment_objects.rda"
    ),
    verbose = TRUE
)
# Loading objects:
#   modeling_results
#  gene_list_05
#  gene_list_01

lobstr::obj_size(modeling_results)
# 19.35 MB
lobstr::obj_size(gene_list_05)
# 360.60 kB
lobstr::obj_size(gene_list_01)
# 114.46 kB

class(modeling_results)
# [1] "list"
names(modeling_results)
# [1] "enrichment" "pairwise"
names(modeling_results$enrichment)
# [1] "t_stat_LS"  "p_value_LS" "fdr_LS"     "ensembl"    "gene"
class(gene_list_05)
# [1] "list"
names(gene_list_05)
# [1] "all"      "positive" "negative"


###########################################################
#### gene_set_enrichment analysis with enrichment data ####
###########################################################

## I only have LS and at some point the function converts my data frame to a vector,
## the only line I changes is L61: colnames(tstats)<-"LS"
gene_set_enrichment_modi <- function(gene_list, fdr_cut = 0.1, modeling_results, model_type = names(modeling_results)[1], reverse = FALSE) {
    model_results <- modeling_results[[model_type]]
    geneList_present <- lapply(gene_list, function(x) {
        x <- x[!is.na(x)]
        x[x %in% model_results$ensembl]
    })
    tstats <- as.data.frame(model_results[, grep("[f|t]_stat_", colnames(model_results))])
    colnames(tstats) <- "LS"
    fdrs <- as.data.frame(model_results[, grep("fdr_", colnames(model_results))])

    enrichTab <- do.call(rbind, lapply(
        seq(along.with = tstats),
        function(i) {
            layer <- tstats[, i] > 0 & fdrs[, i] < fdr_cut
            tabList <- lapply(geneList_present, function(g) {
                table(Set = factor(model_results$ensembl %in%
                    g, c(FALSE, TRUE)), Layer = factor(layer, c(
                    FALSE,
                    TRUE
                )))
            })
            enrichList <- lapply(tabList, fisher.test)
            o <- data.frame(
                OR = vapply(
                    enrichList, "[[", numeric(1),
                    "estimate"
                ), Pval = vapply(
                    enrichList, "[[",
                    numeric(1), "p.value"
                ), test = colnames(tstats)[i],
                NumSig = vapply(tabList, function(x) {
                    x[2, 2]
                }, integer(1)), SetSize = vapply(
                    geneList_present,
                    length, integer(1)
                ), stringsAsFactors = FALSE
            )
            o$ID <- gsub(".odds ratio", "", rownames(o))
            rownames(o) <- NULL
            return(o)
        }
    ))
    enrichTab$model_type <- model_type
    enrichTab$fdr_cut <- fdr_cut
    return(enrichTab)
}

## Running gene_set_enrichment with two sets of DE genes, one with FDR < 0.05 the other with FDR < 0.01

enrichTab_05 <- gene_set_enrichment_modi(gene_list = gene_list_05, modeling_results = modeling_results, model_type = "enrichment")
enrichTab_01 <- gene_set_enrichment_modi(gene_list = gene_list_01, modeling_results = modeling_results, model_type = "enrichment")

names(enrichTab_01) # same for both
# [1] "OR"         "Pval"       "test"       "NumSig"     "SetSize"
# [6] "ID"         "model_type" "fdr_cut"
