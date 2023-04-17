library("dplyr")
library("here")
library("sessioninfo")
library("purrr")

extract_sig_genes <- function(list_tiss, list_set, modeling_results) {
    dims <- dim(modeling_results)[2]
    lapply(list_tiss, function(list_tiss) {
        list_tiss <- cbind(list_tiss, modeling_results[, (dims - 1):dims])
        genes <- list_tiss %>%
            filter(.[, 3] < 0.05) %>%
            select(ensembl, gene)
        if (list_set == "positive") {
            inter_genes <- intersect(as.vector(genes$ensembl), gene_list_FDR01$positive)
        } else {
            inter_genes <- intersect(as.vector(genes$ensembl), gene_list_FDR01$negative)
        }
        genes <- genes %>% filter(ensembl %in% inter_genes)
        return(genes)
    })
}



############### Load objects with inputs from gene_set_enrichment #############

load(
    here(
        "snRNAseq_mouse",
        "processed_data",
        "SCE",
        "gene_set_enrichment_objects.rda"
    ),
    verbose = TRUE
)

###############################################################################



######### Extract gene names from each data set (broad and not-broad) #########

## Select genes by tissue type from broad data
broad_list_sep <- split.default(modeling_results_broad$enrichment[, 1:54], rep(1:18, each = 3))
ct_names <- unique(gsub(names(modeling_results_broad$enrichment)[-(55:56)], pattern = ".+\\_", replacement = ""))

broad_list_genes_neg <- extract_sig_genes(list_tiss = broad_list_sep, list_set = "negative", modeling_results = modeling_results_broad$enrichment)
names(broad_list_genes_neg) <- ct_names

broad_list_genes_pos <- extract_sig_genes(list_tiss = broad_list_sep, list_set = "positive", modeling_results = modeling_results_broad$enrichment)
names(broad_list_genes_pos) <- ct_names

## Select genes by tissue type from not-broad data
list_sep <- split.default(modeling_results$enrichment[, 1:48], rep(1:16, each = 3))
ct_clus_names <- gsub(names(modeling_results$enrichment)[-(49:50)], pattern = "._[a-z]+_", replacement = "")
ct_clus_names <- unique(gsub(ct_clus_names, pattern = "fdr\\_", replacement = ""))

list_genes_neg <- extract_sig_genes(list_tiss = list_sep, list_set = "negative", modeling_results = modeling_results$enrichment)
names(list_genes_neg) <- ct_clus_names

list_genes_pos <- extract_sig_genes(list_tiss = list_sep, list_set = "positive", modeling_results = modeling_results$enrichment)
names(list_genes_pos) <- ct_clus_names




