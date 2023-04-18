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

create_df <- function(list_tiss, csv_name) {
    list_tiss <- Map(cbind, list_tiss, new_clumn = names(list_tiss))
    list_tiss <- bind_rows(list_tiss)
    colnames(list_tiss)[3] <- "cell_type"

    write.csv(
         x = list_tiss,
         file = csv_name,
         quote = FALSE,
         row.names = FALSE
     )

    return(list_tiss)
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
list_sep <- split.default(modeling_results$enrichment[, 1:51], rep(1:17, each = 3))
ct_clus_names <- gsub(names(modeling_results$enrichment)[-(52:53)], pattern = "._[a-z]+_", replacement = "")
ct_clus_names <- unique(gsub(ct_clus_names, pattern = "fdr\\_", replacement = ""))

list_genes_neg <- extract_sig_genes(list_tiss = list_sep, list_set = "negative", modeling_results = modeling_results$enrichment)
names(list_genes_neg) <- ct_clus_names

list_genes_pos <- extract_sig_genes(list_tiss = list_sep, list_set = "positive", modeling_results = modeling_results$enrichment)
names(list_genes_pos) <- ct_clus_names

###############################################################################



########################## Confirm results are right ##########################

unlist(lapply(broad_list_genes_neg, function(list_genes) {
    dim(list_genes)[1]
}))
# Astro       Chol        ChP       Endo  Ependymal        IoC         LS
#     3         17          4          3          4         24         43
# Micro         MS      Mural Neuroblast      Oligo        OPC       Sept
#     1         35          2         11          2         12         38
#   Str       Thal       TNoS   TT.IG.SH
#    47         49         25         40

unlist(lapply(broad_list_genes_pos, function(list_genes) {
    dim(list_genes)[1]
}))
# Astro       Chol        ChP       Endo  Ependymal        IoC         LS
#    21          5         42         70         39          3          4
# Micro         MS      Mural Neuroblast      Oligo        OPC       Sept
#    97          6         18         10         24         15          2
#   Str       Thal       TNoS   TT.IG.SH
#     3          8          3          8

unlist(lapply(list_genes_neg, function(list_genes) {
    dim(list_genes)[1]
}))
#     Chol_Ex.D       LS_In.C       LS_In.D       LS_In.M       LS_In.N
#            18            30            44            28            17
#       LS_In.O       LS_In.P       LS_In.Q       LS_In.R       MS_In.J
#            31            36            21            32            27
#       MS_In.K     Sept_In.G     Sept_In.I     TNoS_Ex.A TT.IG.SH_Ex.C
#            37            24            40            25            34
# TT.IG.SH_Ex.E TT.IG.SH_Ex.F
#            21            39

unlist(lapply(list_genes_pos, function(list_genes) {
    dim(list_genes)[1]
}))
#     Chol_Ex.D       LS_In.C       LS_In.D       LS_In.M       LS_In.N
#             5            10             2             3             1
#       LS_In.O       LS_In.P       LS_In.Q       LS_In.R       MS_In.J
#             1             6             1             3             1
#       MS_In.K     Sept_In.G     Sept_In.I     TNoS_Ex.A TT.IG.SH_Ex.C
#             5             4             2             3             7
# TT.IG.SH_Ex.E TT.IG.SH_Ex.F
#             3             3

###############################################################################



############################### Create data frame #############################

broad_list_genes_neg <- create_df(
    list_tiss = broad_list_genes_neg,
    csv_name = here(
        "snRNAseq_mouse",
        "processed_data",
        "tables",
        "GSEAgenes_glFDR01-negative_broad.csv"
    )
)
broad_list_genes_pos <- create_df(list_tiss = broad_list_genes_pos, csv_name = here(
        "snRNAseq_mouse",
        "processed_data",
        "tables",
        "GSEAgenes_glFDR01-positive_broad.csv"
    )
)
list_genes_neg <- create_df(list_tiss = list_genes_neg, csv_name = here(
        "snRNAseq_mouse",
        "processed_data",
        "tables",
        "GSEAgenes_glFDR01-negative_8clusts.csv"
    )
)
list_genes_pos <- create_df(list_tiss = list_genes_pos, csv_name = here(
        "snRNAseq_mouse",
        "processed_data",
        "tables",
        "GSEAgenes_glFDR01-positive_8clusts.csv"
    )
)

###############################################################################



############################### Save genes in rda #############################

save(broad_list_genes_neg,
    broad_list_genes_pos,
    list_genes_neg,
    list_genes_pos,
    file = here(
        "snRNAseq_mouse",
        "processed_data",
        "SCE",
        "enrichment_genes_glFDR01.rda"
    )
)

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
