## Call libraries
library("here")
library("sessioninfo")
library("dplyr")
library("purrr")


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
#   gene_list

lobstr::obj_size(modeling_results)
# 423.75 MB
lobstr::obj_size(gene_list)
# 19.35 MB

class(modeling_results)
# [1] "list"
names(modeling_results)
# [1] "enrichment" "pairwise"
names(modeling_results$enrichment)
# [1] "t_stat_LS"  "p_value_LS" "fdr_LS"     "ensembl"    "gene"
class(gene_list)
# [1] "list"
names(gene_list)
# [1] "all"      "positive" "negative"





