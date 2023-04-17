library("dplyr")
library("here")
library("sessioninfo")


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



