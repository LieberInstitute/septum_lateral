## Call libraries
library("here")
library("sessioninfo")
library("dplyr")
library("purrr")
library("spatialLIBD")


###########################################################
#### Load objects for gene_set_enrichment_plot_complex ####
###########################################################

## Load gene_set_enrihcment_result_tables.rda. The FDR in the file names refers to the gene set filtering, not the fdr_cut used in gene_set_enrichment which was 0.1.
load(
    here(
        "snRNAseq_mouse",
        "processed_data",
        "SCE",
        "gene_set_enrichment_result_tables.rda"
    ),
    verbose = TRUE
)
# Loading objects:
#  enrichTab_FDR01
#  enrichTab_FDR05
#  prwiseTab_FDR01
#  prwiseTab_FDR05

lobstr::obj_size(enrichTab_FDR01)
# 2.10 kB
lobstr::obj_size(enrichTab_FDR05)
# 2.10 kB
lobstr::obj_size(prwiseTab_FDR01)
# 5.70 kB
lobstr::obj_size(prwiseTab_FDR05)
# 5.70 kB

head(enrichTab_FDR01)
head(prwiseTab_FDR01)


#####################################
#### Reproducibility information ####
#####################################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
