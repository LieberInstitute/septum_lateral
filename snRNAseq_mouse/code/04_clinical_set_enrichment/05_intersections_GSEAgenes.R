library("dplyr")
library("here")
library("sessioninfo")



################### Load lists of genes for each cell types ###################

load(
    here(
        "snRNAseq_mouse",
        "processed_data",
        "SCE",
        "enrichment_genes_glFDR01.rda"
    ),
    verbose = TRUE
)

# Loading objects:
#   broad_list_genes_neg
#   broad_list_genes_pos
#   list_genes_neg
#   list_genes_pos

head(broad_list_genes_neg)
#              ensembl    gene cell_type
# 1 ENSMUSG00000021448    Shc3     Astro
# 2 ENSMUSG00000035566  Pcdh17     Astro
# 3 ENSMUSG00000046447 Camk2n1     Astro
# 4 ENSMUSG00000020297    Nsg2      Chol
# 5 ENSMUSG00000021301   Hecw1      Chol
# 6 ENSMUSG00000022883   Robo1      Chol

unique(broad_list_genes_neg$cell_type)
#  [1] "Astro"      "Chol"       "ChP"        "Endo"       "Ependymal"
#  [6] "IoC"        "LS"         "Micro"      "MS"         "Mural"
# [11] "Neuroblast" "Oligo"      "OPC"        "Sept"       "Str"
# [16] "Thal"       "TNoS"       "TT.IG.SH"

unique(list_genes_neg$cell_type)
#  [1] "Chol_Ex.D"     "LS_In.C"       "LS_In.D"       "LS_In.M"
#  [5] "LS_In.N"       "LS_In.O"       "LS_In.P"       "LS_In.Q"
#  [9] "LS_In.R"       "MS_In.J"       "MS_In.K"       "Sept_In.G"
# [13] "Sept_In.I"     "TNoS_Ex.A"     "TT.IG.SH_Ex.C" "TT.IG.SH_Ex.E"
# [17] "TT.IG.SH_Ex.F"

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
