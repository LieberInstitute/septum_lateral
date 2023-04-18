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



################################ Intersections ################################

## LS vs all broad neuronal markers ##
# Actually the genes that I send you in tables are the genes unique for LS, I'm
# not sure if this is true for all cell types but it's true for LS. So the next
# chink fo code (L57 - L71) is the same as just L57 and L60-64
group_1 <- c("LS")
group_2 <- c("MS", "TNoS", "TT.IG.SH", "Chol")

# I'm selecting from the genes in the negative DE list (broad_list_genes_neg),
# you could change it to the positive (broad_list_genes_pos)
group_1 <- broad_list_genes_neg %>%
    filter(cell_type %in% group_1) %>%
    select(ensembl) %>%
    unique() %>%
    as.vector()
group_2 <- broad_list_genes_neg %>%
    filter(cell_type %in% group_2) %>%
    select(ensembl) %>%
    unique() %>%
    as.vector()

genes_un <- setdiff(group_1, group_2)

broad_list_genes_neg %>%
    filter(cell_type == "LS", ensembl %in% genes_un$ensembl)
#               ensembl    gene cell_type
# 1  ENSMUSG00000005034  Prkacb        LS
# 2  ENSMUSG00000008658  Rbfox1        LS
# 3  ENSMUSG00000018474    Chd3        LS
# 4  ENSMUSG00000020297    Nsg2        LS
# 5  ENSMUSG00000021301   Hecw1        LS
# 6  ENSMUSG00000021373    Cap2        LS
# 7  ENSMUSG00000021448    Shc3        LS
# 8  ENSMUSG00000021700   Rab3c        LS
# 9  ENSMUSG00000022285   Ywhaz        LS
# 10 ENSMUSG00000022883   Robo1        LS
# 11 ENSMUSG00000023236    Scg5        LS
# 12 ENSMUSG00000024268   Celf4        LS
# 13 ENSMUSG00000024897   Apba1        LS
# 14 ENSMUSG00000025104  Hdgfl3        LS
# 15 ENSMUSG00000026585  Kifap3        LS
# 16 ENSMUSG00000027546   Atp9a        LS
# 17 ENSMUSG00000028176   Lrrc7        LS
# 18 ENSMUSG00000028524   Sgip1        LS
# 19 ENSMUSG00000029053   Prkcz        LS
# 20 ENSMUSG00000029405   G3bp2        LS
# 21 ENSMUSG00000029516     Cit        LS
# 22 ENSMUSG00000030839  Sergef        LS
# 23 ENSMUSG00000033676  Gabrb3        LS
# 24 ENSMUSG00000034796   Cpne7        LS
# 25 ENSMUSG00000034958   Atcay        LS
# 26 ENSMUSG00000035566  Pcdh17        LS
# 27 ENSMUSG00000035653   Lrfn5        LS
# 28 ENSMUSG00000036667   Tcaf1        LS
# 29 ENSMUSG00000037386   Rims2        LS
# 30 ENSMUSG00000037492   Zmat4        LS
# 31 ENSMUSG00000039419 Cntnap2        LS
# 32 ENSMUSG00000040037   Negr1        LS
# 33 ENSMUSG00000040785    Ttc3        LS
# 34 ENSMUSG00000041852   Tcf20        LS
# 35 ENSMUSG00000046178   Nxph1        LS
# 36 ENSMUSG00000049583    Grm5        LS
# 37 ENSMUSG00000050587  Lrrc4c        LS
# 38 ENSMUSG00000052698    Tln2        LS
# 39 ENSMUSG00000056755    Grm7        LS
# 40 ENSMUSG00000061911   Myt1l        LS
# 41 ENSMUSG00000062296  Trank1        LS
# 42 ENSMUSG00000064293   Cntn4        LS
# 43 ENSMUSG00000097451    Rian        LS


###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
