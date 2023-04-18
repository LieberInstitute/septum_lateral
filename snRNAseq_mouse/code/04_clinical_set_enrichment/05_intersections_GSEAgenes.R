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

group_1 <- c("Micro")

# I'm selecting from the genes in the negative DE list (broad_list_genes_neg),
# you could change it to the positive (broad_list_genes_pos)
group_1 <- broad_list_genes_pos %>%
    filter(cell_type %in% group_1) %>%
    select(ensembl) %>%
    unique() %>%
    as.vector()
group_2 <- broad_list_genes_pos %>%
    filter(!cell_type %in% group_1) %>% ## when adding ! is the negative of the expression, so == all genes except the ones in group 1, which is Micro
    select(ensembl) %>%
    unique() %>%
    as.vector()

genes_un <- setdiff(group_1, group_2)

broad_list_genes_pos %>%
    filter(cell_type == "Micro", ensembl %in% genes_un$ensembl)
#               ensembl     gene cell_type
# 1  ENSMUSG00000000628      Hk2     Micro
# 2  ENSMUSG00000001995  Sipa1l2     Micro
# 3  ENSMUSG00000002699     Lcp2     Micro
# 4  ENSMUSG00000002985     Apoe     Micro
# 5  ENSMUSG00000003283      Hck     Micro
# 6  ENSMUSG00000007613   Tgfbr1     Micro
# 7  ENSMUSG00000008398     Elk3     Micro
# 8  ENSMUSG00000014361    Mertk     Micro
# 9  ENSMUSG00000015090    Ptgds     Micro
# 10 ENSMUSG00000015243    Abca1     Micro
# 11 ENSMUSG00000015745  Plekho1     Micro
# 12 ENSMUSG00000016087     Fli1     Micro
# 13 ENSMUSG00000017707  Serinc3     Micro
# 14 ENSMUSG00000018008    Cyth4     Micro
# 15 ENSMUSG00000018381     Abi3     Micro
# 16 ENSMUSG00000018654    Ikzf1     Micro
# 17 ENSMUSG00000019978  Epb41l2     Micro
# 18 ENSMUSG00000020101     Vsir     Micro
# 19 ENSMUSG00000020143    Dock2     Micro
# 20 ENSMUSG00000020422     Tns3     Micro
# 21 ENSMUSG00000020709    Adap2     Micro
# 22 ENSMUSG00000021190     Lgmn     Micro
# 23 ENSMUSG00000021423     Ly86     Micro
# 24 ENSMUSG00000021665     Hexb     Micro
# 25 ENSMUSG00000021939     Ctsb     Micro
# 26 ENSMUSG00000022108    Itm2b     Micro
# 27 ENSMUSG00000022148      Fyb     Micro
# 28 ENSMUSG00000022488  Nckap1l     Micro
# 29 ENSMUSG00000022817    Itgb5     Micro
# 30 ENSMUSG00000022952    Runx1     Micro
# 31 ENSMUSG00000024013     Fgd2     Micro
# 32 ENSMUSG00000024300    Myo1f     Micro
# 33 ENSMUSG00000024621    Csf1r     Micro
# 34 ENSMUSG00000024661     Fth1     Micro
# 35 ENSMUSG00000024965   Fermt3     Micro
# 36 ENSMUSG00000025017  Pik3ap1     Micro
# 37 ENSMUSG00000025283     Sat1     Micro
# 38 ENSMUSG00000026031    Cflar     Micro
# 39 ENSMUSG00000026288   Inpp5d     Micro
# 40 ENSMUSG00000026365      Cfh     Micro
# 41 ENSMUSG00000026395    Ptprc     Micro
# 42 ENSMUSG00000026786  Apbb1ip     Micro
# 43 ENSMUSG00000027422    Rrbp1     Micro
# 44 ENSMUSG00000027447     Cst3     Micro
# 45 ENSMUSG00000027695     Pld1     Micro
# 46 ENSMUSG00000027947    Il6ra     Micro
# 47 ENSMUSG00000028163    Nfkb1     Micro
# 48 ENSMUSG00000028382    Ptbp3     Micro
# 49 ENSMUSG00000028859    Csf3r     Micro
# 50 ENSMUSG00000028868    Wasf2     Micro
# 51 ENSMUSG00000029723   Spacdr     Micro
# 52 ENSMUSG00000029919    Hpgds     Micro
# 53 ENSMUSG00000029925   Tbxas1     Micro
# 54 ENSMUSG00000030064   Frmd4b     Micro
# 55 ENSMUSG00000030287    Itpr2     Micro
# 56 ENSMUSG00000030671    Pde3b     Micro
# 57 ENSMUSG00000030737  Slco2b1     Micro
# 58 ENSMUSG00000030786    Itgam     Micro
# 59 ENSMUSG00000030844    Rgs10     Micro
# 60 ENSMUSG00000032440   Tgfbr2     Micro
# 61 ENSMUSG00000033192   Lpcat2     Micro
# 62 ENSMUSG00000033335     Dnm2     Micro
# 63 ENSMUSG00000034330    Plcg2     Micro
# 64 ENSMUSG00000034663    Bmp2k     Micro
# 65 ENSMUSG00000034707      Gns     Micro
# 66 ENSMUSG00000034792    Gna15     Micro
# 67 ENSMUSG00000035158     Mitf     Micro
# 68 ENSMUSG00000035697 Arhgap45     Micro
# 69 ENSMUSG00000036353   P2ry12     Micro
# 70 ENSMUSG00000036833   Pnpla7     Micro
# 71 ENSMUSG00000036896     C1qc     Micro
# 72 ENSMUSG00000036908  Unc93b1     Micro
# 73 ENSMUSG00000038147     Cd84     Micro
# 74 ENSMUSG00000038642     Ctss     Micro
# 75 ENSMUSG00000039936   Pik3cd     Micro
# 76 ENSMUSG00000040940  Arhgef1     Micro
# 77 ENSMUSG00000041515     Irf8     Micro
# 78 ENSMUSG00000041797    Abca9     Micro
# 79 ENSMUSG00000042228      Lyn     Micro
# 80 ENSMUSG00000045962     Wnk1     Micro
# 81 ENSMUSG00000048120   Entpd1     Micro
# 82 ENSMUSG00000048163   Selplg     Micro
# 83 ENSMUSG00000048779    P2ry6     Micro
# 84 ENSMUSG00000052085    Dock8     Micro
# 85 ENSMUSG00000052336   Cx3cr1     Micro
# 86 ENSMUSG00000052384    Nrros     Micro
# 87 ENSMUSG00000055435      Maf     Micro
# 88 ENSMUSG00000055541    Lair1     Micro
# 89 ENSMUSG00000056917    Sipa1     Micro
# 90 ENSMUSG00000059182    Skap2     Micro
# 91 ENSMUSG00000062078      Qki     Micro
# 92 ENSMUSG00000063458    Lrmda     Micro
# 93 ENSMUSG00000075284    Wipf1     Micro
# 94 ENSMUSG00000079110    Capn3     Micro
# 95 ENSMUSG00000079227     Ccr5     Micro
# 96 ENSMUSG00000084796 Mir142hg     Micro
# 97 ENSMUSG00000098112     Bin2     Micro


## In.O vs all neuronal clusters (without O and D) and In.D vs all neuronal clusters (without O and D) ##
group_1 <- c("LS_In.O")
group_2 <- c("LS_In.D")

group_1 <- list_genes_neg %>%
    filter(cell_type %in% group_1) %>%
    select(ensembl) %>%
    unique() %>%
    as.vector()
group_2 <- list_genes_neg %>%
    filter(cell_type %in% group_2) %>%
    select(ensembl) %>%
    unique() %>%
    as.vector()
group_3 <- list_genes_neg %>%
    filter(!cell_type %in% group_1) %>%
    filter(!cell_type %in% group_2) %>%
    select(ensembl) %>%
    unique() %>%
    as.vector()

genes_un <- setdiff(setdiff(group_1$ensembl, group_2$ensembl),group_3)  ## Genes unique for LS-InO
list_genes_neg %>%
    filter(cell_type == "LS_In.O", ensembl %in% genes_un)
#              ensembl   gene cell_type
# 1 ENSMUSG00000009394   Syn2   LS_In.O
# 2 ENSMUSG00000022883  Robo1   LS_In.O
# 3 ENSMUSG00000027419  Pcsk2   LS_In.O
# 4 ENSMUSG00000029405  G3bp2   LS_In.O
# 5 ENSMUSG00000040536 Necab1   LS_In.O

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
