## Call libraries
library("here")
library("sessioninfo")
library("dplyr")
library("purrr")
library("spatialLIBD")


##############################################
#### Load objects for gene_set_enrichment ####
##############################################

## Load gene_set_enrichment_objects.rda
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
#  gene_list_FDR05
#  gene_list_FDR01

lobstr::obj_size(modeling_results)
# 1.56 MB
lobstr::obj_size(gene_list_FDR05)
# 360.60 kB
lobstr::obj_size(gene_list_FDR01)
# 114.46 kB

class(modeling_results)
# [1] "list"
names(modeling_results)
# [1] "enrichment" "pairwise"
dim(modeling_results$enrichment)
# [1] 149   56
head(names(modeling_results$enrichment), 10)
# [1] "t_stat_Astro"  "p_value_Astro" "fdr_Astro"     "t_stat_Chol"
# [5] "p_value_Chol"  "fdr_Chol"      "t_stat_ChP"    "p_value_ChP"
# [9] "fdr_ChP"       "t_stat_Endo"
dim(modeling_results$pairwise)
# [1] 2629   53
head(names(modeling_results$pairwise), 10)
# [1] "t_stat_LS-Astro"  "p_value_LS-Astro" "fdr_LS-Astro"     "t_stat_LS-Chol"
# [5] "p_value_LS-Chol"  "fdr_LS-Chol"      "t_stat_LS-ChP"    "p_value_LS-ChP"
# [9] "fdr_LS-ChP"       "t_stat_LS-Endo"

class(gene_list_FDR01)
# [1] "list"
names(gene_list_FDR01)
# [1] "all"      "positive" "negative"


###########################################################
#### gene_set_enrichment analysis with enrichment data ####
###########################################################

## Running gene_set_enrichment with two sets of DE genes, one with FDR < 0.05 the other with FDR < 0.01
enrichTab_glFDR05_ctFDR05 <- gene_set_enrichment(gene_list = gene_list_FDR05, modeling_results = modeling_results, model_type = "enrichment", fdr_cut = 0.05)
# Warning message:
# In gene_set_enrichment(gene_list = gene_list_FDR05, modeling_results = modeling_results,  :
#   Gene list with n < 25 may have insufficent power for enrichment analysis: all ,positive ,negative
enrichTab_glFDR05_ctFDR05
#            OR         Pval       test NumSig SetSize       ID model_type fdr_cut
# 1  1.22407931 0.3694498417      Astro     14      42      all enrichment    0.05
# 2  0.11737620 0.9990545105      Astro      1      18 positive enrichment    0.05
# 3  3.40178034 0.0066792454      Astro     13      24 negative enrichment    0.05
# 4  1.90759114 0.1174129272       Chol     10      42      all enrichment    0.05
# 5  0.00000000 1.0000000000       Chol      0      18 positive enrichment    0.05
# 6  5.15700377 0.0013063181       Chol     10      24 negative enrichment    0.05
# 7  1.46641134 0.2256063147        ChP     13      42      all enrichment    0.05
# 8  0.15066402 0.9966103492        ChP      1      18 positive enrichment    0.05
# 9  3.76609333 0.0042405208        ChP     12      24 negative enrichment    0.05
# 10 1.39335043 0.2635623928       Endo     13      42      all enrichment    0.05
# 11 0.14515899 0.9971600223       Endo      1      18 positive enrichment    0.05
# 12 3.59187107 0.0055385174       Endo     12      24 negative enrichment    0.05
# 13 1.02839719 0.5442373204  Ependymal     14      42      all enrichment    0.05
# 14 0.10270789 0.9995633441  Ependymal      1      18 positive enrichment    0.05
# 15 2.89828715 0.0161372168  Ependymal     13      24 negative enrichment    0.05
# 16 2.06469347 0.0901907483        IoC     10      42      all enrichment    0.05
# 17 0.00000000 1.0000000000        IoC      0      18 positive enrichment    0.05
# 18 5.57011539 0.0008701416        IoC     10      24 negative enrichment    0.05
# 19 0.98620273 0.5869975099         LS     14      42      all enrichment    0.05
# 20 0.09941449 0.9996419421         LS      1      18 positive enrichment    0.05
# 21 2.78931852 0.0197038223         LS     13      24 negative enrichment    0.05
# 22 1.63138784 0.1583996298      Micro     13      42      all enrichment    0.05
# 23 0.16261617 0.9951960335      Micro      1      18 positive enrichment    0.05
# 24 4.15757615 0.0023926456      Micro     12      24 negative enrichment    0.05
# 25 1.17059966 0.4126933028         MS     14      42      all enrichment    0.05
# 26 0.11346501 0.9992181529         MS      1      18 positive enrichment    0.05
# 27 3.26463396 0.0084384254         MS     13      24 negative enrichment    0.05
# 28 1.54072221 0.2014353986      Mural     12      42      all enrichment    0.05
# 29 0.00000000 1.0000000000      Mural      0      18 positive enrichment    0.05
# 30 4.61938078 0.0012757529      Mural     12      24 negative enrichment    0.05
# 31 1.55287890 0.1797139657 Neuroblast     14      42      all enrichment    0.05
# 32 0.13993582 0.9976247459 Neuroblast      1      18 positive enrichment    0.05
# 33 4.23836204 0.0017802272 Neuroblast     13      24 negative enrichment    0.05
# 34 1.84442168 0.1091358017      Oligo     12      42      all enrichment    0.05
# 35 0.00000000 1.0000000000      Oligo      0      18 positive enrichment    0.05
# 36 5.49178907 0.0004385825      Oligo     12      24 negative enrichment    0.05
# 37 1.84442168 0.1091358017        OPC     12      42      all enrichment    0.05
# 38 0.00000000 1.0000000000        OPC      0      18 positive enrichment    0.05
# 39 5.49178907 0.0004385825        OPC     12      24 negative enrichment    0.05
# 40 1.02839719 0.5442373204       Sept     14      42      all enrichment    0.05
# 41 0.10270789 0.9995633441       Sept      1      18 positive enrichment    0.05
# 42 2.89828715 0.0161372168       Sept     13      24 negative enrichment    0.05
# 43 1.17059966 0.4126933028        Str     14      42      all enrichment    0.05
# 44 0.11346501 0.9992181529        Str      1      18 positive enrichment    0.05
# 45 3.26463396 0.0084384254        Str     13      24 negative enrichment    0.05
# 46 1.17059966 0.4126933028       Thal     14      42      all enrichment    0.05
# 47 0.11346501 0.9992181529       Thal      1      18 positive enrichment    0.05
# 48 3.26463396 0.0084384254       Thal     13      24 negative enrichment    0.05
# 49 1.63227890 0.1671593338       TNoS     12      42      all enrichment    0.05
# 50 0.00000000 1.0000000000       TNoS      0      18 positive enrichment    0.05
# 51 4.88285312 0.0009097511       TNoS     12      24 negative enrichment    0.05
# 52 1.28109833 0.3273696551   TT.IG.SH     14      42      all enrichment    0.05
# 53 0.12146903 0.9988588920   TT.IG.SH      1      18 positive enrichment    0.05
# 54 3.54761387 0.0052365917   TT.IG.SH     13      24 negative enrichment    0.05

enrichTab_glFDR05_ctFDR1 <- gene_set_enrichment(gene_list = gene_list_FDR05, modeling_results = modeling_results, model_type = "enrichment", fdr_cut = 0.1)
# Warning message:
# In gene_set_enrichment(gene_list = gene_list_FDR05, modeling_results = modeling_results,  :
#   Gene list with n < 25 may have insufficent power for enrichment analysis: all ,positive ,negative
enrichTab_glFDR05_ctFDR1
#           OR         Pval       test NumSig SetSize       ID model_type fdr_cut
# 1  1.1705997 0.4126933028      Astro     14      42      all enrichment     0.1
# 2  0.1134650 0.9992181529      Astro      1      18 positive enrichment     0.1
# 3  3.2646340 0.0084384254      Astro     13      24 negative enrichment     0.1
# 4  2.2610736 0.0491844236       Chol     12      42      all enrichment     0.1
# 5  0.0000000 1.0000000000       Chol      0      18 positive enrichment     0.1
# 6  6.6859288 0.0001255451       Chol     12      24 negative enrichment     0.1
# 7  1.5528789 0.1797139657        ChP     14      42      all enrichment     0.1
# 8  0.1399358 0.9976247459        ChP      1      18 positive enrichment     0.1
# 9  4.2383620 0.0017802272        ChP     13      24 negative enrichment     0.1
# 10 1.5528789 0.1797139657       Endo     14      42      all enrichment     0.1
# 11 0.1399358 0.9976247459       Endo      1      18 positive enrichment     0.1
# 12 4.2383620 0.0017802272       Endo     13      24 negative enrichment     0.1
# 13 1.0283972 0.5442373204  Ependymal     14      42      all enrichment     0.1
# 14 0.1027079 0.9995633441  Ependymal      1      18 positive enrichment     0.1
# 15 2.8982871 0.0161372168  Ependymal     13      24 negative enrichment     0.1
# 16 2.1637840 0.0669201760        IoC     11      42      all enrichment     0.1
# 17 0.0000000 1.0000000000        IoC      0      18 positive enrichment     0.1
# 18 6.0960248 0.0003389912        IoC     11      24 negative enrichment     0.1
# 19 1.1165679 0.4550297780         LS     16      42      all enrichment     0.1
# 20 0.1915468 0.9978697019         LS      2      18 positive enrichment     0.1
# 21 2.9513830 0.0141513000         LS     14      24 negative enrichment     0.1
# 22 1.6313878 0.1583996298      Micro     13      42      all enrichment     0.1
# 23 0.1626162 0.9951960335      Micro      1      18 positive enrichment     0.1
# 24 4.1575761 0.0023926456      Micro     12      24 negative enrichment     0.1
# 25 1.0283972 0.5442373204         MS     14      42      all enrichment     0.1
# 26 0.1027079 0.9995633441         MS      1      18 positive enrichment     0.1
# 27 2.8982871 0.0161372168         MS     13      24 negative enrichment     0.1
# 28 1.6313878 0.1583996298      Mural     13      42      all enrichment     0.1
# 29 0.1626162 0.9951960335      Mural      1      18 positive enrichment     0.1
# 30 4.1575761 0.0023926456      Mural     12      24 negative enrichment     0.1
# 31 1.4073899 0.2486104507 Neuroblast     14      42      all enrichment     0.1
# 32 0.1302527 0.9983474914 Neuroblast      1      18 positive enrichment     0.1
# 33 3.8695748 0.0031218604 Neuroblast     13      24 negative enrichment     0.1
# 34 1.6322789 0.1671593338      Oligo     12      42      all enrichment     0.1
# 35 0.0000000 1.0000000000      Oligo      0      18 positive enrichment     0.1
# 36 4.8828531 0.0009097511      Oligo     12      24 negative enrichment     0.1
# 37 1.6322789 0.1671593338        OPC     12      42      all enrichment     0.1
# 38 0.0000000 1.0000000000        OPC      0      18 positive enrichment     0.1
# 39 4.8828531 0.0009097511        OPC     12      24 negative enrichment     0.1
# 40 1.1165679 0.4550297780       Sept     16      42      all enrichment     0.1
# 41 0.1915468 0.9978697019       Sept      2      18 positive enrichment     0.1
# 42 2.9513830 0.0141513000       Sept     14      24 negative enrichment     0.1
# 43 1.1705997 0.4126933028        Str     14      42      all enrichment     0.1
# 44 0.1134650 0.9992181529        Str      1      18 positive enrichment     0.1
# 45 3.2646340 0.0084384254        Str     13      24 negative enrichment     0.1
# 46 1.0283972 0.5442373204       Thal     14      42      all enrichment     0.1
# 47 0.1027079 0.9995633441       Thal      1      18 positive enrichment     0.1
# 48 2.8982871 0.0161372168       Thal     13      24 negative enrichment     0.1
# 49 1.7252536 0.1296543644       TNoS     13      42      all enrichment     0.1
# 50 0.1691186 0.9942952897       TNoS      1      18 positive enrichment     0.1
# 51 4.3785477 0.0017603419       TNoS     12      24 negative enrichment     0.1
# 52 1.1705997 0.4126933028   TT.IG.SH     14      42      all enrichment     0.1
# 53 0.1134650 0.9992181529   TT.IG.SH      1      18 positive enrichment     0.1
# 54 3.2646340 0.0084384254   TT.IG.SH     13      24 negative enrichment     0.1

enrichTab_glFDR01_ctFDR05 <- gene_set_enrichment(gene_list = gene_list_FDR01, modeling_results = modeling_results, model_type = "enrichment", fdr_cut = 0.05)
# Warning message:
# In gene_set_enrichment(gene_list = gene_list_FDR05, modeling_results = modeling_results,  :
#   Gene list with n < 25 may have insufficent power for enrichment analysis: all ,positive ,negative
enrichTab_glFDR01_ctFDR05
#           OR      Pval       test NumSig SetSize       ID model_type fdr_cut
# 1  0.3168288 0.9482450      Astro      1       8      all enrichment    0.05
# 2  0.0000000 1.0000000      Astro      0       7 positive enrichment    0.05
# 3        Inf 0.3020134      Astro      1       1 negative enrichment    0.05
# 4  0.6979478 0.7787950       Chol      1       8      all enrichment    0.05
# 5  0.0000000 1.0000000       Chol      0       7 positive enrichment    0.05
# 6        Inf 0.1677852       Chol      1       1 negative enrichment    0.05
# 7  0.4035145 0.9112833        ChP      1       8      all enrichment    0.05
# 8  0.0000000 1.0000000        ChP      0       7 positive enrichment    0.05
# 9        Inf 0.2550336        ChP      1       1 negative enrichment    0.05
# 10 0.3891796 0.9176773       Endo      1       8      all enrichment    0.05
# 11 0.0000000 1.0000000       Endo      0       7 positive enrichment    0.05
# 12       Inf 0.2617450       Endo      1       1 negative enrichment    0.05
# 13 0.2786209 0.9626091  Ependymal      1       8      all enrichment    0.05
# 14 0.0000000 1.0000000  Ependymal      0       7 positive enrichment    0.05
# 15       Inf 0.3288591  Ependymal      1       1 negative enrichment    0.05
# 16 0.7343107 0.7636698        IoC      1       8      all enrichment    0.05
# 17 0.0000000 1.0000000        IoC      0       7 positive enrichment    0.05
# 18       Inf 0.1610738        IoC      1       1 negative enrichment    0.05
# 19 0.2700412 0.9656004         LS      1       8      all enrichment    0.05
# 20 0.0000000 1.0000000         LS      0       7 positive enrichment    0.05
# 21       Inf 0.3355705         LS      1       1 negative enrichment    0.05
# 22 0.4346089 0.8971797      Micro      1       8      all enrichment    0.05
# 23 0.0000000 1.0000000      Micro      0       7 positive enrichment    0.05
# 24       Inf 0.2416107      Micro      1       1 negative enrichment    0.05
# 25 0.3066418 0.9522262         MS      1       8      all enrichment    0.05
# 26 0.0000000 1.0000000         MS      0       7 positive enrichment    0.05
# 27       Inf 0.3087248         MS      1       1 negative enrichment    0.05
# 28 0.4694621 0.8811520      Mural      1       8      all enrichment    0.05
# 29 0.0000000 1.0000000      Mural      0       7 positive enrichment    0.05
# 30       Inf 0.2281879      Mural      1       1 negative enrichment    0.05
# 31 0.3755773 0.9236644 Neuroblast      1       8      all enrichment    0.05
# 32 0.0000000 1.0000000 Neuroblast      0       7 positive enrichment    0.05
# 33       Inf 0.2684564 Neuroblast      1       1 negative enrichment    0.05
# 34 0.5304896 0.8530143      Oligo      1       8      all enrichment    0.05
# 35 0.0000000 1.0000000      Oligo      0       7 positive enrichment    0.05
# 36       Inf 0.2080537      Oligo      1       1 negative enrichment    0.05
# 37 0.5304896 0.8530143        OPC      1       8      all enrichment    0.05
# 38 0.0000000 1.0000000        OPC      0       7 positive enrichment    0.05
# 39       Inf 0.2080537        OPC      1       1 negative enrichment    0.05
# 40 0.2786209 0.9626091       Sept      1       8      all enrichment    0.05
# 41 0.0000000 1.0000000       Sept      0       7 positive enrichment    0.05
# 42       Inf 0.3288591       Sept      1       1 negative enrichment    0.05
# 43 0.3066418 0.9522262        Str      1       8      all enrichment    0.05
# 44 0.0000000 1.0000000        Str      0       7 positive enrichment    0.05
# 45       Inf 0.3087248        Str      1       1 negative enrichment    0.05
# 46 0.3066418 0.9522262       Thal      1       8      all enrichment    0.05
# 47 0.0000000 1.0000000       Thal      0       7 positive enrichment    0.05
# 48       Inf 0.3087248       Thal      1       1 negative enrichment    0.05
# 49 0.4885179 0.8723484       TNoS      1       8      all enrichment    0.05
# 50 0.0000000 1.0000000       TNoS      0       7 positive enrichment    0.05
# 51       Inf 0.2214765       TNoS      1       1 negative enrichment    0.05
# 52 0.3274883 0.9439766   TT.IG.SH      1       8      all enrichment    0.05
# 53 0.0000000 1.0000000   TT.IG.SH      0       7 positive enrichment    0.05
# 54       Inf 0.2953020   TT.IG.SH      1       1 negative enrichment    0.05

enrichTab_glFDR01_ctFDR1 <- gene_set_enrichment(gene_list = gene_list_FDR01, modeling_results = modeling_results, model_type = "enrichment", fdr_cut = 0.1)
# Warning message:
# In gene_set_enrichment(gene_list = gene_list_FDR05, modeling_results = modeling_results,  :
#   Gene list with n < 25 may have insufficent power for enrichment analysis: all ,positive ,negative
enrichTab_glFDR01_ctFDR1
#           OR      Pval       test NumSig SetSize       ID model_type fdr_cut
# 1  0.3066418 0.9522262      Astro      1       8      all enrichment     0.1
# 2  0.0000000 1.0000000      Astro      0       7 positive enrichment     0.1
# 3        Inf 0.3087248      Astro      1       1 negative enrichment     0.1
# 4  0.6049531 0.8192122       Chol      1       8      all enrichment     0.1
# 5  0.0000000 1.0000000       Chol      0       7 positive enrichment     0.1
# 6        Inf 0.1879195       Chol      1       1 negative enrichment     0.1
# 7  0.3755773 0.9236644        ChP      1       8      all enrichment     0.1
# 8  0.0000000 1.0000000        ChP      0       7 positive enrichment     0.1
# 9        Inf 0.2684564        ChP      1       1 negative enrichment     0.1
# 10 0.3755773 0.9236644       Endo      1       8      all enrichment     0.1
# 11 0.0000000 1.0000000       Endo      0       7 positive enrichment     0.1
# 12       Inf 0.2684564       Endo      1       1 negative enrichment     0.1
# 13 0.2786209 0.9626091  Ependymal      1       8      all enrichment     0.1
# 14 0.0000000 1.0000000  Ependymal      0       7 positive enrichment     0.1
# 15       Inf 0.3288591  Ependymal      1       1 negative enrichment     0.1
# 16 0.6644808 0.7930663        IoC      1       8      all enrichment     0.1
# 17 0.0000000 1.0000000        IoC      0       7 positive enrichment     0.1
# 18       Inf 0.1744966        IoC      1       1 negative enrichment     0.1
# 19 0.5724902 0.8556800         LS      2       8      all enrichment     0.1
# 20 0.2817702 0.9605895         LS      1       7 positive enrichment     0.1
# 21       Inf 0.3624161         LS      1       1 negative enrichment     0.1
# 22 0.4346089 0.8971797      Micro      1       8      all enrichment     0.1
# 23 0.0000000 1.0000000      Micro      0       7 positive enrichment     0.1
# 24       Inf 0.2416107      Micro      1       1 negative enrichment     0.1
# 25 0.2786209 0.9626091         MS      1       8      all enrichment     0.1
# 26 0.0000000 1.0000000         MS      0       7 positive enrichment     0.1
# 27       Inf 0.3288591         MS      1       1 negative enrichment     0.1
# 28 0.4346089 0.8971797      Mural      1       8      all enrichment     0.1
# 29 0.0000000 1.0000000      Mural      0       7 positive enrichment     0.1
# 30       Inf 0.2416107      Mural      1       1 negative enrichment     0.1
# 31 0.3503625 0.9345065 Neuroblast      1       8      all enrichment     0.1
# 32 0.0000000 1.0000000 Neuroblast      0       7 positive enrichment     0.1
# 33       Inf 0.2818792 Neuroblast      1       1 negative enrichment     0.1
# 34 0.4885179 0.8723484      Oligo      1       8      all enrichment     0.1
# 35 0.0000000 1.0000000      Oligo      0       7 positive enrichment     0.1
# 36       Inf 0.2214765      Oligo      1       1 negative enrichment     0.1
# 37 0.4885179 0.8723484        OPC      1       8      all enrichment     0.1
# 38 0.0000000 1.0000000        OPC      0       7 positive enrichment     0.1
# 39       Inf 0.2214765        OPC      1       1 negative enrichment     0.1
# 40 0.5724902 0.8556800       Sept      2       8      all enrichment     0.1
# 41 0.2817702 0.9605895       Sept      1       7 positive enrichment     0.1
# 42       Inf 0.3624161       Sept      1       1 negative enrichment     0.1
# 43 0.3066418 0.9522262        Str      1       8      all enrichment     0.1
# 44 0.0000000 1.0000000        Str      0       7 positive enrichment     0.1
# 45       Inf 0.3087248        Str      1       1 negative enrichment     0.1
# 46 0.2786209 0.9626091       Thal      1       8      all enrichment     0.1
# 47 0.0000000 1.0000000       Thal      0       7 positive enrichment     0.1
# 48       Inf 0.3288591       Thal      1       1 negative enrichment     0.1
# 49 0.4515242 0.8894196       TNoS      1       8      all enrichment     0.1
# 50 0.0000000 1.0000000       TNoS      0       7 positive enrichment     0.1
# 51       Inf 0.2348993       TNoS      1       1 negative enrichment     0.1
# 52 0.3066418 0.9522262   TT.IG.SH      1       8      all enrichment     0.1
# 53 0.0000000 1.0000000   TT.IG.SH      0       7 positive enrichment     0.1
# 54       Inf 0.3087248   TT.IG.SH      1       1 negative enrichment     0.1

###########################################################
#### gene_set_enrichment analysis with pairwise data ######
###########################################################

prwiseTab_glFDR05_ctFDR1 <- spatialLIBD::gene_set_enrichment(gene_list = gene_list_FDR05, modeling_results = modeling_results, fdr_cut = 0.1, model_type = "pairwise", reverse = FALSE)
prwiseTab_glFDR05_ctFDR1
#           OR         Pval          test NumSig SetSize       ID model_type fdr_cut
# 1  1.6408988 1.043270e-05      LS-Astro    549     655      all   pairwise     0.1
# 2  0.1703484 1.000000e+00      LS-Astro     33      84 positive   pairwise     0.1
# 3  3.2199728 3.240272e-18      LS-Astro    516     571 negative   pairwise     0.1
# 4  2.0054110 1.159890e-12       LS-Chol    241     655      all   pairwise     0.1
# 5  0.3318498 9.998850e-01       LS-Chol      9      84 positive   pairwise     0.1
# 6  2.4238205 2.883148e-18       LS-Chol    232     571 negative   pairwise     0.1
# 7  3.2252861 3.122345e-37        LS-ChP    391     655      all   pairwise     0.1
# 8  0.4256770 9.997809e-01        LS-ChP     18      84 positive   pairwise     0.1
# 9  4.1807278 3.797932e-49        LS-ChP    373     571 negative   pairwise     0.1
# 10 1.8517178 5.303671e-10       LS-Endo    500     655      all   pairwise     0.1
# 11 0.1666896 1.000000e+00       LS-Endo     22      84 positive   pairwise     0.1
# 12 3.1488872 1.720776e-24       LS-Endo    478     571 negative   pairwise     0.1
# 13 2.2160806 6.435243e-16  LS-Ependymal    500     655      all   pairwise     0.1
# 14 0.3409012 9.999996e-01  LS-Ependymal     32      84 positive   pairwise     0.1
# 15 3.2344669 3.800595e-27  LS-Ependymal    468     571 negative   pairwise     0.1
# 16 1.0408639 3.477428e-01        LS-IoC    383     655      all   pairwise     0.1
# 17 0.7608991 9.105472e-01        LS-IoC     43      84 positive   pairwise     0.1
# 18 1.0994627 1.740393e-01        LS-IoC    340     571 negative   pairwise     0.1
# 19 2.0711850 1.855156e-13      LS-Micro    501     655      all   pairwise     0.1
# 20 0.1929115 1.000000e+00      LS-Micro     23      84 positive   pairwise     0.1
# 21 3.4655868 7.091726e-29      LS-Micro    478     571 negative   pairwise     0.1
# 22 1.3992335 6.521933e-04         LS-MS    188     655      all   pairwise     0.1
# 23 0.7416690 8.861004e-01         LS-MS     16      84 positive   pairwise     0.1
# 24 1.5099353 7.178174e-05         LS-MS    172     571 negative   pairwise     0.1
# 25 1.9378542 1.124539e-10      LS-Mural    518     655      all   pairwise     0.1
# 26 0.1660181 1.000000e+00      LS-Mural     24      84 positive   pairwise     0.1
# 27 3.5174666 1.643899e-26      LS-Mural    494     571 negative   pairwise     0.1
# 28 2.1140241 1.306323e-14 LS-Neuroblast    490     655      all   pairwise     0.1
# 29 0.3211484 9.999999e-01 LS-Neuroblast     30      84 positive   pairwise     0.1
# 30 3.0640016 1.192549e-25 LS-Neuroblast    460     571 negative   pairwise     0.1
# 31 1.7360135 3.264765e-06      LS-Oligo    564     655      all   pairwise     0.1
# 32 0.2452151 1.000000e+00      LS-Oligo     43      84 positive   pairwise     0.1
# 33 3.1084634 5.619943e-16      LS-Oligo    521     571 negative   pairwise     0.1
# 34 1.7025483 2.846850e-08        LS-OPC    477     655      all   pairwise     0.1
# 35 0.2133901 1.000000e+00        LS-OPC     24      84 positive   pairwise     0.1
# 36 2.5783878 5.163492e-19        LS-OPC    453     571 negative   pairwise     0.1
# 37 1.2813397 4.355875e-03       LS-Sept    265     655      all   pairwise     0.1
# 38 0.5077274 9.976023e-01       LS-Sept     19      84 positive   pairwise     0.1
# 39 1.4586553 5.944022e-05       LS-Sept    246     571 negative   pairwise     0.1
# 40 0.5690231 1.000000e+00        LS-Str    315     655      all   pairwise     0.1
# 41 1.0446203 4.702155e-01        LS-Str     50      84 positive   pairwise     0.1
# 42 0.5341508 1.000000e+00        LS-Str    265     571 negative   pairwise     0.1
# 43 1.1027752 1.666060e-01       LS-Thal    216     655      all   pairwise     0.1
# 44 0.9242873 6.674294e-01       LS-Thal     25      84 positive   pairwise     0.1
# 45 1.1288800 1.245951e-01       LS-Thal    191     571 negative   pairwise     0.1
# 46 1.2572525 6.364512e-03       LS-TNoS    353     655      all   pairwise     0.1
# 47 0.5253335 9.983759e-01       LS-TNoS     29      84 positive   pairwise     0.1
# 48 1.4426951 6.826609e-05       LS-TNoS    324     571 negative   pairwise     0.1
# 49 0.6849095 9.999856e-01   LS-TT.IG.SH    251     655      all   pairwise     0.1
# 50 1.2168030 2.189068e-01   LS-TT.IG.SH     42      84 positive   pairwise     0.1
# 51 0.6339620 9.999991e-01   LS-TT.IG.SH    209     571 negative   pairwise     0.1

prwiseTab_glFDR05_ctFDR05 <- spatialLIBD::gene_set_enrichment(gene_list = gene_list_FDR05, modeling_results = modeling_results, fdr_cut = 0.05, model_type = "pairwise", reverse = FALSE)
prwiseTab_glFDR05_ctFDR05
#           OR         Pval          test NumSig SetSize       ID model_type fdr_cut
# 1  1.6720379 4.083137e-06      LS-Astro    547     655      all   pairwise    0.05
# 2  0.1773023 1.000000e+00      LS-Astro     33      84 positive   pairwise    0.05
# 3  3.2146725 1.282344e-18      LS-Astro    514     571 negative   pairwise    0.05
# 4  2.0339727 1.273491e-12       LS-Chol    226     655      all   pairwise    0.05
# 5  0.3702909 9.994959e-01       LS-Chol      9      84 positive   pairwise    0.05
# 6  2.4258728 1.232848e-17       LS-Chol    217     571 negative   pairwise    0.05
# 7  3.3463961 2.852401e-39        LS-ChP    384     655      all   pairwise    0.05
# 8  0.3918667 9.999086e-01        LS-ChP     16      84 positive   pairwise    0.05
# 9  4.3710952 5.291124e-52        LS-ChP    368     571 negative   pairwise    0.05
# 10 1.9286251 2.134752e-11       LS-Endo    494     655      all   pairwise    0.05
# 11 0.1593333 1.000000e+00       LS-Endo     20      84 positive   pairwise    0.05
# 12 3.2749406 7.405669e-27       LS-Endo    474     571 negative   pairwise    0.05
# 13 2.2301360 1.962634e-16  LS-Ependymal    494     655      all   pairwise    0.05
# 14 0.3605583 9.999985e-01  LS-Ependymal     32      84 positive   pairwise    0.05
# 15 3.1839315 3.384533e-27  LS-Ependymal    462     571 negative   pairwise    0.05
# 16 1.0679374 2.488543e-01        LS-IoC    368     655      all   pairwise    0.05
# 17 0.8140377 8.509758e-01        LS-IoC     42      84 positive   pairwise    0.05
# 18 1.1165255 1.338501e-01        LS-IoC    326     571 negative   pairwise    0.05
# 19 2.1609781 4.026090e-15      LS-Micro    497     655      all   pairwise    0.05
# 20 0.2067753 1.000000e+00      LS-Micro     23      84 positive   pairwise    0.05
# 21 3.5415652 1.280173e-30      LS-Micro    474     571 negative   pairwise    0.05
# 22 1.4101003 6.600753e-04         LS-MS    176     655      all   pairwise    0.05
# 23 0.7550295 8.679032e-01         LS-MS     15      84 positive   pairwise    0.05
# 24 1.5175599 8.509076e-05         LS-MS    161     571 negative   pairwise    0.05
# 25 1.9750288 1.224594e-11      LS-Mural    508     655      all   pairwise    0.05
# 26 0.1636968 1.000000e+00      LS-Mural     22      84 positive   pairwise    0.05
# 27 3.4885322 9.811314e-28      LS-Mural    486     571 negative   pairwise    0.05
# 28 2.0333014 1.212110e-13 LS-Neuroblast    478     655      all   pairwise    0.05
# 29 0.3428660 9.999994e-01 LS-Neuroblast     30      84 positive   pairwise    0.05
# 30 2.8409374 1.793606e-23 LS-Neuroblast    448     571 negative   pairwise    0.05
# 31 1.6993031 5.160260e-06      LS-Oligo    559     655      all   pairwise    0.05
# 32 0.2446743 1.000000e+00      LS-Oligo     42      84 positive   pairwise    0.05
# 33 2.9748566 1.505499e-15      LS-Oligo    517     571 negative   pairwise    0.05
# 34 1.7302799 7.224640e-09        LS-OPC    468     655      all   pairwise    0.05
# 35 0.2320416 1.000000e+00        LS-OPC     24      84 positive   pairwise    0.05
# 36 2.5441555 3.562047e-19        LS-OPC    444     571 negative   pairwise    0.05
# 37 1.3000827 3.240659e-03       LS-Sept    243     655      all   pairwise    0.05
# 38 0.5527004 9.924084e-01       LS-Sept     18      84 positive   pairwise    0.05
# 39 1.4604157 7.626810e-05       LS-Sept    225     571 negative   pairwise    0.05
# 40 0.6060305 1.000000e+00        LS-Str    310     655      all   pairwise    0.05
# 41 1.0742956 4.198723e-01        LS-Str     49      84 positive   pairwise    0.05
# 42 0.5691716 1.000000e+00        LS-Str    261     571 negative   pairwise    0.05
# 43 1.1555048 7.697660e-02       LS-Thal    206     655      all   pairwise    0.05
# 44 0.8039214 8.357082e-01       LS-Thal     21      84 positive   pairwise    0.05
# 45 1.2153616 3.182957e-02       LS-Thal    185     571 negative   pairwise    0.05
# 46 1.2984920 2.189601e-03       LS-TNoS    339     655      all   pairwise    0.05
# 47 0.4698648 9.996086e-01       LS-TNoS     25      84 positive   pairwise    0.05
# 48 1.5169831 6.822972e-06       LS-TNoS    314     571 negative   pairwise    0.05
# 49 0.7300100 9.997195e-01   LS-TT.IG.SH    237     655      all   pairwise    0.05
# 50 1.1534992 2.970611e-01   LS-TT.IG.SH     38      84 positive   pairwise    0.05
# 51 0.6870400 9.999522e-01   LS-TT.IG.SH    199     571 negative   pairwise    0.05

prwiseTab_glFDR01_ctFDR1 <- spatialLIBD::gene_set_enrichment(gene_list = gene_list_FDR01, modeling_results = modeling_results, fdr_cut = 0.1, model_type = "pairwise", reverse = FALSE)
# Warning message:
# In spatialLIBD::gene_set_enrichment(gene_list = gene_list_FDR01,  :
#   Gene list with n < 25 may have insufficent power for enrichment analysis: all ,positive ,negative
prwiseTab_glFDR01_ctFDR1
#            OR         Pval          test NumSig SetSize       ID model_type fdr_cut
# 1  1.09675227 4.481858e-01      LS-Astro     54      68      all   pairwise     0.1
# 2  0.15240225 9.999749e-01      LS-Astro      6      17 positive   pairwise     0.1
# 3  4.62236598 1.616703e-03      LS-Astro     48      51 negative   pairwise     0.1
# 4  2.44208843 3.337523e-04       LS-Chol     31      68      all   pairwise     0.1
# 5  0.00000000 1.000000e+00       LS-Chol      0      17 positive   pairwise     0.1
# 6  4.55640578 1.354870e-07       LS-Chol     31      51 negative   pairwise     0.1
# 7  2.82468583 2.513777e-05        LS-ChP     43      68      all   pairwise     0.1
# 8  0.09901993 9.997503e-01        LS-ChP      1      17 positive   pairwise     0.1
# 9  7.73127519 1.134842e-10        LS-ChP     42      51 negative   pairwise     0.1
# 10 1.29468118 2.084664e-01       LS-Endo     49      68      all   pairwise     0.1
# 11 0.10545830 9.999959e-01       LS-Endo      3      17 positive   pairwise     0.1
# 12 4.68363533 8.736880e-05       LS-Endo     46      51 negative   pairwise     0.1
# 13 1.89374219 1.495786e-02  LS-Ependymal     52      68      all   pairwise     0.1
# 14 0.31089689 9.956325e-01  LS-Ependymal      6      17 positive   pairwise     0.1
# 15 5.40139764 1.361399e-05  LS-Ependymal     46      51 negative   pairwise     0.1
# 16 0.72584176 9.235044e-01        LS-IoC     34      68      all   pairwise     0.1
# 17 0.22326843 9.991495e-01        LS-IoC      4      17 positive   pairwise     0.1
# 18 1.04645580 4.968944e-01        LS-IoC     30      51 negative   pairwise     0.1
# 19 1.13290877 3.684898e-01      LS-Micro     46      68      all   pairwise     0.1
# 20 0.07098842 9.999994e-01      LS-Micro      2      17 positive   pairwise     0.1
# 21 3.45719445 5.246602e-04      LS-Micro     44      51 negative   pairwise     0.1
# 22 1.64853001 3.991368e-02         LS-MS     23      68      all   pairwise     0.1
# 23 0.67992944 8.113465e-01         LS-MS      3      17 positive   pairwise     0.1
# 24 2.08529010 1.021112e-02         LS-MS     20      51 negative   pairwise     0.1
# 25 1.14403483 3.653126e-01      LS-Mural     49      68      all   pairwise     0.1
# 26 0.13425369 9.999849e-01      LS-Mural      4      17 positive   pairwise     0.1
# 27 3.37325668 1.327711e-03      LS-Mural     45      51 negative   pairwise     0.1
# 28 1.26200653 2.243614e-01 LS-Neuroblast     46      68      all   pairwise     0.1
# 29 0.12707781 9.999779e-01 LS-Neuroblast      3      17 positive   pairwise     0.1
# 30 3.28425819 5.254161e-04 LS-Neuroblast     43      51 negative   pairwise     0.1
# 31 0.73877902 8.871274e-01      LS-Oligo     51      68      all   pairwise     0.1
# 32 0.13313139 9.999920e-01      LS-Oligo      6      17 positive   pairwise     0.1
# 33 1.88095612 9.318584e-02      LS-Oligo     45      51 negative   pairwise     0.1
# 34 0.84850953 7.843558e-01        LS-OPC     41      68      all   pairwise     0.1
# 35 0.11875970 9.999878e-01        LS-OPC      3      17 positive   pairwise     0.1
# 36 1.65481507 7.474945e-02        LS-OPC     38      51 negative   pairwise     0.1
# 37 0.96472861 6.013315e-01       LS-Sept     24      68      all   pairwise     0.1
# 38 0.00000000 1.000000e+00       LS-Sept      0      17 positive   pairwise     0.1
# 39 1.58824829 6.862925e-02       LS-Sept     24      51 negative   pairwise     0.1
# 40 0.51851693 9.973631e-01        LS-Str     29      68      all   pairwise     0.1
# 41 0.49437999 9.545581e-01        LS-Str      7      17 positive   pairwise     0.1
# 42 0.53157885 9.911892e-01        LS-Str     22      51 negative   pairwise     0.1
# 43 0.78254970 8.455010e-01       LS-Thal     18      68      all   pairwise     0.1
# 44 0.67131887 8.313114e-01       LS-Thal      4      17 positive   pairwise     0.1
# 45 0.82448029 7.748364e-01       LS-Thal     14      51 negative   pairwise     0.1
# 46 1.89252726 7.901548e-03       LS-TNoS     44      68      all   pairwise     0.1
# 47 0.21604931 9.987329e-01       LS-TNoS      3      17 positive   pairwise     0.1
# 48 4.26633626 4.735321e-06       LS-TNoS     41      51 negative   pairwise     0.1
# 49 1.07691383 4.281123e-01   LS-TT.IG.SH     32      68      all   pairwise     0.1
# 50 0.65796536 8.587166e-01   LS-TT.IG.SH      6      17 positive   pairwise     0.1
# 51 1.26324800 2.457745e-01   LS-TT.IG.SH     26      51 negative   pairwise     0.1

prwiseTab_glFDR01_ctFDR05 <- spatialLIBD::gene_set_enrichment(gene_list = gene_list_FDR01, modeling_results = modeling_results, fdr_cut = 0.05, model_type = "pairwise", reverse = FALSE)
# Warning message:
# In spatialLIBD::gene_set_enrichment(gene_list = gene_list_FDR01,  :
#   Gene list with n < 25 may have insufficent power for enrichment analysis: all ,positive ,negative
prwiseTab_glFDR01_ctFDR05
#            OR         Pval          test NumSig SetSize       ID model_type fdr_cut
# 1  1.13933250 3.975817e-01      LS-Astro     54      68      all   pairwise    0.05
# 2  0.15825973 9.999657e-01      LS-Astro      6      17 positive   pairwise    0.05
# 3  4.79908272 1.166629e-03      LS-Astro     48      51 negative   pairwise    0.05
# 4  2.26723044 1.153739e-03       LS-Chol     28      68      all   pairwise    0.05
# 5  0.00000000 1.000000e+00       LS-Chol      0      17 positive   pairwise    0.05
# 6  3.97569804 1.675261e-06       LS-Chol     28      51 negative   pairwise    0.05
# 7  3.02533872 7.434545e-06        LS-ChP     43      68      all   pairwise    0.05
# 8  0.10584738 9.996169e-01        LS-ChP      1      17 positive   pairwise    0.05
# 9  8.27747150 2.391434e-11        LS-ChP     42      51 negative   pairwise    0.05
# 10 1.40680841 1.292798e-01       LS-Endo     49      68      all   pairwise    0.05
# 11 0.11446212 9.999913e-01       LS-Endo      3      17 positive   pairwise    0.05
# 12 5.08488523 3.075951e-05       LS-Endo     46      51 negative   pairwise    0.05
# 13 2.00078247 8.605933e-03  LS-Ependymal     52      68      all   pairwise    0.05
# 14 0.32821953 9.939744e-01  LS-Ependymal      6      17 positive   pairwise    0.05
# 15 5.70412655 6.291871e-06  LS-Ependymal     46      51 negative   pairwise    0.05
# 16 0.76732328 8.854943e-01        LS-IoC     33      68      all   pairwise    0.05
# 17 0.17398535 9.997218e-01        LS-IoC      3      17 positive   pairwise    0.05
# 18 1.17408288 3.398001e-01        LS-IoC     30      51 negative   pairwise    0.05
# 19 1.21236673 2.741420e-01      LS-Micro     46      68      all   pairwise    0.05
# 20 0.07589712 9.999989e-01      LS-Micro      2      17 positive   pairwise    0.05
# 21 3.69718541 2.403801e-04      LS-Micro     44      51 negative   pairwise    0.05
# 22 1.70075241 3.314542e-02         LS-MS     22      68      all   pairwise    0.05
# 23 0.46518818 9.188591e-01         LS-MS      2      17 positive   pairwise    0.05
# 24 2.30295000 4.241741e-03         LS-MS     20      51 negative   pairwise    0.05
# 25 1.18172669 3.161229e-01      LS-Mural     48      68      all   pairwise    0.05
# 26 0.10364280 9.999965e-01      LS-Mural      3      17 positive   pairwise    0.05
# 27 3.74635776 4.385622e-04      LS-Mural     45      51 negative   pairwise    0.05
# 28 1.34592566 1.560135e-01 LS-Neuroblast     46      68      all   pairwise    0.05
# 29 0.13538691 9.999620e-01 LS-Neuroblast      3      17 positive   pairwise    0.05
# 30 3.50047807 2.463002e-04 LS-Neuroblast     43      51 negative   pairwise    0.05
# 31 0.77365681 8.544023e-01      LS-Oligo     51      68      all   pairwise    0.05
# 32 0.13932599 9.999882e-01      LS-Oligo      6      17 positive   pairwise    0.05
# 33 1.96819034 7.423731e-02      LS-Oligo     45      51 negative   pairwise    0.05
# 34 0.86525656 7.615912e-01        LS-OPC     40      68      all   pairwise    0.05
# 35 0.12874977 9.999753e-01        LS-OPC      3      17 positive   pairwise    0.05
# 36 1.62049720 7.949879e-02        LS-OPC     37      51 negative   pairwise    0.05
# 37 0.91870802 6.693420e-01       LS-Sept     21      68      all   pairwise    0.05
# 38 0.00000000 1.000000e+00       LS-Sept      0      17 positive   pairwise    0.05
# 39 1.45327769 1.246151e-01       LS-Sept     21      51 negative   pairwise    0.05
# 40 0.56088459 9.932597e-01        LS-Str     29      68      all   pairwise    0.05
# 41 0.53387669 9.371493e-01        LS-Str      7      17 positive   pairwise    0.05
# 42 0.57465873 9.820413e-01        LS-Str     22      51 negative   pairwise    0.05
# 43 0.80496628 8.156985e-01       LS-Thal     17      68      all   pairwise    0.05
# 44 0.51843543 9.124611e-01       LS-Thal      3      17 positive   pairwise    0.05
# 45 0.91708341 6.595054e-01       LS-Thal     14      51 negative   pairwise    0.05
# 46 1.86062792 8.809218e-03       LS-TNoS     42      68      all   pairwise    0.05
# 47 0.14989251 9.996684e-01       LS-TNoS      2      17 positive   pairwise    0.05
# 48 4.22635273 3.287026e-06       LS-TNoS     40      51 negative   pairwise    0.05
# 49 1.16929317 3.034623e-01   LS-TT.IG.SH     31      68      all   pairwise    0.05
# 50 0.75693409 7.848540e-01   LS-TT.IG.SH      6      17 positive   pairwise    0.05
# 51 1.34419309 1.821273e-01   LS-TT.IG.SH     25      51 negative   pairwise    0.05

###################################################
#### Save gene_set_enrichment results to rda ######
###################################################

save(enrichTab_glFDR01_ctFDR05, enrichTab_glFDR01_ctFDR1, enrichTab_glFDR05_ctFDR05, enrichTab_glFDR05_ctFDR1, file = here(
    "snRNAseq_mouse",
    "processed_data",
    "SCE",
    "gene_set_enrichment_1vsAll_result_tables.rda"
))

save(prwiseTab_glFDR01_ctFDR05, prwiseTab_glFDR01_ctFDR1, prwiseTab_glFDR05_ctFDR05, prwiseTab_glFDR05_ctFDR1,
file = here(
    "snRNAseq_mouse",
    "processed_data",
    "SCE",
    "gene_set_enrichment_1vs1_result_tables.rda"
))


#####################################
#### Reproducibility information ####
#####################################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
