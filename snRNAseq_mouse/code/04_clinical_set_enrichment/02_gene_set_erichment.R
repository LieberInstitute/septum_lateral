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
# 32.48 MB
lobstr::obj_size(gene_list_FDR05)
# 360.60 kB
lobstr::obj_size(gene_list_FDR01)
# 114.46 kB

class(modeling_results)
# [1] "list"
names(modeling_results)
# [1] "enrichment" "pairwise"
dim(modeling_results$enrichment)
# [1] 32180    56
head(names(modeling_results$enrichment), 10)
# [1] "t_stat_Astro"  "p_value_Astro" "fdr_Astro"     "t_stat_Chol"
# [5] "p_value_Chol"  "fdr_Chol"      "t_stat_ChP"    "p_value_ChP"
# [9] "fdr_ChP"       "t_stat_Endo"
dim(modeling_results$pairwise)
# [1] 32180    53
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
enrichTab_FDR05 <- gene_set_enrichment(gene_list = gene_list_FDR05, modeling_results = modeling_results, model_type = "enrichment")
enrichTab_FDR05
#           OR          Pval       test NumSig SetSize       ID model_type fdr_cut
# 1  2.2599759 3.964764e-105      Astro   1514    3533      all enrichment     0.1
# 2  0.9869134  8.065248e-01      Astro    608    2282 positive enrichment     0.1
# 3  7.8611122 7.427975e-258      Astro    906    1251 negative enrichment     0.1
# 4  2.2279869 4.567026e-103       Chol   1559    3533      all enrichment     0.1
# 5  1.0089024  8.657723e-01       Chol    646    2282 positive enrichment     0.1
# 6  7.5596947 1.386367e-247       Chol    913    1251 negative enrichment     0.1
# 7  2.1268993  2.684898e-95        ChP   1704    3533      all enrichment     0.1
# 8  1.0423101  3.779766e-01        ChP    759    2282 positive enrichment     0.1
# 9  6.9812490 7.710855e-224        ChP    945    1251 negative enrichment     0.1
# 10 2.3007210 9.864797e-108       Endo   1467    3533      all enrichment     0.1
# 11 0.9677463  5.336649e-01       Endo    570    2282 positive enrichment     0.1
# 12 8.1630831 5.574128e-268       Endo    897    1251 negative enrichment     0.1
# 13 2.1339702  4.799598e-96  Ependymal   1696    3533      all enrichment     0.1
# 14 1.0422327  3.769876e-01  Ependymal    753    2282 positive enrichment     0.1
# 15 7.0071861 4.497475e-225  Ependymal    943    1251 negative enrichment     0.1
# 16 2.4200139 5.331306e-112        IoC   1300    3533      all enrichment     0.1
# 17 0.8903044  3.384069e-02        IoC    446    2282 positive enrichment     0.1
# 18 8.9365620 2.126570e-292        IoC    854    1251 negative enrichment     0.1
# 19 1.8184189  1.423792e-62         LS   1947    3533      all enrichment     0.1
# 20 1.0023921  9.648964e-01         LS    958    2282 positive enrichment     0.1
# 21 5.5627476 1.127876e-164         LS    989    1251 negative enrichment     0.1
# 22 2.3516968 1.482725e-111      Micro   1431    3533      all enrichment     0.1
# 23 0.9533194  3.629240e-01      Micro    539    2282 positive enrichment     0.1
# 24 8.5402218 1.566928e-279      Micro    892    1251 negative enrichment     0.1
# 25 1.8695259  2.541247e-68         MS   1925    3533      all enrichment     0.1
# 26 1.0186414  6.746207e-01         MS    939    2282 positive enrichment     0.1
# 27 5.7723514 6.555770e-173         MS    986    1251 negative enrichment     0.1
# 28 2.3110266 1.847097e-108      Mural   1459    3533      all enrichment     0.1
# 29 0.9617313  4.536927e-01      Mural    562    2282 positive enrichment     0.1
# 30 8.2796801 1.871784e-271      Mural    897    1251 negative enrichment     0.1
# 31 2.2190432 9.008214e-104 Neuroblast   1617    3533      all enrichment     0.1
# 32 1.0427470  3.785865e-01 Neuroblast    693    2282 positive enrichment     0.1
# 33 7.3517036 1.635762e-239 Neuroblast    924    1251 negative enrichment     0.1
# 34 2.5120112 4.006545e-113      Oligo   1178    3533      all enrichment     0.1
# 35 0.8549169  7.172847e-03      Oligo    373    2282 positive enrichment     0.1
# 36 9.0744614 1.100062e-294      Oligo    805    1251 negative enrichment     0.1
# 37 2.2408986 5.830909e-104        OPC   1543    3533      all enrichment     0.1
# 38 1.0032370  9.419695e-01        OPC    633    2282 positive enrichment     0.1
# 39 7.6533113 6.494587e-251        OPC    910    1251 negative enrichment     0.1
# 40 1.7168822  1.542163e-51       Sept   1989    3533      all enrichment     0.1
# 41 0.9645176  4.186754e-01       Sept    993    2282 positive enrichment     0.1
# 42 5.1954316 1.384587e-149       Sept    996    1251 negative enrichment     0.1
# 43 2.1759320 2.632373e-100        Str   1678    3533      all enrichment     0.1
# 44 1.0569148  2.326726e-01        Str    741    2282 positive enrichment     0.1
# 45 7.0985327 3.635928e-229        Str    937    1251 negative enrichment     0.1
# 46 1.9065337  2.335670e-72       Thal   1885    3533      all enrichment     0.1
# 47 1.0211059  6.405325e-01       Thal    906    2282 positive enrichment     0.1
# 48 5.9575403 3.752961e-181       Thal    979    1251 negative enrichment     0.1
# 49 2.2330526 8.822001e-105       TNoS   1598    3533      all enrichment     0.1
# 50 1.0352911  4.728048e-01       TNoS    677    2282 positive enrichment     0.1
# 51 7.4701133 9.526117e-244       TNoS    921    1251 negative enrichment     0.1
# 52 1.9626172  3.980205e-78   TT.IG.SH   1803    3533      all enrichment     0.1
# 53 1.0092863  8.391389e-01   TT.IG.SH    837    2282 positive enrichment     0.1
# 54 6.3411898 1.123542e-197   TT.IG.SH    966    1251 negative enrichment     0.1

enrichTab_FDR01 <- gene_set_enrichment(gene_list = gene_list_FDR01, modeling_results = modeling_results, model_type = "enrichment")
enrichTab_FDR01
#            OR         Pval       test NumSig SetSize       ID model_type fdr_cut
# 1   1.3534001 4.348035e-06      Astro    371    1125      all enrichment     0.1
# 2   1.1032190 1.659318e-01      Astro    300    1042 positive enrichment     0.1
# 3  16.2206756 6.225013e-29      Astro     71      83 negative enrichment     0.1
# 4   1.3820776 6.678158e-07       Chol    392    1125      all enrichment     0.1
# 5   1.1418408 5.408819e-02       Chol    321    1042 positive enrichment     0.1
# 6  15.2222278 1.323857e-27       Chol     71      83 negative enrichment     0.1
# 7   1.4406733 5.983318e-09        ChP    456    1125      all enrichment     0.1
# 8   1.2304695 1.740953e-03        ChP    385    1042 positive enrichment     0.1
# 9  12.4171902 1.516124e-23        ChP     71      83 negative enrichment     0.1
# 10  1.3012535 8.316965e-05       Endo    345    1125      all enrichment     0.1
# 11  1.0461185 5.393693e-01       Endo    275    1042 positive enrichment     0.1
# 12 15.8052452 3.333872e-29       Endo     70      83 negative enrichment     0.1
# 13  1.4473239 3.786815e-09  Ependymal    454    1125      all enrichment     0.1
# 14  1.2350054 1.355341e-03  Ependymal    383    1042 positive enrichment     0.1
# 15 12.5653289 8.973503e-24  Ependymal     71      83 negative enrichment     0.1
# 16  1.1498483 5.388232e-02        IoC    266    1125      all enrichment     0.1
# 17  0.8794310 1.147331e-01        IoC    201    1042 positive enrichment     0.1
# 18 13.4561645 1.959431e-28        IoC     65      83 negative enrichment     0.1
# 19  1.4953138 4.094660e-11         LS    580    1125      all enrichment     0.1
# 20  1.3300771 6.626039e-06         LS    508    1042 positive enrichment     0.1
# 21  9.1084941 2.915129e-17         LS     72      83 negative enrichment     0.1
# 22  1.2676300 5.307665e-04      Micro    325    1125      all enrichment     0.1
# 23  1.0076736 9.124802e-01      Micro    256    1042 positive enrichment     0.1
# 24 15.3675360 2.694636e-29      Micro     69      83 negative enrichment     0.1
# 25  1.4890048 6.396861e-11         MS    565    1125      all enrichment     0.1
# 26  1.3185259 1.281756e-05         MS    493    1042 positive enrichment     0.1
# 27  9.5694790 4.413039e-18         MS     72      83 negative enrichment     0.1
# 28  1.2904166 1.605560e-04      Mural    340    1125      all enrichment     0.1
# 29  1.0333222 6.379640e-01      Mural    270    1042 positive enrichment     0.1
# 30 16.0125682 1.763804e-29      Mural     70      83 negative enrichment     0.1
# 31  1.4395349 1.195636e-08 Neuroblast    420    1125      all enrichment     0.1
# 32  1.2081022 5.160946e-03 Neuroblast    349    1042 positive enrichment     0.1
# 33 14.2049518 3.424967e-26 Neuroblast     71      83 negative enrichment     0.1
# 34  1.1161084 1.477695e-01      Oligo    226    1125      all enrichment     0.1
# 35  0.8149441 1.844479e-02      Oligo    163    1042 positive enrichment     0.1
# 36 14.0613127 6.611338e-30      Oligo     63      83 negative enrichment     0.1
# 37  1.3926681 3.473426e-07        OPC    388    1125      all enrichment     0.1
# 38  1.1476981 4.489338e-02        OPC    317    1042 positive enrichment     0.1
# 39 15.5789023 4.376100e-28        OPC     71      83 negative enrichment     0.1
# 40  1.4698348 2.335569e-10       Sept    603    1125      all enrichment     0.1
# 41  1.3108947 1.836559e-05       Sept    530    1042 positive enrichment     0.1
# 42  9.2044773 1.522604e-16       Sept     73      83 negative enrichment     0.1
# 43  1.4582935 2.383836e-09        Str    446    1125      all enrichment     0.1
# 44  1.2397315 1.251018e-03        Str    375    1042 positive enrichment     0.1
# 45 13.0396028 1.715627e-24        Str     71      83 negative enrichment     0.1
# 46  1.4755281 2.049728e-10       Thal    545    1125      all enrichment     0.1
# 47  1.2983954 4.135567e-05       Thal    473    1042 positive enrichment     0.1
# 48 10.1866602 3.764042e-19       Thal     72      83 negative enrichment     0.1
# 49  1.4111640 8.276242e-08       TNoS    408    1125      all enrichment     0.1
# 50  1.1766988 1.662458e-02       TNoS    337    1042 positive enrichment     0.1
# 51 14.5918962 9.767752e-27       TNoS     71      83 negative enrichment     0.1
# 52  1.4532111 1.298745e-09   TT.IG.SH    508    1125      all enrichment     0.1
# 53  1.2625289 3.160090e-04   TT.IG.SH    436    1042 positive enrichment     0.1
# 54 11.4618552 3.133324e-21   TT.IG.SH     72      83 negative enrichment     0.1


###########################################################
#### gene_set_enrichment analysis with pairwise data ######
###########################################################

prwiseTab_FDR05 <- spatialLIBD::gene_set_enrichment(gene_list = gene_list_05, modeling_results = modeling_results, fdr_cut = 0.1, model_type = "pairwise", reverse = FALSE)
prwiseTab_FDR01 <- spatialLIBD::gene_set_enrichment(gene_list = gene_list_01, modeling_results = modeling_results, fdr_cut = 0.1, model_type = "pairwise", reverse = FALSE)

names(prwiseTab_FDR01) # same for both
# [1] "OR"         "Pval"       "test"       "NumSig"     "SetSize"
# [6] "ID"         "model_type" "fdr_cut"


###################################################
#### Save gene_set_enrichment results to rda ######
###################################################

save(enrichTab_FDR01, enrichTab_FDR05, prwiseTab_FDR01, prwiseTab_FDR05, file = here(
    "snRNAseq_mouse",
    "processed_data",
    "SCE",
    "gene_set_enrichment_result_tables.rda"
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
