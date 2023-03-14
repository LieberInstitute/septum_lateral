## Call libraries
library("here")
library("sessioninfo")
library("dplyr")
library("purrr")
library("spatialLIBD")
library("ComplexHeatmap")
library("circlize")


###########################################################
#### Load objects for gene_set_enrichment_plot_complex ####
###########################################################

## Load gene_set_enrichment_objects.rda which includes the gene list (FDR < 0.05 and FDR < 0.01) and modeling results
load(
    here(
        "snRNAseq_mouse",
        "processed_data",
        "SCE",
        "gene_set_enrichment_objects.rda"
    ),
    verbose = TRUE
)

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
# 5.88 kB
lobstr::obj_size(enrichTab_FDR05)
# 5.88 kB
lobstr::obj_size(prwiseTab_FDR01)
# 5.70 kB
lobstr::obj_size(prwiseTab_FDR05)
# 5.70 kB

enrichTab_FDR01
#            OR         Pval       test NumSig SetSize       ID model_type fdr_cut
# 1   1.3534001 2.446558e-06      Astro    371    1125      all enrichment     0.1
# 2   1.1032190 8.501599e-02      Astro    300    1042 positive enrichment     0.1
# 3  16.2206756 6.225013e-29      Astro     71      83 negative enrichment     0.1
# 4   1.3820776 3.580013e-07       Chol    392    1125      all enrichment     0.1
# 5   1.1418408 2.888658e-02       Chol    321    1042 positive enrichment     0.1
# 6  15.2222278 1.323857e-27       Chol     71      83 negative enrichment     0.1
# 7   1.4406733 3.488106e-09        ChP    456    1125      all enrichment     0.1
# 8   1.2304695 9.125863e-04        ChP    385    1042 positive enrichment     0.1
# 9  12.4171902 1.516124e-23        ChP     71      83 negative enrichment     0.1
# 10  1.3012535 4.848377e-05       Endo    345    1125      all enrichment     0.1
# 11  1.0461185 2.748428e-01       Endo    275    1042 positive enrichment     0.1
# 12 15.8052452 3.333872e-29       Endo     70      83 negative enrichment     0.1
# 13  1.4473239 2.339093e-09  Ependymal    454    1125      all enrichment     0.1
# 14  1.2350054 7.658126e-04  Ependymal    383    1042 positive enrichment     0.1
# 15 12.5653289 8.973503e-24  Ependymal     71      83 negative enrichment     0.1
# 16  1.1498483 2.874094e-02        IoC    266    1125      all enrichment     0.1
# 17  0.8794310 9.519240e-01        IoC    201    1042 positive enrichment     0.1
# 18 13.4561645 1.959431e-28        IoC     65      83 negative enrichment     0.1
# 19  1.4953138 2.257400e-11         LS    580    1125      all enrichment     0.1
# 20  1.3300771 3.674878e-06         LS    508    1042 positive enrichment     0.1
# 21  9.1084941 2.769417e-17         LS     72      83 negative enrichment     0.1
# 22  1.2676300 2.812132e-04      Micro    325    1125      all enrichment     0.1
# 23  1.0076736 4.706433e-01      Micro    256    1042 positive enrichment     0.1
# 24 15.3675360 2.694636e-29      Micro     69      83 negative enrichment     0.1
# 25  1.4890048 3.797694e-11         MS    565    1125      all enrichment     0.1
# 26  1.3185259 7.280313e-06         MS    493    1042 positive enrichment     0.1
# 27  9.5694790 4.283736e-18         MS     72      83 negative enrichment     0.1
# 28  1.2904166 8.540136e-05      Mural    340    1125      all enrichment     0.1
# 29  1.0333222 3.356818e-01      Mural    270    1042 positive enrichment     0.1
# 30 16.0125682 1.763804e-29      Mural     70      83 negative enrichment     0.1
# 31  1.4395349 7.077805e-09 Neuroblast    420    1125      all enrichment     0.1
# 32  1.2081022 2.786138e-03 Neuroblast    349    1042 positive enrichment     0.1
# 33 14.2049518 3.424967e-26 Neuroblast     71      83 negative enrichment     0.1
# 34  1.1161084 8.037190e-02      Oligo    226    1125      all enrichment     0.1
# 35  0.8149441 9.928885e-01      Oligo    163    1042 positive enrichment     0.1
# 36 14.0613127 6.611338e-30      Oligo     63      83 negative enrichment     0.1
# 37  1.3926681 2.120236e-07        OPC    388    1125      all enrichment     0.1
# 38  1.1476981 2.478349e-02        OPC    317    1042 positive enrichment     0.1
# 39 15.5789023 4.376100e-28        OPC     71      83 negative enrichment     0.1
# 40  1.4698348 1.411437e-10       Sept    603    1125      all enrichment     0.1
# 41  1.3108947 1.007998e-05       Sept    530    1042 positive enrichment     0.1
# 42  9.2044773 1.179255e-16       Sept     73      83 negative enrichment     0.1
# 43  1.4582935 1.313811e-09        Str    446    1125      all enrichment     0.1
# 44  1.2397315 6.650515e-04        Str    375    1042 positive enrichment     0.1
# 45 13.0396028 1.715627e-24        Str     71      83 negative enrichment     0.1
# 46  1.4755281 1.131733e-10       Thal    545    1125      all enrichment     0.1
# 47  1.2983954 2.297553e-05       Thal    473    1042 positive enrichment     0.1
# 48 10.1866602 3.764042e-19       Thal     72      83 negative enrichment     0.1
# 49  1.4111640 4.988569e-08       TNoS    408    1125      all enrichment     0.1
# 50  1.1766988 9.073720e-03       TNoS    337    1042 positive enrichment     0.1
# 51 14.5918962 9.767752e-27       TNoS     71      83 negative enrichment     0.1
# 52  1.4532111 7.512207e-10   TT.IG.SH    508    1125      all enrichment     0.1
# 53  1.2625289 1.627375e-04   TT.IG.SH    436    1042 positive enrichment     0.1
# 54 11.4618552 3.133324e-21   TT.IG.SH     72      83 negative enrichment     0.1

enrichTab_FDR05
#           OR          Pval       test NumSig SetSize       ID model_type fdr_cut
# 1  2.2599759 2.764776e-105      Astro   1514    3533      all enrichment     0.1
# 2  0.9869134  6.139565e-01      Astro    608    2282 positive enrichment     0.1
# 3  7.8611122 7.427975e-258      Astro    906    1251 negative enrichment     0.1
# 4  2.2279869 2.772218e-103       Chol   1559    3533      all enrichment     0.1
# 5  1.0089024  4.354104e-01       Chol    646    2282 positive enrichment     0.1
# 6  7.5596947 1.386367e-247       Chol    913    1251 negative enrichment     0.1
# 7  2.1268993  1.430482e-95        ChP   1704    3533      all enrichment     0.1
# 8  1.0423101  1.904662e-01        ChP    759    2282 positive enrichment     0.1
# 9  6.9812490 7.710855e-224        ChP    945    1251 negative enrichment     0.1
# 10 2.3007210 6.990228e-108       Endo   1467    3533      all enrichment     0.1
# 11 0.9677463  7.506128e-01       Endo    570    2282 positive enrichment     0.1
# 12 8.1630831 5.574128e-268       Endo    897    1251 negative enrichment     0.1
# 13 2.1339702  3.220579e-96  Ependymal   1696    3533      all enrichment     0.1
# 14 1.0422327  1.914060e-01  Ependymal    753    2282 positive enrichment     0.1
# 15 7.0071861 4.497475e-225  Ependymal    943    1251 negative enrichment     0.1
# 16 2.4200139 3.855943e-112        IoC   1300    3533      all enrichment     0.1
# 17 0.8903044  9.850269e-01        IoC    446    2282 positive enrichment     0.1
# 18 8.9365620 2.126570e-292        IoC    854    1251 negative enrichment     0.1
# 19 1.8184189  7.871715e-63         LS   1947    3533      all enrichment     0.1
# 20 1.0023921  4.867174e-01         LS    958    2282 positive enrichment     0.1
# 21 5.5627476 8.433297e-165         LS    989    1251 negative enrichment     0.1
# 22 2.3516968 1.108816e-111      Micro   1431    3533      all enrichment     0.1
# 23 0.9533194  8.317802e-01      Micro    539    2282 positive enrichment     0.1
# 24 8.5402218 1.566928e-279      Micro    892    1251 negative enrichment     0.1
# 25 1.8695259  1.699108e-68         MS   1925    3533      all enrichment     0.1
# 26 1.0186414  3.456050e-01         MS    939    2282 positive enrichment     0.1
# 27 5.7723514 4.551572e-173         MS    986    1251 negative enrichment     0.1
# 28 2.3110266 1.213202e-108      Mural   1459    3533      all enrichment     0.1
# 29 0.9617313  7.878475e-01      Mural    562    2282 positive enrichment     0.1
# 30 8.2796801 1.871784e-271      Mural    897    1251 negative enrichment     0.1
# 31 2.2190432 6.430304e-104 Neuroblast   1617    3533      all enrichment     0.1
# 32 1.0427470  1.941295e-01 Neuroblast    693    2282 positive enrichment     0.1
# 33 7.3517036 1.635762e-239 Neuroblast    924    1251 negative enrichment     0.1
# 34 2.5120112 2.881395e-113      Oligo   1178    3533      all enrichment     0.1
# 35 0.8549169  9.969535e-01      Oligo    373    2282 positive enrichment     0.1
# 36 9.0744614 1.100062e-294      Oligo    805    1251 negative enrichment     0.1
# 37 2.2408986 3.467477e-104        OPC   1543    3533      all enrichment     0.1
# 38 1.0032370  4.818770e-01        OPC    633    2282 positive enrichment     0.1
# 39 7.6533113 6.494587e-251        OPC    910    1251 negative enrichment     0.1
# 40 1.7168822  9.593074e-52       Sept   1989    3533      all enrichment     0.1
# 41 0.9645176  8.013396e-01       Sept    993    2282 positive enrichment     0.1
# 42 5.1954316 1.156915e-149       Sept    996    1251 negative enrichment     0.1
# 43 2.1759320 1.454390e-100        Str   1678    3533      all enrichment     0.1
# 44 1.0569148  1.212916e-01        Str    741    2282 positive enrichment     0.1
# 45 7.0985327 3.635928e-229        Str    937    1251 negative enrichment     0.1
# 46 1.9065337  1.522027e-72       Thal   1885    3533      all enrichment     0.1
# 47 1.0211059  3.266636e-01       Thal    906    2282 positive enrichment     0.1
# 48 5.9575403 2.290225e-181       Thal    979    1251 negative enrichment     0.1
# 49 2.2330526 6.032571e-105       TNoS   1598    3533      all enrichment     0.1
# 50 1.0352911  2.398953e-01       TNoS    677    2282 positive enrichment     0.1
# 51 7.4701133 9.526117e-244       TNoS    921    1251 negative enrichment     0.1
# 52 1.9626172  2.388376e-78   TT.IG.SH   1803    3533      all enrichment     0.1
# 53 1.0092863  4.269354e-01   TT.IG.SH    837    2282 positive enrichment     0.1
# 54 6.3411898 8.155726e-198   TT.IG.SH    966    1251 negative enrichment     0.1

prwiseTab_FDR01
#           OR         Pval          test NumSig SetSize       ID model_type fdr_cut
# 1  0.8441017 9.976160e-01      LS-Astro    550    1125      all   pairwise     0.1
# 2  0.7452277 9.999987e-01      LS-Astro    478    1042 positive   pairwise     0.1
# 3  5.8306375 6.223408e-11      LS-Astro     72      83 negative   pairwise     0.1
# 4  1.4490954 1.237158e-09       LS-Chol    703    1125      all   pairwise     0.1
# 5  1.3873343 1.894707e-07       LS-Chol    641    1042 positive   pairwise     0.1
# 6  2.5412572 7.073768e-05       LS-Chol     62      83 negative   pairwise     0.1
# 7  0.9890193 5.838666e-01        LS-ChP    540    1125      all   pairwise     0.1
# 8  0.8630534 9.908986e-01        LS-ChP    466    1042 positive   pairwise     0.1
# 9  8.8498671 3.549271e-15        LS-ChP     74      83 negative   pairwise     0.1
# 10 0.5061589 1.000000e+00       LS-Endo    437    1125      all   pairwise     0.1
# 11 0.4319198 1.000000e+00       LS-Endo    367    1042 positive   pairwise     0.1
# 12 4.4072244 1.444519e-08       LS-Endo     70      83 negative   pairwise     0.1
# 13 0.8129481 9.996947e-01  LS-Ependymal    492    1125      all   pairwise     0.1
# 14 0.7034108 1.000000e+00  LS-Ependymal    420    1042 positive   pairwise     0.1
# 15 6.9222103 3.625275e-13  LS-Ependymal     72      83 negative   pairwise     0.1
# 16 1.1020732 6.368978e-02        LS-IoC    701    1125      all   pairwise     0.1
# 17 1.0685494 1.606085e-01        LS-IoC    642    1042 positive   pairwise     0.1
# 18 1.6350796 2.473962e-02        LS-IoC     59      83 negative   pairwise     0.1
# 19 0.2013188 1.000000e+00      LS-Micro    241    1125      all   pairwise     0.1
# 20 0.1435482 1.000000e+00      LS-Micro    170    1042 positive   pairwise     0.1
# 21 4.6136307 1.064952e-08      LS-Micro     71      83 negative   pairwise     0.1
# 22 1.5213887 3.196622e-12         LS-MS    584    1125      all   pairwise     0.1
# 23 1.4843475 2.363982e-10         LS-MS    535    1042 positive   pairwise     0.1
# 24 2.0046850 1.187222e-03         LS-MS     49      83 negative   pairwise     0.1
# 25 0.5164756 1.000000e+00      LS-Mural    443    1125      all   pairwise     0.1
# 26 0.4457561 1.000000e+00      LS-Mural    375    1042 positive   pairwise     0.1
# 27 3.6987841 2.504991e-07      LS-Mural     68      83 negative   pairwise     0.1
# 28 0.6193723 1.000000e+00 LS-Neuroblast    434    1125      all   pairwise     0.1
# 29 0.5367588 1.000000e+00 LS-Neuroblast    368    1042 positive   pairwise     0.1
# 30 3.9034934 2.537315e-08 LS-Neuroblast     66      83 negative   pairwise     0.1
# 31 0.4454148 1.000000e+00      LS-Oligo    454    1125      all   pairwise     0.1
# 32 0.3799063 1.000000e+00      LS-Oligo    382    1042 positive   pairwise     0.1
# 33 4.4486383 6.161884e-08      LS-Oligo     72      83 negative   pairwise     0.1
# 34 0.6891235 1.000000e+00        LS-OPC    499    1125      all   pairwise     0.1
# 35 0.6154754 1.000000e+00        LS-OPC    434    1042 positive   pairwise     0.1
# 36 3.1707107 2.025812e-06        LS-OPC     65      83 negative   pairwise     0.1
# 37 0.9646414 7.286979e-01       LS-Sept    441    1125      all   pairwise     0.1
# 38 0.9118632 9.273803e-01       LS-Sept    395    1042 positive   pairwise     0.1
# 39 1.8654508 3.230239e-03       LS-Sept     46      83 negative   pairwise     0.1
# 40 1.2230371 5.525295e-04        LS-Str    647    1125      all   pairwise     0.1
# 41 1.2769614 6.821932e-05        LS-Str    610    1042 positive   pairwise     0.1
# 42 0.7211399 9.446000e-01        LS-Str     37      83 negative   pairwise     0.1
# 43 1.0034017 4.895806e-01       LS-Thal    498    1125      all   pairwise     0.1
# 44 1.0475984 2.408912e-01       LS-Thal    472    1042 positive   pairwise     0.1
# 45 0.5754190 9.939541e-01       LS-Thal     26      83 negative   pairwise     0.1
# 46 1.1854529 3.069895e-03       LS-TNoS    661    1125      all   pairwise     0.1
# 47 1.1181899 4.215921e-02       LS-TNoS    598    1042 positive   pairwise     0.1
# 48 2.6115074 5.303828e-05       LS-TNoS     63      83 negative   pairwise     0.1
# 49 1.2118497 8.738978e-04   LS-TT.IG.SH    610    1125      all   pairwise     0.1
# 50 1.2406630 3.560530e-04   LS-TT.IG.SH    571    1042 positive   pairwise     0.1
# 51 0.9005673 7.207808e-01   LS-TT.IG.SH     39      83 negative   pairwise     0.1

prwiseTab_FDR05
#           OR          Pval          test NumSig SetSize       ID model_type fdr_cut
# 1  1.0199036  2.969269e-01      LS-Astro   1887    3533      all   pairwise     0.1
# 2  0.5011011  1.000000e+00      LS-Astro    849    2282 positive   pairwise     0.1
# 3  4.5413165 2.101168e-114      LS-Astro   1038    1251 negative   pairwise     0.1
# 4  1.2625743  5.400500e-11       LS-Chol   2081    3533      all   pairwise     0.1
# 5  1.1263701  3.460976e-03       LS-Chol   1290    2282 positive   pairwise     0.1
# 6  1.4998078  3.708477e-12       LS-Chol    791    1251 negative   pairwise     0.1
# 7  1.0914319  7.439879e-03        LS-ChP   1774    3533      all   pairwise     0.1
# 8  0.6390076  1.000000e+00        LS-ChP    869    2282 positive   pairwise     0.1
# 9  2.9150396  1.212752e-69        LS-ChP    905    1251 negative   pairwise     0.1
# 10 0.7941486  1.000000e+00       LS-Endo   1765    3533      all   pairwise     0.1
# 11 0.3857778  1.000000e+00       LS-Endo    766    2282 positive   pairwise     0.1
# 12 3.3681802  2.636660e-78       LS-Endo    999    1251 negative   pairwise     0.1
# 13 1.0382606  1.504260e-01  LS-Ependymal   1750    3533      all   pairwise     0.1
# 14 0.4966355  1.000000e+00  LS-Ependymal    755    2282 positive   pairwise     0.1
# 15 4.3045553 3.785889e-116  LS-Ependymal    995    1251 negative   pairwise     0.1
# 16 1.2470697  1.242622e-09        LS-IoC   2286    3533      all   pairwise     0.1
# 17 1.1309870  3.156287e-03        LS-IoC   1433    2282 positive   pairwise     0.1
# 18 1.4433661  8.885726e-10        LS-IoC    853    1251 negative   pairwise     0.1
# 19 0.6050466  1.000000e+00      LS-Micro   1595    3533      all   pairwise     0.1
# 20 0.2307025  1.000000e+00      LS-Micro    563    2282 positive   pairwise     0.1
# 21 3.8241365  1.276602e-89      LS-Micro   1032    1251 negative   pairwise     0.1
# 22 1.3056766  5.708529e-14         LS-MS   1686    3533      all   pairwise     0.1
# 23 1.3256064  5.909821e-11         LS-MS   1103    2282 positive   pairwise     0.1
# 24 1.2213668  3.129114e-04         LS-MS    583    1251 negative   pairwise     0.1
# 25 0.7895550  1.000000e+00      LS-Mural   1763    3533      all   pairwise     0.1
# 26 0.3616345  1.000000e+00      LS-Mural    737    2282 positive   pairwise     0.1
# 27 3.8763838  6.092780e-93      LS-Mural   1026    1251 negative   pairwise     0.1
# 28 0.9585093  8.860057e-01 LS-Neuroblast   1731    3533      all   pairwise     0.1
# 29 0.4698151  1.000000e+00 LS-Neuroblast    754    2282 positive   pairwise     0.1
# 30 3.7412968  5.287397e-97 LS-Neuroblast    977    1251 negative   pairwise     0.1
# 31 0.8418348  9.999992e-01      LS-Oligo   1974    3533      all   pairwise     0.1
# 32 0.4057361  1.000000e+00      LS-Oligo    890    2282 positive   pairwise     0.1
# 33 4.6023610 2.975224e-101      LS-Oligo   1084    1251 negative   pairwise     0.1
# 34 0.9293898  9.807460e-01        LS-OPC   1826    3533      all   pairwise     0.1
# 35 0.4752699  1.000000e+00        LS-OPC    830    2282 positive   pairwise     0.1
# 36 3.5697189  6.769411e-87        LS-OPC    996    1251 negative   pairwise     0.1
# 37 1.1353371  2.403345e-04       LS-Sept   1511    3533      all   pairwise     0.1
# 38 0.8831922  9.974623e-01       LS-Sept    851    2282 positive   pairwise     0.1
# 39 1.7092135  1.252007e-20       LS-Sept    660    1251 negative   pairwise     0.1
# 40 1.2917757  6.060623e-13        LS-Str   2061    3533      all   pairwise     0.1
# 41 1.3838779  8.512955e-14        LS-Str   1372    2282 positive   pairwise     0.1
# 42 1.1042766  4.606445e-02        LS-Str    689    1251 negative   pairwise     0.1
# 43 1.0635334  4.442401e-02       LS-Thal   1609    3533      all   pairwise     0.1
# 44 0.9344307  9.415181e-01       LS-Thal    973    2282 positive   pairwise     0.1
# 45 1.3206202  8.407784e-07       LS-Thal    636    1251 negative   pairwise     0.1
# 46 1.1947641  4.315168e-07       LS-TNoS   2071    3533      all   pairwise     0.1
# 47 0.9795512  6.900746e-01       LS-TNoS   1238    2282 positive   pairwise     0.1
# 48 1.6807665  2.153656e-18       LS-TNoS    833    1251 negative   pairwise     0.1
# 49 1.2304863  3.528960e-09   LS-TT.IG.SH   1915    3533      all   pairwise     0.1
# 50 1.4081411  2.677881e-15   LS-TT.IG.SH   1312    2282 positive   pairwise     0.1
# 51 0.9435630  8.497005e-01   LS-TT.IG.SH    603    1251 negative   pairwise     0.1


######################################################################
#### Plot using as base gene_set_enrichment_plot from spatialLIBD ####
######################################################################

# gene_set_enrichment_plot_mod <- function(
#         enrichment,
#         cols = colorRamp2(
#             c(0, 12),
#             c("white", "red")),
#     path_to_plot,
#     plot_name) {
#     enrichment$log10_P_thresh <- round(-log10(enrichment$Pval), 2)
#     enrichment$log10_P_thresh[which(enrichment$log10_P_thresh > 12)] <- 12
#     enrichment$OR_char <- as.character(round(enrichment$OR, 2))
#     enrichment$OR_char[enrichment$log10_P_thresh < 3] <- ""
#     make_wide <- function(var = "OR_char") {
#         res <- reshape(enrichment,
#             idvar = "ID", timevar = "test",
#             direction = "wide", drop = colnames(enrichment)[!colnames(enrichment) %in%
#                 c("ID", "test", var)], sep = "_mypattern_"
#         )[, -1, drop = FALSE]
#         colnames(res) <- gsub(".*_mypattern_", "", colnames(res))
#         rownames(res) <- unique(enrichment$ID)
#         res <- res[, levels(as.factor(enrichment$test))]
#         t(res)
#     }
#     wide_or <- make_wide("OR_char")
#     wide_p <- make_wide("log10_P_thresh")
#
#     pdf(path_to_plot, height = 4, width = 6)
#     plot_gs <- Heatmap(wide_p,
#         row_title_side = "left",
#         rect_gp = gpar(col = "white", lwd = 2),
#         cluster_rows = FALSE,
#         cluster_columns = FALSE,
#         cell_fun = function(j, i, x, y, width, height, fill) {
#             grid.text(wide_or[i, j], x, y, gp = gpar(fontsize = 10))
#         },
#         col = cols,
#         column_title = plot_name,
#         heatmap_legend_param = list(
#             title = "-log10(p-val)",
#             at = c(0, 2, 4, 6, 8, 10, 12)
#         )
#     )
#     print(plot_gs)
#     dev.off()
# }
#
# gene_set_enrichment_plot_mod(enrichment = enrichTab_FDR01, path_to_plot = here("snRNAseq_mouse/", "plots/", "04_clinical_set_enrichment/", "Gene_set_enrichment_FDR01.pdf"), plot_name = "Enrichment FDR < 0.01")
# gene_set_enrichment_plot_mod(enrichment = enrichTab_FDR05, path_to_plot = here("snRNAseq_mouse/", "plots/", "04_clinical_set_enrichment/", "Gene_set_enrichment_FDR05.pdf"), plot_name = "Enrichment FDR < 0.05")
#
# gene_set_enrichment_plot_mod(enrichment = prwiseTab_FDR01, path_to_plot = here("snRNAseq_mouse/", "plots/", "04_clinical_set_enrichment/", "Gene_set_pairwise_FDR01.pdf"), plot_name = "Pairwise FDR < 0.01")
# gene_set_enrichment_plot_mod(enrichment = prwiseTab_FDR05, path_to_plot = here("snRNAseq_mouse/", "plots/", "04_clinical_set_enrichment/", "Gene_set_pairwise_FDR05.pdf"), plot_name = "Pairwise FDR < 0.05")


#####################################################
#### Plot using gene_set_enrichment_plot_complex ####
#####################################################

source(
    here(
        "snRNAseq_mouse",
        "code",
        "04_clinical_set_enrichment",
        "gene_set_enrichment_plot_complex.R"
        )
    )

use_gsepc<-function(modeling_results, gene_list, enrichTab, path_to_plot){
    gene_enrichment_count <- get_gene_enrichment_count(model_results = modeling_results, bayes_anno = NULL)
    gene_list_count <- get_gene_list_count(gene_list)

    gse_plot<-gene_set_enrichment_plot_complex(
        enrichment = enrichTab,
        gene_count_col = gene_list_count,
        gene_count_row = gene_enrichment_count,
        anno_title_col = "DE\nGenes",
        anno_title_row = "Cluster\nGenes"
        )

    pdf(path_to_plot, height = 4, width = 6)
    print(gse_plot)
    dev.off()
}


use_gsepc(
    modeling_results = modeling_results,
    gene_list = gene_list_FDR01,
    enrichTab = enrichTab_FDR01,
    path_to_plot = here(
        "snRNAseq_mouse",
        "plots",
        "04_clinical_set_enrichment",
        "Gene_set_enrichment_FDR01.pdf"
    )
)

use_gsepc(
    modeling_results = modeling_results,
    gene_list = gene_list_FDR05,
    enrichTab = enrichTab_FDR05,
    path_to_plot = here(
        "snRNAseq_mouse",
        "plots",
        "04_clinical_set_enrichment",
        "Gene_set_enrichment_FDR05.pdf"
    )
)


#####################################
#### Reproducibility information ####
#####################################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
