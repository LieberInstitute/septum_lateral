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
# 4.24 MB
lobstr::obj_size(gene_list_FDR05)
# 360.60 kB
lobstr::obj_size(gene_list_FDR01)
# 114.46 kB

class(modeling_results)
# [1] "list"
names(modeling_results)
# [1] "enrichment" "pairwise"
dim(modeling_results$enrichment)
# [1] 5312   56
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
enrichTab_FDR05 <- gene_set_enrichment(gene_list = gene_list_FDR05, modeling_results = modeling_results, model_type = "enrichment")
enrichTab_FDR05
#             OR         Pval       test NumSig SetSize       ID model_type fdr_cut
# 1   1.54034364 9.757121e-02      Astro    186    1360      all enrichment     0.1
# 2   0.10027807 1.000000e+00      Astro     17     602 positive enrichment     0.1
# 3  16.84447941 5.657800e-09      Astro    169     758 negative enrichment     0.1
# 4   1.72459518 1.703007e-04       Chol    229    1360      all enrichment     0.1
# 5   0.16480026 1.000000e+00       Chol     17     602 positive enrichment     0.1
# 6   4.39771152 1.671269e-16       Chol    212     758 negative enrichment     0.1
# 7   1.51283336 1.038501e-02        ChP    306    1360      all enrichment     0.1
# 8   0.25063180 1.000000e+00        ChP     48     602 positive enrichment     0.1
# 9   5.77654398 8.409755e-13        ChP    258     758 negative enrichment     0.1
# 10  1.42540773 7.865324e-02       Endo    214    1360      all enrichment     0.1
# 11  0.16304906 1.000000e+00       Endo     22     602 positive enrichment     0.1
# 12  5.03239275 1.970301e-07       Endo    192     758 negative enrichment     0.1
# 13  1.35077159 8.036386e-02  Ependymal    336    1360      all enrichment     0.1
# 14  0.20157046 1.000000e+00  Ependymal     43     602 positive enrichment     0.1
# 15  4.24106150 3.375277e-08  Ependymal    293     758 negative enrichment     0.1
# 16  1.96931716 5.997075e-04        IoC    229    1360      all enrichment     0.1
# 17  0.07682601 1.000000e+00        IoC      7     602 positive enrichment     0.1
# 18  4.97365076 2.263235e-11        IoC    222     758 negative enrichment     0.1
# 19  1.38691482 3.970426e-02         LS    412    1360      all enrichment     0.1
# 20  0.18608799 1.000000e+00         LS     62     602 positive enrichment     0.1
# 21  9.27867319 3.089488e-15         LS    350     758 negative enrichment     0.1
# 22  1.22752839 3.244889e-01      Micro    139    1360      all enrichment     0.1
# 23  0.15139779 9.999753e-01      Micro     13     602 positive enrichment     0.1
# 24  3.35896035 2.292910e-03      Micro    126     758 negative enrichment     0.1
# 25  1.05525613 4.022036e-01         MS    355    1360      all enrichment     0.1
# 26  0.14581960 1.000000e+00         MS     48     602 positive enrichment     0.1
# 27  5.17511489 1.343914e-12         MS    307     758 negative enrichment     0.1
# 28  0.88811319 7.134830e-01      Mural     73    1360      all enrichment     0.1
# 29  0.18548549 9.999899e-01      Mural     11     602 positive enrichment     0.1
# 30  2.48116853 7.898352e-03      Mural     62     758 negative enrichment     0.1
# 31  0.99987793 5.532114e-01 Neuroblast    167    1360      all enrichment     0.1
# 32  0.12064729 1.000000e+00 Neuroblast     21     602 positive enrichment     0.1
# 33  5.63102004 1.246862e-06 Neuroblast    146     758 negative enrichment     0.1
# 34  1.24468966 3.153025e-01      Oligo    161    1360      all enrichment     0.1
# 35  0.04490009 1.000000e+00      Oligo      9     602 positive enrichment     0.1
# 36 28.68913790 2.649219e-08      Oligo    152     758 negative enrichment     0.1
# 37  1.73028203 2.418160e-03        OPC    217    1360      all enrichment     0.1
# 38  0.10834999 1.000000e+00        OPC     15     602 positive enrichment     0.1
# 39  8.22328251 1.816495e-16        OPC    202     758 negative enrichment     0.1
# 40  1.02154597 4.924148e-01       Sept    342    1360      all enrichment     0.1
# 41  0.16328271 1.000000e+00       Sept     52     602 positive enrichment     0.1
# 42  4.40843978 3.378711e-09       Sept    290     758 negative enrichment     0.1
# 43  1.49430181 1.576533e-02        Str    397    1360      all enrichment     0.1
# 44  0.15289460 1.000000e+00        Str     47     602 positive enrichment     0.1
# 45 10.26658481 8.285197e-17        Str    350     758 negative enrichment     0.1
# 46  0.82651614 9.140142e-01       Thal    417    1360      all enrichment     0.1
# 47  0.11704084 1.000000e+00       Thal     66     602 positive enrichment     0.1
# 48  7.02626554 8.589305e-17       Thal    351     758 negative enrichment     0.1
# 49  2.60034046 2.480758e-06       TNoS    258    1360      all enrichment     0.1
# 50  0.20022890 1.000000e+00       TNoS     25     602 positive enrichment     0.1
# 51 19.14047243 1.445502e-19       TNoS    233     758 negative enrichment     0.1
# 52  0.98590003 5.703323e-01   TT.IG.SH    379    1360      all enrichment     0.1
# 53  0.16775205 1.000000e+00   TT.IG.SH     59     602 positive enrichment     0.1
# 54  4.11247981 6.808819e-13   TT.IG.SH    320     758 negative enrichment     0.1

enrichTab_FDR01 <- gene_set_enrichment(gene_list = gene_list_FDR01, modeling_results = modeling_results, model_type = "enrichment")
enrichTab_FDR01
#            OR        Pval       test NumSig SetSize       ID model_type fdr_cut
# 1  0.53479814 0.930268366      Astro     22     254      all enrichment     0.1
# 2  0.04692113 0.999989938      Astro      2     191 positive enrichment     0.1
# 3         Inf 0.087419500      Astro     20      63 negative enrichment     0.1
# 4  0.55119476 0.979142887       Chol     22     254      all enrichment     0.1
# 5  0.11591987 0.999998834       Chol      4     191 positive enrichment     0.1
# 6  3.14268956 0.039185888       Chol     18      63 negative enrichment     0.1
# 7  0.70852545 0.878580607        ChP     37     254      all enrichment     0.1
# 8  0.22451009 0.999853338        ChP     12     191 positive enrichment     0.1
# 9         Inf 0.007500901        ChP     25      63 negative enrichment     0.1
# 10 0.53269456 0.966729901       Endo     34     254      all enrichment     0.1
# 11 0.12100832 0.999996425       Endo      8     191 positive enrichment     0.1
# 12        Inf 0.018313137       Endo     26      63 negative enrichment     0.1
# 13 0.70341168 0.866078547  Ependymal     42     254      all enrichment     0.1
# 14 0.17903602 0.999871257  Ependymal     11     191 positive enrichment     0.1
# 15        Inf 0.018058780  Ependymal     31      63 negative enrichment     0.1
# 16 0.76688622 0.814619341        IoC     28     254      all enrichment     0.1
# 17 0.06860305 0.999983648        IoC      2     191 positive enrichment     0.1
# 18 3.35105550 0.057707994        IoC     26      63 negative enrichment     0.1
# 19 0.93861635 0.654304231         LS     60     254      all enrichment     0.1
# 20 0.37817319 0.992182707         LS     22     191 positive enrichment     0.1
# 21 5.49379505 0.035964185         LS     38      63 negative enrichment     0.1
# 22 0.94531794 0.682025009      Micro     23     254      all enrichment     0.1
# 23 0.15747098 0.996444503      Micro      4     191 positive enrichment     0.1
# 24        Inf 0.103712825      Micro     19      63 negative enrichment     0.1
# 25 0.54384982 0.983629947         MS     48     254      all enrichment     0.1
# 26 0.17730163 0.999998359         MS     15     191 positive enrichment     0.1
# 27 6.32318788 0.020427878         MS     33      63 negative enrichment     0.1
# 28 0.81247130 0.756903392      Mural     15     254      all enrichment     0.1
# 29 0.31533056 0.982835718      Mural      5     191 positive enrichment     0.1
# 30 3.40647979 0.195164850      Mural     10      63 negative enrichment     0.1
# 31 0.70872651 0.841474952 Neuroblast     25     254      all enrichment     0.1
# 32 0.21933481 0.998254450 Neuroblast      8     191 positive enrichment     0.1
# 33        Inf 0.062983237 Neuroblast     17      63 negative enrichment     0.1
# 34 0.66165674 0.847149522      Oligo     18     254      all enrichment     0.1
# 35 0.07117660 0.999571238      Oligo      2     191 positive enrichment     0.1
# 36        Inf 0.174352434      Oligo     16      63 negative enrichment     0.1
# 37 0.77106259 0.818992476        OPC     28     254      all enrichment     0.1
# 38 0.06526470 0.999996323        OPC      2     191 positive enrichment     0.1
# 39 4.13445649 0.023030031        OPC     26      63 negative enrichment     0.1
# 40 0.39012262 0.999180651       Sept     46     254      all enrichment     0.1
# 41 0.13238321 0.999999928       Sept     14     191 positive enrichment     0.1
# 42 2.51379623 0.145428873       Sept     32      63 negative enrichment     0.1
# 43 0.95143399 0.640370331        Str     57     254      all enrichment     0.1
# 44 0.32899272 0.996083446        Str     18     191 positive enrichment     0.1
# 45 6.03510501 0.024206103        Str     39      63 negative enrichment     0.1
# 46 0.52279926 0.990936359       Thal     53     254      all enrichment     0.1
# 47 0.17288704 0.999999838       Thal     18     191 positive enrichment     0.1
# 48        Inf 0.002506266       Thal     35      63 negative enrichment     0.1
# 49 1.06488095 0.552333980       TNoS     28     254      all enrichment     0.1
# 50 0.22166262 0.997699630       TNoS      6     191 positive enrichment     0.1
# 51        Inf 0.010324449       TNoS     22      63 negative enrichment     0.1
# 52 0.51986024 0.995211478   TT.IG.SH     54     254      all enrichment     0.1
# 53 0.18684794 0.999999936   TT.IG.SH     19     191 positive enrichment     0.1
# 54 8.23628234 0.005043533   TT.IG.SH     35      63 negative enrichment     0.1


###########################################################
#### gene_set_enrichment analysis with pairwise data ######
###########################################################

prwiseTab_FDR05 <- spatialLIBD::gene_set_enrichment(gene_list = gene_list_FDR05, modeling_results = modeling_results, fdr_cut = 0.1, model_type = "pairwise", reverse = FALSE)
prwiseTab_FDR05
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

prwiseTab_FDR01 <- spatialLIBD::gene_set_enrichment(gene_list = gene_list_FDR01, modeling_results = modeling_results, fdr_cut = 0.1, model_type = "pairwise", reverse = FALSE)
prwiseTab_FDR01
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
