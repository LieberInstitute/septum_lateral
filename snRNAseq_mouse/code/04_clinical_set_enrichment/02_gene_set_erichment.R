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
# 4.30 MB
lobstr::obj_size(gene_list_FDR05)
# 360.60 kB
lobstr::obj_size(gene_list_FDR01)
# 114.46 kB

class(modeling_results)
# [1] "list"
names(modeling_results)
# [1] "enrichment" "pairwise"
dim(modeling_results$enrichment)
# [1] 5427   56
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
enrichTab_glFDR05_ctFDR05
#            OR         Pval       test NumSig SetSize       ID model_type fdr_cut
# 1  1.24565962 2.151620e-02      Astro    143    1382      all enrichment    0.05
# 2  2.37748438 1.115005e-11      Astro    103     605 positive enrichment    0.05
# 3  0.51165136 9.999931e-01      Astro     40     777 negative enrichment    0.05
# 4  1.28964309 2.752548e-03       Chol    209    1382      all enrichment    0.05
# 5  0.26485393 1.000000e+00       Chol     25     605 positive enrichment    0.05
# 6  2.48529246 1.778859e-19       Chol    184     777 negative enrichment    0.05
# 7  0.69213816 9.999997e-01        ChP    285    1382      all enrichment    0.05
# 8  1.42012550 1.305363e-04        ChP    193     605 positive enrichment    0.05
# 9  0.34726565 1.000000e+00        ChP     92     777 negative enrichment    0.05
# 10 1.31414319 4.901460e-04       Endo    266    1382      all enrichment    0.05
# 11 2.89073218 5.013470e-26       Endo    197     605 positive enrichment    0.05
# 12 0.45659874 1.000000e+00       Endo     69     777 negative enrichment    0.05
# 13 0.67725116 1.000000e+00  Ependymal    291    1382      all enrichment    0.05
# 14 1.48553540 1.458011e-05  Ependymal    204     605 positive enrichment    0.05
# 15 0.30923474 1.000000e+00  Ependymal     87     777 negative enrichment    0.05
# 16 1.76035052 8.795986e-10        IoC    217    1382      all enrichment    0.05
# 17 0.23487052 1.000000e+00        IoC     19     605 positive enrichment    0.05
# 18 3.57372093 3.232524e-35        IoC    198     777 negative enrichment    0.05
# 19 1.06459941 1.754786e-01         LS    504    1382      all enrichment    0.05
# 20 0.07798136 1.000000e+00         LS     29     605 positive enrichment    0.05
# 21 3.48435098 1.784204e-56         LS    475     777 negative enrichment    0.05
# 22 2.08516898 1.350182e-12      Micro    180    1382      all enrichment    0.05
# 23 5.66303486 2.177373e-48      Micro    161     605 positive enrichment    0.05
# 24 0.24477976 1.000000e+00      Micro     19     777 negative enrichment    0.05
# 25 1.07944476 1.310706e-01         MS    464    1382      all enrichment    0.05
# 26 0.10054995 1.000000e+00         MS     32     605 positive enrichment    0.05
# 27 3.15173632 8.120123e-48         MS    432     777 negative enrichment    0.05
# 28 1.39191930 1.161113e-02      Mural     80    1382      all enrichment    0.05
# 29 2.93137389 5.113218e-11      Mural     64     605 positive enrichment    0.05
# 30 0.39504744 9.999826e-01      Mural     16     777 negative enrichment    0.05
# 31 0.72256806 9.992409e-01 Neuroblast    122    1382      all enrichment    0.05
# 32 0.53887065 9.999678e-01 Neuroblast     40     605 positive enrichment    0.05
# 33 0.94119173 7.040582e-01 Neuroblast     82     777 negative enrichment    0.05
# 34 1.63624166 1.022488e-05      Oligo    135    1382      all enrichment    0.05
# 35 3.65853377 1.020910e-22      Oligo    110     605 positive enrichment    0.05
# 36 0.39502078 9.999998e-01      Oligo     25     777 negative enrichment    0.05
# 37 1.28345677 4.661934e-03        OPC    188    1382      all enrichment    0.05
# 38 1.18351157 1.072731e-01        OPC     80     605 positive enrichment    0.05
# 39 1.27657340 1.973523e-02        OPC    108     777 negative enrichment    0.05
# 40 1.21133604 2.926342e-03       Sept    422    1382      all enrichment    0.05
# 41 0.06919295 1.000000e+00       Sept     18     605 positive enrichment    0.05
# 42 3.51513139 3.375648e-55       Sept    404     777 negative enrichment    0.05
# 43 1.52810942 1.384467e-10        Str    499    1382      all enrichment    0.05
# 44 0.06335258 1.000000e+00        Str     18     605 positive enrichment    0.05
# 45 5.18051666 7.990711e-94        Str    481     777 negative enrichment    0.05
# 46 1.16389814 1.049393e-02       Thal    519    1382      all enrichment    0.05
# 47 0.11622860 1.000000e+00       Thal     41     605 positive enrichment    0.05
# 48 3.63910992 3.376953e-60       Thal    478     777 negative enrichment    0.05
# 49 1.54325244 9.136557e-08       TNoS    276    1382      all enrichment    0.05
# 50 0.14071643 1.000000e+00       TNoS     17     605 positive enrichment    0.05
# 51 3.50763229 2.606076e-42       TNoS    259     777 negative enrichment    0.05
# 52 1.43293418 2.069199e-08   TT.IG.SH    539    1382      all enrichment    0.05
# 53 0.12117846 1.000000e+00   TT.IG.SH     39     605 positive enrichment    0.05
# 54 4.71469854 1.527651e-84   TT.IG.SH    500     777 negative enrichment    0.05

enrichTab_glFDR05_ctFDR1 <- gene_set_enrichment(gene_list = gene_list_FDR05, modeling_results = modeling_results, model_type = "enrichment", fdr_cut = 0.1)
enrichTab_glFDR05_ctFDR1
#            OR         Pval       test NumSig SetSize       ID model_type fdr_cut
# 1  1.23386582 2.637220e-02      Astro    143    1382      all enrichment     0.1
# 2  2.35744194 1.656153e-11      Astro    103     605 positive enrichment     0.1
# 3  0.50787085 9.999946e-01      Astro     40     777 negative enrichment     0.1
# 4  1.24618697 6.606005e-03       Chol    223    1382      all enrichment     0.1
# 5  0.24853620 1.000000e+00       Chol     26     605 positive enrichment     0.1
# 6  2.44532737 6.620177e-20       Chol    197     777 negative enrichment     0.1
# 7  0.66288040 1.000000e+00        ChP    300    1382      all enrichment     0.1
# 8  1.30238321 2.787095e-03        ChP    196     605 positive enrichment     0.1
# 9  0.36285927 1.000000e+00        ChP    104     777 negative enrichment     0.1
# 10 1.29542990 8.395842e-04       Endo    269    1382      all enrichment     0.1
# 11 2.85686723 1.000954e-25       Endo    199     605 positive enrichment     0.1
# 12 0.45241561 1.000000e+00       Endo     70     777 negative enrichment     0.1
# 13 0.67257025 1.000000e+00  Ependymal    295    1382      all enrichment     0.1
# 14 1.45967777 3.115760e-05  Ependymal    205     605 positive enrichment     0.1
# 15 0.31435130 1.000000e+00  Ependymal     90     777 negative enrichment     0.1
# 16 1.79063924 1.402835e-10        IoC    227    1382      all enrichment     0.1
# 17 0.22475959 1.000000e+00        IoC     19     605 positive enrichment     0.1
# 18 3.69035804 4.210083e-38        IoC    208     777 negative enrichment     0.1
# 19 1.05057370 2.326938e-01         LS    506    1382      all enrichment     0.1
# 20 0.07663882 1.000000e+00         LS     29     605 positive enrichment     0.1
# 21 3.45969083 7.044539e-56         LS    477     777 negative enrichment     0.1
# 22 2.04418264 3.028463e-12      Micro    183    1382      all enrichment     0.1
# 23 5.47042023 3.331978e-47      Micro    162     605 positive enrichment     0.1
# 24 0.26383450 1.000000e+00      Micro     21     777 negative enrichment     0.1
# 25 1.03522279 3.104177e-01         MS    474    1382      all enrichment     0.1
# 26 0.10033408 1.000000e+00         MS     34     605 positive enrichment     0.1
# 27 3.06144059 8.483070e-46         MS    440     777 negative enrichment     0.1
# 28 1.34068238 1.978366e-02      Mural     83    1382      all enrichment     0.1
# 29 2.81435994 1.083458e-10      Mural     66     605 positive enrichment     0.1
# 30 0.39373049 9.999904e-01      Mural     17     777 negative enrichment     0.1
# 31 0.73637602 9.987406e-01 Neuroblast    126    1382      all enrichment     0.1
# 32 0.52712790 9.999834e-01 Neuroblast     40     605 positive enrichment     0.1
# 33 0.97787814 5.916910e-01 Neuroblast     86     777 negative enrichment     0.1
# 34 1.62932060 1.189016e-05      Oligo    135    1382      all enrichment     0.1
# 35 3.64454045 1.288082e-22      Oligo    110     605 positive enrichment     0.1
# 36 0.39383782 9.999998e-01      Oligo     25     777 negative enrichment     0.1
# 37 1.24223378 1.138930e-02        OPC    190    1382      all enrichment     0.1
# 38 1.13671687 1.744433e-01        OPC     80     605 positive enrichment     0.1
# 39 1.25514952 2.656294e-02        OPC    110     777 negative enrichment     0.1
# 40 1.16802398 1.253143e-02       Sept    426    1382      all enrichment     0.1
# 41 0.06631181 1.000000e+00       Sept     18     605 positive enrichment     0.1
# 42 3.42317812 2.496889e-53       Sept    408     777 negative enrichment     0.1
# 43 1.52910140 1.162902e-10        Str    504    1382      all enrichment     0.1
# 44 0.06598492 1.000000e+00        Str     19     605 positive enrichment     0.1
# 45 5.21463501 1.267507e-94        Str    485     777 negative enrichment     0.1
# 46 1.12389898 3.738446e-02       Thal    531    1382      all enrichment     0.1
# 47 0.11747029 1.000000e+00       Thal     44     605 positive enrichment     0.1
# 48 3.56789667 1.769920e-58       Thal    487     777 negative enrichment     0.1
# 49 1.49741742 4.595370e-07       TNoS    281    1382      all enrichment     0.1
# 50 0.13455194 1.000000e+00       TNoS     17     605 positive enrichment     0.1
# 51 3.43321095 1.052338e-41       TNoS    264     777 negative enrichment     0.1
# 52 1.41764252 4.361350e-08   TT.IG.SH    551    1382      all enrichment     0.1
# 53 0.11892101 1.000000e+00   TT.IG.SH     40     605 positive enrichment     0.1
# 54 4.79855221 2.866402e-86   TT.IG.SH    511     777 negative enrichment     0.1

enrichTab_glFDR01_ctFDR05 <- gene_set_enrichment(gene_list = gene_list_FDR01, modeling_results = modeling_results, model_type = "enrichment", fdr_cut = 0.05)
enrichTab_glFDR01_ctFDR05
#             OR         Pval       test NumSig SetSize       ID model_type fdr_cut
# 1   1.04495705 4.538821e-01      Astro     24     258      all enrichment    0.05
# 2   1.25970284 1.947541e-01      Astro     21     192 positive enrichment    0.05
# 3   0.48097153 9.432725e-01      Astro      3      66 negative enrichment    0.05
# 4   0.61753034 9.905177e-01       Chol     22     258      all enrichment    0.05
# 5   0.17469295 9.999999e-01       Chol      5     192 positive enrichment    0.05
# 6   2.37559189 3.376388e-03       Chol     17      66 negative enrichment    0.05
# 7   0.61817814 9.990513e-01        ChP     46     258      all enrichment    0.05
# 8   0.80823393 9.025152e-01        ChP     42     192 positive enrichment    0.05
# 9   0.18524218 9.999935e-01        ChP      4      66 negative enrichment    0.05
# 10  2.11078737 5.201970e-07       Endo     73     258      all enrichment    0.05
# 11  3.10183412 4.775736e-12       Endo     70     192 positive enrichment    0.05
# 12  0.24120395 9.993031e-01       Endo      3      66 negative enrichment    0.05
# 13  0.54325925 9.999480e-01  Ependymal     43     258      all enrichment    0.05
# 14  0.70173459 9.813932e-01  Ependymal     39     192 positive enrichment    0.05
# 15  0.17738266 9.999965e-01  Ependymal      4      66 negative enrichment    0.05
# 16  0.93020306 6.663091e-01        IoC     27     258      all enrichment    0.05
# 17  0.12241128 1.000000e+00        IoC      3     192 positive enrichment    0.05
# 18  4.70821552 5.617345e-08        IoC     24      66 negative enrichment    0.05
# 19  0.39171162 1.000000e+00         LS     47     258      all enrichment    0.05
# 20  0.03683643 1.000000e+00         LS      4     192 positive enrichment    0.05
# 21  3.46651137 7.221305e-07         LS     43      66 negative enrichment    0.05
# 22  8.35002945 7.610828e-43      Micro     98     258      all enrichment    0.05
# 23 14.06644460 2.870195e-56      Micro     97     192 positive enrichment    0.05
# 24  0.16792538 9.968561e-01      Micro      1      66 negative enrichment    0.05
# 25  0.38124642 1.000000e+00         MS     41     258      all enrichment    0.05
# 26  0.06436602 1.000000e+00         MS      6     192 positive enrichment    0.05
# 27  2.39164719 3.693040e-04         MS     35      66 negative enrichment    0.05
# 28  1.79604775 1.535170e-02      Mural     20     258      all enrichment    0.05
# 29  2.22045817 3.131909e-03      Mural     18     192 positive enrichment    0.05
# 30  0.64160042 8.172853e-01      Mural      2      66 negative enrichment    0.05
# 31  0.70247435 9.534768e-01 Neuroblast     21     258      all enrichment    0.05
# 32  0.43262421 9.986712e-01 Neuroblast     10     192 positive enrichment    0.05
# 33  1.62019181 1.067801e-01 Neuroblast     11      66 negative enrichment    0.05
# 34  1.49691494 4.296142e-02      Oligo     26     258      all enrichment    0.05
# 35  1.92278581 4.467972e-03      Oligo     24     192 positive enrichment    0.05
# 36  0.40507595 9.544347e-01      Oligo      2      66 negative enrichment    0.05
# 37  0.88510139 7.500794e-01        OPC     27     258      all enrichment    0.05
# 38  0.63667267 9.685785e-01        OPC     15     192 positive enrichment    0.05
# 39  1.70529253 7.489293e-02        OPC     12      66 negative enrichment    0.05
# 40  0.46663073 9.999991e-01       Sept     40     258      all enrichment    0.05
# 41  0.02629090 1.000000e+00       Sept      2     192 positive enrichment    0.05
# 42  3.62173908 2.674831e-07       Sept     38      66 negative enrichment    0.05
# 43  0.56599395 9.999275e-01        Str     50     258      all enrichment    0.05
# 44  0.03646337 1.000000e+00        Str      3     192 positive enrichment    0.05
# 45  6.11231279 1.815393e-12        Str     47      66 negative enrichment    0.05
# 46  0.51312508 9.999985e-01       Thal     57     258      all enrichment    0.05
# 47  0.07703142 1.000000e+00       Thal      8     192 positive enrichment    0.05
# 48  5.47743001 6.687739e-11       Thal     49      66 negative enrichment    0.05
# 49  0.65422746 9.883455e-01       TNoS     28     258      all enrichment    0.05
# 50  0.08354024 1.000000e+00       TNoS      3     192 positive enrichment    0.05
# 51  3.40500755 6.892093e-06       TNoS     25      66 negative enrichment    0.05
# 52  0.45089130 1.000000e+00   TT.IG.SH     48     258      all enrichment    0.05
# 53  0.08448297 1.000000e+00   TT.IG.SH      8     192 positive enrichment    0.05
# 54  3.18185209 3.220426e-06   TT.IG.SH     40      66 negative enrichment    0.05

enrichTab_glFDR01_ctFDR1 <- gene_set_enrichment(gene_list = gene_list_FDR01, modeling_results = modeling_results, model_type = "enrichment", fdr_cut = 0.1)
enrichTab_glFDR01_ctFDR1
#             OR         Pval       test NumSig SetSize       ID model_type fdr_cut
# 1   1.03755221 4.665715e-01      Astro     24     258      all enrichment     0.1
# 2   1.25084259 2.026048e-01      Astro     21     192 positive enrichment     0.1
# 3   0.47770894 9.448835e-01      Astro      3      66 negative enrichment     0.1
# 4   0.58490847 9.961900e-01       Chol     23     258      all enrichment     0.1
# 5   0.15770858 1.000000e+00       Chol      5     192 positive enrichment     0.1
# 6   2.31939984 3.477964e-03       Chol     18      66 negative enrichment     0.1
# 7   0.55817547 9.999243e-01        ChP     46     258      all enrichment     0.1
# 8   0.73026515 9.709733e-01        ChP     42     192 positive enrichment     0.1
# 9   0.16779765 9.999985e-01        ChP      4      66 negative enrichment     0.1
# 10  2.05661613 1.158432e-06       Endo     73     258      all enrichment     0.1
# 11  3.02263073 1.277681e-11       Endo     70     192 positive enrichment     0.1
# 12  0.23544055 9.994404e-01       Endo      3      66 negative enrichment     0.1
# 13  0.53065420 9.999730e-01  Ependymal     43     258      all enrichment     0.1
# 14  0.68550229 9.868026e-01  Ependymal     39     192 positive enrichment     0.1
# 15  0.17338983 9.999975e-01  Ependymal      4      66 negative enrichment     0.1
# 16  0.92879284 6.707176e-01        IoC     28     258      all enrichment     0.1
# 17  0.11731439 1.000000e+00        IoC      3     192 positive enrichment     0.1
# 18  4.81809004 2.425674e-08        IoC     25      66 negative enrichment     0.1
# 19  0.38522393 1.000000e+00         LS     47     258      all enrichment     0.1
# 20  0.03623640 1.000000e+00         LS      4     192 positive enrichment     0.1
# 21  3.41030213 9.988181e-07         LS     43      66 negative enrichment     0.1
# 22  8.03173463 1.143999e-41      Micro     98     258      all enrichment     0.1
# 23 13.52843482 4.943273e-55      Micro     97     192 positive enrichment     0.1
# 24  0.16277862 9.973577e-01      Micro      1      66 negative enrichment     0.1
# 25  0.37856999 1.000000e+00         MS     43     258      all enrichment     0.1
# 26  0.06034002 1.000000e+00         MS      6     192 positive enrichment     0.1
# 27  2.53845329 1.471976e-04         MS     37      66 negative enrichment     0.1
# 28  1.67434670 2.796742e-02      Mural     20     258      all enrichment     0.1
# 29  2.07109541 6.005992e-03      Mural     18     192 positive enrichment     0.1
# 30  0.60097697 8.436194e-01      Mural      2      66 negative enrichment     0.1
# 31  0.72491027 9.406012e-01 Neuroblast     22     258      all enrichment     0.1
# 32  0.42370104 9.989873e-01 Neuroblast     10     192 positive enrichment     0.1
# 33  1.76641642 6.213668e-02 Neuroblast     12      66 negative enrichment     0.1
# 34  1.49245913 4.411887e-02      Oligo     26     258      all enrichment     0.1
# 35  1.91710094 4.616492e-03      Oligo     24     192 positive enrichment     0.1
# 36  0.40394273 9.549262e-01      Oligo      2      66 negative enrichment     0.1
# 37  0.85291845 8.048851e-01        OPC     27     258      all enrichment     0.1
# 38  0.61395386 9.778904e-01        OPC     15     192 positive enrichment     0.1
# 39  1.64488325 9.007088e-02        OPC     12      66 negative enrichment     0.1
# 40  0.46133587 9.999995e-01       Sept     41     258      all enrichment     0.1
# 41  0.02523967 1.000000e+00       Sept      2     192 positive enrichment     0.1
# 42  3.70321886 1.649355e-07       Sept     39      66 negative enrichment     0.1
# 43  0.55719630 9.999534e-01        Str     50     258      all enrichment     0.1
# 44  0.03590904 1.000000e+00        Str      3     192 positive enrichment     0.1
# 45  6.01891081 2.777209e-12        Str     47      66 negative enrichment     0.1
# 46  0.49201513 9.999997e-01       Thal     58     258      all enrichment     0.1
# 47  0.08179356 1.000000e+00       Thal      9     192 positive enrichment     0.1
# 48  5.13902261 3.421831e-10       Thal     49      66 negative enrichment     0.1
# 49  0.65175854 9.898531e-01       TNoS     29     258      all enrichment     0.1
# 50  0.07998573 1.000000e+00       TNoS      3     192 positive enrichment     0.1
# 51  3.47760709 3.831914e-06       TNoS     26      66 negative enrichment     0.1
# 52  0.44225314 1.000000e+00   TT.IG.SH     49     258      all enrichment     0.1
# 53  0.08078098 1.000000e+00   TT.IG.SH      8     192 positive enrichment     0.1
# 54  3.24641631 2.276398e-06   TT.IG.SH     41      66 negative enrichment     0.1

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
