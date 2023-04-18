library("here")
library("sessioninfo")
library("dplyr")
library("purrr")
library("spatialLIBD")



##################### Load objects for gene_set_enrichment ####################

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
#   modeling_results_broad
#   modeling_results
#   gene_list_FDR05
#   gene_list_FDR01

lobstr::obj_size(modeling_results_broad)
# 6.62 MB
lobstr::obj_size(modeling_results)
# 34.57 MB
lobstr::obj_size(gene_list_FDR05)
# 360.60 kB
lobstr::obj_size(gene_list_FDR01)
# 114.46 kB

class(modeling_results_broad)
# [1] "list"
names(modeling_results_broad)
# [1] "enrichment" "pairwise"
dim(modeling_results_broad$enrichment)
# [1] 5436   56
head(names(modeling_results_broad$enrichment), 10)
# [1] "t_stat_Astro"  "p_value_Astro" "fdr_Astro"     "t_stat_Chol"
# [5] "p_value_Chol"  "fdr_Chol"      "t_stat_ChP"    "p_value_ChP"
# [9] "fdr_ChP"       "t_stat_Endo"
dim(modeling_results_broad$'pairwise_ct-LS')
# [1] 5436   53
head(names(modeling_results_broad$'pairwise_ct-LS'), 10)
# [1] "t_stat_Astro-LS"  "p_value_Astro-LS" "fdr_Astro-LS"     "t_stat_Chol-LS"
# [5] "p_value_Chol-LS"  "fdr_Chol-LS"      "t_stat_ChP-LS"    "p_value_ChP-LS"
# [9] "fdr_ChP-LS"       "t_stat_Endo-LS"

class(modeling_results)
# [1] "list"
names(modeling_results)
# [1] "enrichment" "pairwise"
dim(modeling_results$enrichment)
# [1] 5126   53
head(names(modeling_results$enrichment), 10)
#  [1] "t_stat_Chol_Ex.D"  "p_value_Chol_Ex.D" "fdr_Chol_Ex.D"
#  [4] "t_stat_LS_In.C"    "p_value_LS_In.C"   "fdr_LS_In.C"
#  [7] "t_stat_LS_In.D"    "p_value_LS_In.D"   "fdr_LS_In.D"
# [10] "t_stat_LS_In.M"
dim(modeling_results$pairwise)
# [1] 5126  770
head(names(modeling_results$pairwise), 10)
# [1] "t_stat_Chol_Ex.D-LS_In.C"  "p_value_Chol_Ex.D-LS_In.C"
# [3] "fdr_Chol_Ex.D-LS_In.C"     "t_stat_Chol_Ex.D-LS_In.D"
# [5] "p_value_Chol_Ex.D-LS_In.D" "fdr_Chol_Ex.D-LS_In.D"
# [7] "t_stat_Chol_Ex.D-LS_In.M"  "p_value_Chol_Ex.D-LS_In.M"
# [9] "fdr_Chol_Ex.D-LS_In.M"     "t_stat_Chol_Ex.D-LS_In.N"

class(gene_list_FDR01)
# [1] "list"
names(gene_list_FDR01)
# [1] "all"      "positive" "negative"
lengths(gene_list_FDR01)
#  all positive negative
# 1186     1089       97
lengths(gene_list_FDR05)
#  all positive negative
# 3750     2402     1348

###############################################################################



########### gene_set_enrichment analysis with 1vsAll broad data ###############

## Running gene_set_enrichment with two sets of DE genes, one with FDR < 0.05 the other with FDR < 0.01
enrichTab_broad_glFDR05_ctFDR05 <- gene_set_enrichment(gene_list = gene_list_FDR05, modeling_results = modeling_results_broad, model_type = "enrichment", fdr_cut = 0.05)
enrichTab_broad_glFDR05_ctFDR05
#            OR         Pval       test NumSig SetSize       ID model_type fdr_cut
# 1  1.23944911 2.393786e-02      Astro    143    1383      all enrichment    0.05
# 2  2.36887710 1.313460e-11      Astro    103     605 positive enrichment    0.05
# 3  0.50940332 9.999940e-01      Astro     40     778 negative enrichment    0.05
# 4  1.29144385 2.624287e-03       Chol    209    1383      all enrichment    0.05
# 5  0.26542863 1.000000e+00       Chol     25     605 positive enrichment    0.05
# 6  2.48590954 1.722941e-19       Chol    184     778 negative enrichment    0.05
# 7  0.69080639 9.999998e-01        ChP    285    1383      all enrichment    0.05
# 8  1.41891666 1.348204e-04        ChP    193     605 positive enrichment    0.05
# 9  0.34647512 1.000000e+00        ChP     92     778 negative enrichment    0.05
# 10 1.30856565 5.832012e-04       Endo    266    1383      all enrichment    0.05
# 11 2.88240130 6.468584e-26       Endo    197     605 positive enrichment    0.05
# 12 0.45488187 1.000000e+00       Endo     69     778 negative enrichment    0.05
# 13 0.67814176 1.000000e+00  Ependymal    292    1383      all enrichment    0.05
# 14 1.48116978 1.662673e-05  Ependymal    204     605 positive enrichment    0.05
# 15 0.31224039 1.000000e+00  Ependymal     88     778 negative enrichment    0.05
# 16 1.76225500 7.570930e-10        IoC    218    1383      all enrichment    0.05
# 17 0.23400310 1.000000e+00        IoC     19     605 positive enrichment    0.05
# 18 3.57922337 1.949497e-35        IoC    199     778 negative enrichment    0.05
# 19 1.06648293 1.683624e-01         LS    505    1383      all enrichment    0.05
# 20 0.07794962 1.000000e+00         LS     29     605 positive enrichment    0.05
# 21 3.48988441 1.102491e-56         LS    476     778 negative enrichment    0.05
# 22 2.07146704 1.976763e-12      Micro    180    1383      all enrichment    0.05
# 23 5.63312760 3.489287e-48      Micro    161     605 positive enrichment    0.05
# 24 0.24367688 1.000000e+00      Micro     19     778 negative enrichment    0.05
# 25 1.07772121 1.361508e-01         MS    464    1383      all enrichment    0.05
# 26 0.10056934 1.000000e+00         MS     32     605 positive enrichment    0.05
# 27 3.13966285 1.374640e-47         MS    432     778 negative enrichment    0.05
# 28 1.37690438 1.395668e-02      Mural     80    1383      all enrichment    0.05
# 29 2.90474050 6.892085e-11      Mural     64     605 positive enrichment    0.05
# 30 0.39173153 9.999857e-01      Mural     16     778 negative enrichment    0.05
# 31 0.72019727 9.993249e-01 Neuroblast    122    1383      all enrichment    0.05
# 32 0.53783494 9.999696e-01 Neuroblast     40     605 positive enrichment    0.05
# 33 0.93758104 7.146944e-01 Neuroblast     82     778 negative enrichment    0.05
# 34 1.62457414 1.313340e-05      Oligo    135    1383      all enrichment    0.05
# 35 3.63783229 1.420102e-22      Oligo    110     605 positive enrichment    0.05
# 36 0.39287161 9.999998e-01      Oligo     25     778 negative enrichment    0.05
# 37 1.27551006 5.585087e-03        OPC    188    1383      all enrichment    0.05
# 38 1.17874454 1.129359e-01        OPC     80     605 positive enrichment    0.05
# 39 1.26892109 2.224630e-02        OPC    108     778 negative enrichment    0.05
# 40 1.21287294 2.743706e-03       Sept    423    1383      all enrichment    0.05
# 41 0.06910982 1.000000e+00       Sept     18     605 positive enrichment    0.05
# 42 3.51914771 2.197870e-55       Sept    405     778 negative enrichment    0.05
# 43 1.53531824 8.584377e-11        Str    500    1383      all enrichment    0.05
# 44 0.06346818 1.000000e+00        Str     18     605 positive enrichment    0.05
# 45 5.20301530 2.084610e-94        Str    482     778 negative enrichment    0.05
# 46 1.16220093 1.110450e-02       Thal    519    1383      all enrichment    0.05
# 47 0.11627579 1.000000e+00       Thal     41     605 positive enrichment    0.05
# 48 3.62493396 6.161571e-60       Thal    478     778 negative enrichment    0.05
# 49 1.53904899 1.072133e-07       TNoS    276    1383      all enrichment    0.05
# 50 0.14062050 1.000000e+00       TNoS     17     605 positive enrichment    0.05
# 51 3.49398631 4.120942e-42       TNoS    259     778 negative enrichment    0.05
# 52 1.43636952 1.646763e-08   TT.IG.SH    540    1383      all enrichment    0.05
# 53 0.12120697 1.000000e+00   TT.IG.SH     39     605 positive enrichment    0.05
# 54 4.72520867 7.030047e-85   TT.IG.SH    501     778 negative enrichment    0.05

enrichTab_broad_glFDR01_ctFDR05 <- gene_set_enrichment(gene_list = gene_list_FDR01, modeling_results = modeling_results_broad, model_type = "enrichment", fdr_cut = 0.05)
enrichTab_broad_glFDR01_ctFDR05
#             OR         Pval       test NumSig SetSize       ID model_type fdr_cut
# 1   1.04199962 4.589249e-01      Astro     24     258      all enrichment    0.05
# 2   1.25614955 1.978612e-01      Astro     21     192 positive enrichment    0.05
# 3   0.47967555 9.439160e-01      Astro      3      66 negative enrichment    0.05
# 4   0.61876760 9.902613e-01       Chol     22     258      all enrichment    0.05
# 5   0.17503919 9.999999e-01       Chol      5     192 positive enrichment    0.05
# 6   2.38046280 3.317965e-03       Chol     17      66 negative enrichment    0.05
# 7   0.61776773 9.990657e-01        ChP     46     258      all enrichment    0.05
# 8   0.80768247 9.032073e-01        ChP     42     192 positive enrichment    0.05
# 9   0.18512068 9.999935e-01        ChP      4      66 negative enrichment    0.05
# 10  2.10593713 5.581527e-07       Endo     73     258      all enrichment    0.05
# 11  3.09468703 5.207355e-12       Endo     70     192 positive enrichment    0.05
# 12  0.24071028 9.993158e-01       Endo      3      66 negative enrichment    0.05
# 13  0.54188671 9.999515e-01  Ependymal     43     258      all enrichment    0.05
# 14  0.69995763 9.820621e-01  Ependymal     39     192 positive enrichment    0.05
# 15  0.17694444 9.999966e-01  Ependymal      4      66 negative enrichment    0.05
# 16  0.92660006 6.731898e-01        IoC     27     258      all enrichment    0.05
# 17  0.12196097 1.000000e+00        IoC      3     192 positive enrichment    0.05
# 18  4.68988752 5.992276e-08        IoC     24      66 negative enrichment    0.05
# 19  0.39147038 1.000000e+00         LS     47     258      all enrichment    0.05
# 20  0.03681511 1.000000e+00         LS      4     192 positive enrichment    0.05
# 21  3.46411142 7.319815e-07         LS     43      66 negative enrichment    0.05
# 22  8.31505926 1.006369e-42      Micro     98     258      all enrichment    0.05
# 23 14.00641851 3.850388e-56      Micro     97     192 positive enrichment    0.05
# 24  0.16742062 9.969080e-01      Micro      1      66 negative enrichment    0.05
# 25  0.38124054 1.000000e+00         MS     41     258      all enrichment    0.05
# 26  0.06436628 1.000000e+00         MS      6     192 positive enrichment    0.05
# 27  2.39141458 3.697636e-04         MS     35      66 negative enrichment    0.05
# 28  1.78318125 1.635292e-02      Mural     20     258      all enrichment    0.05
# 29  2.20441802 3.352150e-03      Mural     18     192 positive enrichment    0.05
# 30  0.63735978 8.200481e-01      Mural      2      66 negative enrichment    0.05
# 31  0.70112396 9.543637e-01 Neuroblast     21     258      all enrichment    0.05
# 32  0.43181239 9.987032e-01 Neuroblast     10     192 positive enrichment    0.05
# 33  1.61707847 1.077312e-01 Neuroblast     11      66 negative enrichment    0.05
# 34  1.49081339 4.454521e-02      Oligo     26     258      all enrichment    0.05
# 35  1.91497612 4.671830e-03      Oligo     24     192 positive enrichment    0.05
# 36  0.40354394 9.550993e-01      Oligo      2      66 negative enrichment    0.05
# 37  0.88187568 7.558013e-01        OPC     27     258      all enrichment    0.05
# 38  0.63440511 9.696153e-01        OPC     15     192 positive enrichment    0.05
# 39  1.69920297 7.629654e-02        OPC     12      66 negative enrichment    0.05
# 40  0.46598207 9.999992e-01       Sept     40     258      all enrichment    0.05
# 41  0.02625590 1.000000e+00       Sept      2     192 positive enrichment    0.05
# 42  3.61648461 2.752880e-07       Sept     38      66 negative enrichment    0.05
# 43  0.56687392 9.999243e-01        Str     50     258      all enrichment    0.05
# 44  0.03652034 1.000000e+00        Str      3     192 positive enrichment    0.05
# 45  6.12115417 1.743086e-12        Str     47      66 negative enrichment    0.05
# 46  0.51321304 9.999985e-01       Thal     57     258      all enrichment    0.05
# 47  0.07704715 1.000000e+00       Thal      8     192 positive enrichment    0.05
# 48  5.47789527 6.669855e-11       Thal     49      66 negative enrichment    0.05
# 49  0.65366633 9.884896e-01       TNoS     28     258      all enrichment    0.05
# 50  0.08347345 1.000000e+00       TNoS      3     192 positive enrichment    0.05
# 51  3.40189384 6.985695e-06       TNoS     25      66 negative enrichment    0.05
# 52  0.45090207 1.000000e+00   TT.IG.SH     48     258      all enrichment    0.05
# 53  0.08448714 1.000000e+00   TT.IG.SH      8     192 positive enrichment    0.05
# 54  3.18167226 3.223033e-06   TT.IG.SH     40      66 negative enrichment    0.05

###############################################################################



############# gene_set_enrichment analysis with 1vs1 broad data ###############

modeling_results_broad$pairwise <- modeling_results_broad$'pairwise_LS-ct'

prwiseTab_LS_ct_broad_glFDR05_ctFDR05 <- spatialLIBD::gene_set_enrichment(gene_list = gene_list_FDR05, modeling_results = modeling_results_broad, fdr_cut = 0.05, model_type = "pairwise", reverse = FALSE)
prwiseTab_LS_ct_broad_glFDR05_ctFDR05
#           OR         Pval          test NumSig SetSize       ID model_type fdr_cut
# 1  1.6785369 3.420903e-06      LS-Astro    548     656      all   pairwise    0.05
# 2  0.1775944 1.000000e+00      LS-Astro     33      84 positive   pairwise    0.05
# 3  3.2264532 9.465992e-19      LS-Astro    515     572 negative   pairwise    0.05
# 4  2.0294245 1.445336e-12       LS-Chol    226     656      all   pairwise    0.05
# 5  0.3706590 9.994897e-01       LS-Chol      9      84 positive   pairwise    0.05
# 6  2.4190778 1.476913e-17       LS-Chol    217     572 negative   pairwise    0.05
# 7  3.3461145 2.531124e-39        LS-ChP    384     656      all   pairwise    0.05
# 8  0.3933447 9.999027e-01        LS-ChP     16      84 positive   pairwise    0.05
# 9  4.3646166 5.295388e-52        LS-ChP    368     572 negative   pairwise    0.05
# 10 1.9369395 1.534973e-11       LS-Endo    495     656      all   pairwise    0.05
# 11 0.1596054 1.000000e+00       LS-Endo     20      84 positive   pairwise    0.05
# 12 3.2884310 4.541545e-27       LS-Endo    475     572 negative   pairwise    0.05
# 13 2.2204905 2.510115e-16  LS-Ependymal    494     656      all   pairwise    0.05
# 14 0.3616411 9.999984e-01  LS-Ependymal     32      84 positive   pairwise    0.05
# 15 3.1603361 5.400609e-27  LS-Ependymal    462     572 negative   pairwise    0.05
# 16 1.0658084 2.556219e-01        LS-IoC    368     656      all   pairwise    0.05
# 17 0.8157253 8.488034e-01        LS-IoC     42      84 positive   pairwise    0.05
# 18 1.1135504 1.396903e-01        LS-IoC    326     572 negative   pairwise    0.05
# 19 2.1651608 3.317616e-15      LS-Micro    498     656      all   pairwise    0.05
# 20 0.2067310 1.000000e+00      LS-Micro     23      84 positive   pairwise    0.05
# 21 3.5483216 9.688183e-31      LS-Micro    475     572 negative   pairwise    0.05
# 22 1.4226325 4.832903e-04         LS-MS    177     656      all   pairwise    0.05
# 23 0.7556116 8.673036e-01         LS-MS     15      84 positive   pairwise    0.05
# 24 1.5316499 5.908507e-05         LS-MS    162     572 negative   pairwise    0.05
# 25 1.9797000 1.013670e-11      LS-Mural    509     656      all   pairwise    0.05
# 26 0.1637282 1.000000e+00      LS-Mural     22      84 positive   pairwise    0.05
# 27 3.4964706 7.343321e-28      LS-Mural    487     572 negative   pairwise    0.05
# 28 2.0369336 1.019080e-13 LS-Neuroblast    479     656      all   pairwise    0.05
# 29 0.3427011 9.999995e-01 LS-Neuroblast     30      84 positive   pairwise    0.05
# 30 2.8462031 1.412188e-23 LS-Neuroblast    449     572 negative   pairwise    0.05
# 31 1.7066244 4.282795e-06      LS-Oligo    560     656      all   pairwise    0.05
# 32 0.2451727 1.000000e+00      LS-Oligo     42      84 positive   pairwise    0.05
# 33 2.9871177 1.121422e-15      LS-Oligo    518     572 negative   pairwise    0.05
# 34 1.7410766 4.790348e-09        LS-OPC    469     656      all   pairwise    0.05
# 35 0.2327456 1.000000e+00        LS-OPC     24      84 positive   pairwise    0.05
# 36 2.5594973 1.975488e-19        LS-OPC    445     572 negative   pairwise    0.05
# 37 1.2986478 3.339801e-03       LS-Sept    243     656      all   pairwise    0.05
# 38 0.5536633 9.922595e-01       LS-Sept     18      84 positive   pairwise    0.05
# 39 1.4579993 8.075664e-05       LS-Sept    225     572 negative   pairwise    0.05
# 40 0.6042620 1.000000e+00        LS-Str    310     656      all   pairwise    0.05
# 41 1.0749729 4.187636e-01        LS-Str     49      84 positive   pairwise    0.05
# 42 0.5673232 1.000000e+00        LS-Str    261     572 negative   pairwise    0.05
# 43 1.1565722 7.524598e-02       LS-Thal    207     656      all   pairwise    0.05
# 44 0.8005248 8.399262e-01       LS-Thal     21      84 positive   pairwise    0.05
# 45 1.2173112 3.047749e-02       LS-Thal    186     572 negative   pairwise    0.05
# 46 1.2977328 2.216868e-03       LS-TNoS    339     656      all   pairwise    0.05
# 47 0.4712289 9.995901e-01       LS-TNoS     25      84 positive   pairwise    0.05
# 48 1.5147567 7.176966e-06       LS-TNoS    314     572 negative   pairwise    0.05
# 49 0.7270523 9.997642e-01   LS-TT.IG.SH    237     656      all   pairwise    0.05
# 50 1.1525774 2.982904e-01   LS-TT.IG.SH     38      84 positive   pairwise    0.05
# 51 0.6841056 9.999607e-01   LS-TT.IG.SH    199     572 negative   pairwise    0.05

prwiseTab_LS_ct_broad_glFDR01_ctFDR05 <- spatialLIBD::gene_set_enrichment(gene_list = gene_list_FDR01, modeling_results = modeling_results_broad, fdr_cut = 0.05, model_type = "pairwise", reverse = FALSE)
# Warning message:
# In spatialLIBD::gene_set_enrichment(gene_list = gene_list_FDR01,  :
#   Gene list with n < 25 may have insufficent power for enrichment analysis: all ,positive ,negative
prwiseTab_LS_ct_broad_glFDR01_ctFDR05
#            OR         Pval          test NumSig SetSize       ID model_type fdr_cut
# 1  1.14092543 3.957399e-01      LS-Astro     54      68      all   pairwise    0.05
# 2  0.15848535 9.999653e-01      LS-Astro      6      17 positive   pairwise    0.05
# 3  4.80552376 1.152620e-03      LS-Astro     48      51 negative   pairwise    0.05
# 4  2.26926643 1.140716e-03       LS-Chol     28      68      all   pairwise    0.05
# 5  0.00000000 1.000000e+00       LS-Chol      0      17 positive   pairwise    0.05
# 6  3.97917812 1.651032e-06       LS-Chol     28      51 negative   pairwise    0.05
# 7  3.03645186 6.948114e-06        LS-ChP     43      68      all   pairwise    0.05
# 8  0.10623399 9.996081e-01        LS-ChP      1      17 positive   pairwise    0.05
# 9  8.30755641 2.195998e-11        LS-ChP     42      51 negative   pairwise    0.05
# 10 1.40892581 1.280565e-01       LS-Endo     49      68      all   pairwise    0.05
# 11 0.11463742 9.999912e-01       LS-Endo      3      17 positive   pairwise    0.05
# 12 5.09230325 3.016477e-05       LS-Endo     46      51 negative   pairwise    0.05
# 13 2.00644967 8.353149e-03  LS-Ependymal     52      68      all   pairwise    0.05
# 14 0.32915342 9.938740e-01  LS-Ependymal      6      17 positive   pairwise    0.05
# 15 5.72005776 6.040400e-06  LS-Ependymal     46      51 negative   pairwise    0.05
# 16 0.76890454 8.838678e-01        LS-IoC     33      68      all   pairwise    0.05
# 17 0.17433785 9.997176e-01        LS-IoC      3      17 positive   pairwise    0.05
# 18 1.17645741 3.371859e-01        LS-IoC     30      51 negative   pairwise    0.05
# 19 1.21195579 2.745757e-01      LS-Micro     46      68      all   pairwise    0.05
# 20 0.07587499 9.999989e-01      LS-Micro      2      17 positive   pairwise    0.05
# 21 3.69582327 2.413810e-04      LS-Micro     44      51 negative   pairwise    0.05
# 22 1.70197827 3.295297e-02         LS-MS     22      68      all   pairwise    0.05
# 23 0.46553366 9.186818e-01         LS-MS      2      17 positive   pairwise    0.05
# 24 2.30458207 4.213162e-03         LS-MS     20      51 negative   pairwise    0.05
# 25 1.18176912 3.160668e-01      LS-Mural     48      68      all   pairwise    0.05
# 26 0.10365073 9.999965e-01      LS-Mural      3      17 positive   pairwise    0.05
# 27 3.74636564 4.384457e-04      LS-Mural     45      51 negative   pairwise    0.05
# 28 1.34515653 1.565456e-01 LS-Neuroblast     46      68      all   pairwise    0.05
# 29 0.13531651 9.999622e-01 LS-Neuroblast      3      17 positive   pairwise    0.05
# 30 3.49840128 2.480291e-04 LS-Neuroblast     43      51 negative   pairwise    0.05
# 31 0.77508051 8.529769e-01      LS-Oligo     51      68      all   pairwise    0.05
# 32 0.13958188 9.999880e-01      LS-Oligo      6      17 positive   pairwise    0.05
# 33 1.97167161 7.355492e-02      LS-Oligo     45      51 negative   pairwise    0.05
# 34 0.86774997 7.580277e-01        LS-OPC     40      68      all   pairwise    0.05
# 35 0.12911640 9.999746e-01        LS-OPC      3      17 positive   pairwise    0.05
# 36 1.62507660 7.808967e-02        LS-OPC     37      51 negative   pairwise    0.05
# 37 0.92027507 6.669874e-01       LS-Sept     21      68      all   pairwise    0.05
# 38 0.00000000 1.000000e+00       LS-Sept      0      17 positive   pairwise    0.05
# 39 1.45571536 1.234461e-01       LS-Sept     21      51 negative   pairwise    0.05
# 40 0.56125915 9.932096e-01        LS-Str     29      68      all   pairwise    0.05
# 41 0.53421221 9.369888e-01        LS-Str      7      17 positive   pairwise    0.05
# 42 0.57503406 9.819407e-01        LS-Str     22      51 negative   pairwise    0.05
# 43 0.80158209 8.197779e-01       LS-Thal     17      68      all   pairwise    0.05
# 44 0.51629378 9.136335e-01       LS-Thal      3      17 positive   pairwise    0.05
# 45 0.91323878 6.644294e-01       LS-Thal     14      51 negative   pairwise    0.05
# 46 1.86584438 8.539124e-03       LS-TNoS     42      68      all   pairwise    0.05
# 47 0.15031227 9.996619e-01       LS-TNoS      2      17 positive   pairwise    0.05
# 48 4.23803830 3.140119e-06       LS-TNoS     40      51 negative   pairwise    0.05
# 49 1.16836602 3.045756e-01   LS-TT.IG.SH     31      68      all   pairwise    0.05
# 50 0.75635595 7.853078e-01   LS-TT.IG.SH      6      17 positive   pairwise    0.05
# 51 1.34312875 1.828590e-01   LS-TT.IG.SH     25      51 negative   pairwise    0.05


modeling_results_broad$pairwise <- modeling_results_broad$'pairwise_ct-LS'

prwiseTab_ct_LS_broad_glFDR05_ctFDR05 <- spatialLIBD::gene_set_enrichment(gene_list = gene_list_FDR05, modeling_results = modeling_results_broad, fdr_cut = 0.05, model_type = "pairwise", reverse = FALSE)
prwiseTab_ct_LS_broad_glFDR05_ctFDR05
#           OR         Pval          test NumSig SetSize       ID model_type fdr_cut
# 1  1.2366735 2.836170e-02      Astro-LS    135    1383      all   pairwise    0.05
# 2  2.6008860 1.155442e-13      Astro-LS    104     605 positive   pairwise    0.05
# 3  0.4080870 1.000000e+00      Astro-LS     31     778 negative   pairwise    0.05
# 4  0.9376364 7.307271e-01       Chol-LS    110    1383      all   pairwise    0.05
# 5  0.4843298 9.999684e-01       Chol-LS     27     605 positive   pairwise    0.05
# 6  1.3880134 7.467919e-03       Chol-LS     83     778 negative   pairwise    0.05
# 7  0.8210546 9.966265e-01        ChP-LS    302    1383      all   pairwise    0.05
# 8  1.8659468 1.875212e-11        ChP-LS    217     605 positive   pairwise    0.05
# 9  0.3359293 1.000000e+00        ChP-LS     85     778 negative   pairwise    0.05
# 10 1.4391333 4.333630e-06       Endo-LS    281    1383      all   pairwise    0.05
# 11 3.2681800 1.077878e-32       Endo-LS    211     605 positive   pairwise    0.05
# 12 0.4621250 1.000000e+00       Endo-LS     70     778 negative   pairwise    0.05
# 13 0.7610177 9.999167e-01  Ependymal-LS    301    1383      all   pairwise    0.05
# 14 1.8410198 3.310315e-11  Ependymal-LS    223     605 positive   pairwise    0.05
# 15 0.2854387 1.000000e+00  Ependymal-LS     78     778 negative   pairwise    0.05
# 16 1.3389502 3.596622e-03        IoC-LS    145    1383      all   pairwise    0.05
# 17 0.5390196 9.998236e-01        IoC-LS     31     605 positive   pairwise    0.05
# 18 2.0680889 1.444708e-09        IoC-LS    114     778 negative   pairwise    0.05
# 19 2.2399606 3.700237e-15      Micro-LS    190    1383      all   pairwise    0.05
# 20 6.1383385 3.458968e-54      Micro-LS    170     605 positive   pairwise    0.05
# 21 0.2536136 1.000000e+00      Micro-LS     20     778 negative   pairwise    0.05
# 22 1.1301099 9.646736e-02         MS-LS    195    1383      all   pairwise    0.05
# 23 0.3073806 1.000000e+00         MS-LS     29     605 positive   pairwise    0.05
# 24 2.0552247 1.875723e-12         MS-LS    166     778 negative   pairwise    0.05
# 25 1.5421628 8.353388e-04      Mural-LS     92    1383      all   pairwise    0.05
# 26 3.1431560 2.814878e-13      Mural-LS     72     605 positive   pairwise    0.05
# 27 0.4633215 9.998944e-01      Mural-LS     20     778 negative   pairwise    0.05
# 28 0.6924054 9.996007e-01 Neuroblast-LS    103    1383      all   pairwise    0.05
# 29 0.8262016 9.058688e-01 Neuroblast-LS     50     605 positive   pairwise    0.05
# 30 0.6483822 9.989786e-01 Neuroblast-LS     53     778 negative   pairwise    0.05
# 31 1.6614799 7.250853e-06      Oligo-LS    132    1383      all   pairwise    0.05
# 32 3.9605396 3.260983e-25      Oligo-LS    112     605 positive   pairwise    0.05
# 33 0.3208398 1.000000e+00      Oligo-LS     20     778 negative   pairwise    0.05
# 34 1.0880881 2.082749e-01        OPC-LS    160    1383      all   pairwise    0.05
# 35 1.5896940 1.514745e-04        OPC-LS     94     605 positive   pairwise    0.05
# 36 0.7235699 9.938430e-01        OPC-LS     66     778 negative   pairwise    0.05
# 37 1.1158278 1.809432e-01       Sept-LS    117    1383      all   pairwise    0.05
# 38 0.3919631 9.999991e-01       Sept-LS     21     605 positive   pairwise    0.05
# 39 1.8397793 1.438079e-06       Sept-LS     96     778 negative   pairwise    0.05
# 40 1.7166759 1.754274e-11        Str-LS    299    1383      all   pairwise    0.05
# 41 0.2607547 1.000000e+00        Str-LS     31     605 positive   pairwise    0.05
# 42 3.6081765 5.374774e-45        Str-LS    268     778 negative   pairwise    0.05
# 43 1.4039647 1.118778e-05       Thal-LS    295    1383      all   pairwise    0.05
# 44 0.3393501 1.000000e+00       Thal-LS     44     605 positive   pairwise    0.05
# 45 2.6926299 5.437591e-28       Thal-LS    251     778 negative   pairwise    0.05
# 46 1.0115469 4.778049e-01       TNoS-LS    120    1383      all   pairwise    0.05
# 47 0.5430901 9.997892e-01       TNoS-LS     31     605 positive   pairwise    0.05
# 48 1.4582757 1.994894e-03       TNoS-LS     89     778 negative   pairwise    0.05
# 49 1.4988115 1.964167e-07   TT.IG.SH-LS    303    1383      all   pairwise    0.05
# 50 0.2656862 1.000000e+00   TT.IG.SH-LS     35     605 positive   pairwise    0.05
# 51 3.1053385 1.219563e-36   TT.IG.SH-LS    268     778 negative   pairwise    0.05

prwiseTab_ct_LS_broad_glFDR01_ctFDR05 <- spatialLIBD::gene_set_enrichment(gene_list = gene_list_FDR01, modeling_results = modeling_results_broad, fdr_cut = 0.05, model_type = "pairwise", reverse = FALSE)
prwiseTab_ct_LS_broad_glFDR01_ctFDR05
#            OR         Pval          test NumSig SetSize       ID model_type fdr_cut
# 1   1.0591726 4.324529e-01      Astro-LS     23     258      all   pairwise    0.05
# 2   1.3407440 1.341664e-01      Astro-LS     21     192 positive   pairwise    0.05
# 3   0.3343983 9.800004e-01      Astro-LS      2      66 negative   pairwise    0.05
# 4   0.3863778 9.995568e-01       Chol-LS      9     258      all   pairwise    0.05
# 5   0.2278091 9.999515e-01       Chol-LS      4     192 positive   pairwise    0.05
# 6   0.9027500 6.521156e-01       Chol-LS      5      66 negative   pairwise    0.05
# 7   0.8896118 7.989676e-01        ChP-LS     58     258      all   pairwise    0.05
# 8   1.2155715 1.343096e-01        ChP-LS     54     192 positive   pairwise    0.05
# 9   0.1965960 9.999852e-01        ChP-LS      4      66 negative   pairwise    0.05
# 10  2.3261115 1.002929e-08       Endo-LS     78     258      all   pairwise    0.05
# 11  3.3970519 4.746356e-14       Endo-LS     74     192 positive   pairwise    0.05
# 12  0.3261162 9.968726e-01       Endo-LS      4      66 negative   pairwise    0.05
# 13  0.6913568 9.926936e-01  Ependymal-LS     50     258      all   pairwise    0.05
# 14  0.9179458 7.168672e-01  Ependymal-LS     46     192 positive   pairwise    0.05
# 15  0.1862045 9.999930e-01  Ependymal-LS      4      66 negative   pairwise    0.05
# 16  1.0333397 4.755168e-01        IoC-LS     23     258      all   pairwise    0.05
# 17  0.5701177 9.756354e-01        IoC-LS     10     192 positive   pairwise    0.05
# 18  2.6299429 3.744903e-03        IoC-LS     13      66 negative   pairwise    0.05
# 19  8.6547323 5.880038e-45      Micro-LS    101     258      all   pairwise    0.05
# 20 14.7695725 4.951622e-59      Micro-LS    100     192 positive   pairwise    0.05
# 21  0.1650259 9.971458e-01      Micro-LS      1      66 negative   pairwise    0.05
# 22  0.5475514 9.977842e-01         MS-LS     20     258      all   pairwise    0.05
# 23  0.2084043 9.999997e-01         MS-LS      6     192 positive   pairwise    0.05
# 24  1.8107675 4.241741e-02         MS-LS     14      66 negative   pairwise    0.05
# 25  2.0471840 2.110131e-03      Mural-LS     24     258      all   pairwise    0.05
# 26  2.5952609 1.931483e-04      Mural-LS     22     192 positive   pairwise    0.05
# 27  0.5926215 8.489833e-01      Mural-LS      2      66 negative   pairwise    0.05
# 28  0.7340579 9.214306e-01 Neuroblast-LS     19     258      all   pairwise    0.05
# 29  0.7866806 8.431622e-01 Neuroblast-LS     15     192 positive   pairwise    0.05
# 30  0.6005045 8.929055e-01 Neuroblast-LS      4      66 negative   pairwise    0.05
# 31  1.5552777 3.023449e-02      Oligo-LS     26     258      all   pairwise    0.05
# 32  2.0992844 1.349224e-03      Oligo-LS     25     192 positive   pairwise    0.05
# 33  0.2061365 9.912088e-01      Oligo-LS      1      66 negative   pairwise    0.05
# 34  0.9077740 7.063797e-01        OPC-LS     26     258      all   pairwise    0.05
# 35  0.9441894 6.297706e-01        OPC-LS     20     192 positive   pairwise    0.05
# 36  0.8117649 7.439865e-01        OPC-LS      6      66 negative   pairwise    0.05
# 37  0.8745140 7.383614e-01       Sept-LS     18     258      all   pairwise    0.05
# 38  0.4346426 9.949135e-01       Sept-LS      7     192 positive   pairwise    0.05
# 39  2.3811902 1.282305e-02       Sept-LS     11      66 negative   pairwise    0.05
# 40  0.7439150 9.518216e-01        Str-LS     32     258      all   pairwise    0.05
# 41  0.1372808 1.000000e+00        Str-LS      5     192 positive   pairwise    0.05
# 42  3.7694079 7.926208e-07        Str-LS     27      66 negative   pairwise    0.05
# 43  0.9132051 7.257208e-01       Thal-LS     42     258      all   pairwise    0.05
# 44  0.1983375 1.000000e+00       Thal-LS      8     192 positive   pairwise    0.05
# 45  5.1573936 2.429621e-10       Thal-LS     34      66 negative   pairwise    0.05
# 46  0.6439004 9.661009e-01       TNoS-LS     15     258      all   pairwise    0.05
# 47  0.5127461 9.879961e-01       TNoS-LS      9     192 positive   pairwise    0.05
# 48  1.0623479 5.078197e-01       TNoS-LS      6      66 negative   pairwise    0.05
# 49  0.4986322 9.998606e-01   TT.IG.SH-LS     25     258      all   pairwise    0.05
# 50  0.0727878 1.000000e+00   TT.IG.SH-LS      3     192 positive   pairwise    0.05
# 51  2.4179742 1.144005e-03   TT.IG.SH-LS     22      66 negative   pairwise    0.05

###############################################################################



################ gene_set_enrichment analysis with 1vsAll data ################

## Running gene_set_enrichment with two sets of DE genes, one with FDR < 0.05 the other with FDR < 0.01
enrichTab_glFDR05_ctFDR05 <- gene_set_enrichment(gene_list = gene_list_FDR05, modeling_results = modeling_results, model_type = "enrichment", fdr_cut = 0.05)
enrichTab_glFDR05_ctFDR05
#           OR         Pval          test NumSig SetSize       ID model_type fdr_cut
# 1  1.6759137 1.528135e-08     Chol_Ex.D    211    1099      all enrichment    0.05
# 2  0.6724000 9.783978e-01     Chol_Ex.D     25     252 positive enrichment    0.05
# 3  2.0117813 1.024594e-12     Chol_Ex.D    186     847 negative enrichment    0.05
# 4  1.8448670 4.420322e-16       LS_In.C    369    1099      all enrichment    0.05
# 5  0.4603787 9.999971e-01       LS_In.C     33     252 positive enrichment    0.05
# 6  2.4718043 8.949557e-29       LS_In.C    336     847 negative enrichment    0.05
# 7  1.7342849 3.089941e-15       LS_In.D    479    1099      all enrichment    0.05
# 8  0.1880946 1.000000e+00       LS_In.D     23     252 positive enrichment    0.05
# 9  2.7812018 1.088146e-40       LS_In.D    456     847 negative enrichment    0.05
# 10 1.4459703 2.171798e-06       LS_In.M    296    1099      all enrichment    0.05
# 11 0.3313156 1.000000e+00       LS_In.M     22     252 positive enrichment    0.05
# 12 1.9574572 1.980449e-15       LS_In.M    274     847 negative enrichment    0.05
# 13 1.9503481 3.721572e-10       LS_In.N    157    1099      all enrichment    0.05
# 14 0.3927903 9.996721e-01       LS_In.N     10     252 positive enrichment    0.05
# 15 2.5373747 1.234385e-16       LS_In.N    147     847 negative enrichment    0.05
# 16 2.4064508 1.006164e-20       LS_In.O    228    1099      all enrichment    0.05
# 17 0.3489028 9.999877e-01       LS_In.O     12     252 positive enrichment    0.05
# 18 3.2556307 7.603599e-33       LS_In.O    216     847 negative enrichment    0.05
# 19 1.5897604 4.748146e-11       LS_In.P    440    1099      all enrichment    0.05
# 20 0.2343743 1.000000e+00       LS_In.P     26     252 positive enrichment    0.05
# 21 2.4051696 5.063967e-30       LS_In.P    414     847 negative enrichment    0.05
# 22 1.9611287 4.140115e-12       LS_In.Q    195    1099      all enrichment    0.05
# 23 0.5039550 9.985001e-01       LS_In.Q     16     252 positive enrichment    0.05
# 24 2.4944015 7.015233e-19       LS_In.Q    179     847 negative enrichment    0.05
# 25 1.4896438 1.006550e-07       LS_In.R    341    1099      all enrichment    0.05
# 26 0.4093834 9.999999e-01       LS_In.R     31     252 positive enrichment    0.05
# 27 1.9822364 4.446531e-17       LS_In.R    310     847 negative enrichment    0.05
# 28 1.5918534 7.113863e-10       MS_In.J    345    1099      all enrichment    0.05
# 29 0.3761087 1.000000e+00       MS_In.J     28     252 positive enrichment    0.05
# 30 2.1623963 3.521923e-21       MS_In.J    317     847 negative enrichment    0.05
# 31 1.6411978 2.677483e-12       MS_In.K    437    1099      all enrichment    0.05
# 32 0.2220070 1.000000e+00       MS_In.K     24     252 positive enrichment    0.05
# 33 2.5016237 1.955792e-32       MS_In.K    413     847 negative enrichment    0.05
# 34 1.8270050 5.597970e-12     Sept_In.G    244    1099      all enrichment    0.05
# 35 0.7338369 9.546976e-01     Sept_In.G     30     252 positive enrichment    0.05
# 36 2.1817661 9.598863e-17     Sept_In.G    214     847 negative enrichment    0.05
# 37 1.7552004 2.334696e-15     Sept_In.I    446    1099      all enrichment    0.05
# 38 0.2047548 1.000000e+00     Sept_In.I     22     252 positive enrichment    0.05
# 39 2.7265192 1.953125e-38     Sept_In.I    424     847 negative enrichment    0.05
# 40 2.0392981 3.884863e-17     TNoS_Ex.A    277    1099      all enrichment    0.05
# 41 0.3519989 9.999993e-01     TNoS_Ex.A     17     252 positive enrichment    0.05
# 42 2.7796804 5.809717e-30     TNoS_Ex.A    260     847 negative enrichment    0.05
# 43 2.1561354 1.589634e-25 TT.IG.SH_Ex.C    422    1099      all enrichment    0.05
# 44 0.4481697 9.999992e-01 TT.IG.SH_Ex.C     35     252 positive enrichment    0.05
# 45 2.9957601 2.521425e-43 TT.IG.SH_Ex.C    387     847 negative enrichment    0.05
# 46 1.8572679 2.131629e-10 TT.IG.SH_Ex.E    189    1099      all enrichment    0.05
# 47 0.3359908 9.999884e-01 TT.IG.SH_Ex.E     11     252 positive enrichment    0.05
# 48 2.4701799 1.754582e-18 TT.IG.SH_Ex.E    178     847 negative enrichment    0.05
# 49 2.2286815 1.344878e-30 TT.IG.SH_Ex.F    533    1099      all enrichment    0.05
# 50 0.2525851 1.000000e+00 TT.IG.SH_Ex.F     30     252 positive enrichment    0.05
# 51 3.6401316 1.595154e-63 TT.IG.SH_Ex.F    503     847 negative enrichment    0.05

enrichTab_glFDR01_ctFDR05 <- gene_set_enrichment(gene_list = gene_list_FDR01, modeling_results = modeling_results, model_type = "enrichment", fdr_cut = 0.05)
enrichTab_glFDR01_ctFDR05
#            OR         Pval          test NumSig SetSize       ID model_type fdr_cut
# 1  1.60632716 3.487201e-02     Chol_Ex.D     23     113      all enrichment    0.05
# 2  0.73742687 8.003561e-01     Chol_Ex.D      5      47 positive enrichment    0.05
# 3  2.36256260 2.940469e-03     Chol_Ex.D     18      66 negative enrichment    0.05
# 4  1.75045571 4.199874e-03       LS_In.C     40     113      all enrichment    0.05
# 5  0.85033138 7.277400e-01       LS_In.C     10      47 positive enrichment    0.05
# 6  2.66535808 1.083248e-04       LS_In.C     30      66 negative enrichment    0.05
# 7  1.36940818 6.478599e-02       LS_In.D     46     113      all enrichment    0.05
# 8  0.08696972 9.999999e-01       LS_In.D      2      47 positive enrichment    0.05
# 9  4.03696867 3.291828e-08       LS_In.D     44      66 negative enrichment    0.05
# 10 1.37177438 8.768624e-02       LS_In.M     31     113      all enrichment    0.05
# 11 0.24355565 9.990665e-01       LS_In.M      3      47 positive enrichment    0.05
# 12 2.69561368 1.178266e-04       LS_In.M     28      66 negative enrichment    0.05
# 13 1.89327816 1.485447e-02       LS_In.N     18     113      all enrichment    0.05
# 14 0.21172761 9.897652e-01       LS_In.N      1      47 positive enrichment    0.05
# 15 3.49307044 6.870110e-05       LS_In.N     17      66 negative enrichment    0.05
# 16 2.95505522 2.290595e-06       LS_In.O     32     113      all enrichment    0.05
# 17 0.15580260 9.978009e-01       LS_In.O      1      47 positive enrichment    0.05
# 18 6.68031082 2.116354e-12       LS_In.O     31      66 negative enrichment    0.05
# 19 1.27457189 1.292836e-01       LS_In.P     42     113      all enrichment    0.05
# 20 0.31111189 9.993550e-01       LS_In.P      6      47 positive enrichment    0.05
# 21 2.60636329 1.004923e-04       LS_In.P     36      66 negative enrichment    0.05
# 22 1.87662813 9.357558e-03       LS_In.Q     22     113      all enrichment    0.05
# 23 0.16448477 9.970218e-01       LS_In.Q      1      47 positive enrichment    0.05
# 24 3.65289801 8.851459e-06       LS_In.Q     21      66 negative enrichment    0.05
# 25 1.36524364 8.194583e-02       LS_In.R     35     113      all enrichment    0.05
# 26 0.20409905 9.998131e-01       LS_In.R      3      47 positive enrichment    0.05
# 27 2.88937172 2.548509e-05       LS_In.R     32      66 negative enrichment    0.05
# 28 1.02859487 4.859274e-01       MS_In.J     28     113      all enrichment    0.05
# 29 0.06710463 9.999980e-01       MS_In.J      1      47 positive enrichment    0.05
# 30 2.18578845 1.959949e-03       MS_In.J     27      66 negative enrichment    0.05
# 31 1.32160109 9.528010e-02       MS_In.K     42     113      all enrichment    0.05
# 32 0.26202841 9.997790e-01       MS_In.K      5      47 positive enrichment    0.05
# 33 2.87515319 1.993241e-05       MS_In.K     37      66 negative enrichment    0.05
# 34 1.84310375 5.680203e-03     Sept_In.G     28     113      all enrichment    0.05
# 35 0.50966625 9.450010e-01     Sept_In.G      4      47 positive enrichment    0.05
# 36 3.21213731 2.112987e-05     Sept_In.G     24      66 negative enrichment    0.05
# 37 1.34406549 8.193471e-02     Sept_In.I     42     113      all enrichment    0.05
# 38 0.09917487 9.999993e-01     Sept_In.I      2      47 positive enrichment    0.05
# 39 3.53526251 4.244952e-07     Sept_In.I     40      66 negative enrichment    0.05
# 40 1.68421432 1.521883e-02     TNoS_Ex.A     28     113      all enrichment    0.05
# 41 0.34168317 9.895069e-01     TNoS_Ex.A      3      47 positive enrichment    0.05
# 42 3.13820727 2.312554e-05     TNoS_Ex.A     25      66 negative enrichment    0.05
# 43 1.65360326 8.541093e-03 TT.IG.SH_Ex.C     41     113      all enrichment    0.05
# 44 0.49941032 9.763931e-01 TT.IG.SH_Ex.C      7      47 positive enrichment    0.05
# 45 3.10101365 6.467075e-06 TT.IG.SH_Ex.C     34      66 negative enrichment    0.05
# 46 2.10156270 2.106565e-03 TT.IG.SH_Ex.E     24     113      all enrichment    0.05
# 47 0.51781889 9.221497e-01 TT.IG.SH_Ex.E      3      47 positive enrichment    0.05
# 48 3.65289801 8.851459e-06 TT.IG.SH_Ex.E     21      66 negative enrichment    0.05
# 49 1.16622002 2.462071e-01 TT.IG.SH_Ex.F     42     113      all enrichment    0.05
# 50 0.13247811 9.999989e-01 TT.IG.SH_Ex.F      3      47 positive enrichment    0.05
# 51 2.87969675 1.908361e-05 TT.IG.SH_Ex.F     39      66 negative enrichment    0.05

###############################################################################



################# gene_set_enrichment analysis with 1vs1 data #################

prwiseTab_glFDR05_ctFDR05 <- spatialLIBD::gene_set_enrichment(gene_list = gene_list_FDR05, modeling_results = modeling_results, fdr_cut = 0.05, model_type = "pairwise", reverse = FALSE)

prwiseTab_glFDR01_ctFDR05 <- spatialLIBD::gene_set_enrichment(gene_list = gene_list_FDR01, modeling_results = modeling_results, fdr_cut = 0.05, model_type = "pairwise", reverse = FALSE)

###############################################################################



################### Save gene_set_enrichment results to rda ###################

save(enrichTab_broad_glFDR01_ctFDR05,
    enrichTab_broad_glFDR05_ctFDR05,
    enrichTab_glFDR01_ctFDR05,
    enrichTab_glFDR05_ctFDR05,
    file = here(
        "snRNAseq_mouse",
        "processed_data",
        "SCE",
        "gene_set_enrichment_1vsAll_result_tables.rda"
    )
)

save(prwiseTab_LS_ct_broad_glFDR01_ctFDR05,
    prwiseTab_LS_ct_broad_glFDR05_ctFDR05,
    prwiseTab_ct_LS_broad_glFDR01_ctFDR05,
    prwiseTab_ct_LS_broad_glFDR05_ctFDR05,
    prwiseTab_glFDR01_ctFDR05,
    prwiseTab_glFDR05_ctFDR05,
    file = here(
        "snRNAseq_mouse",
        "processed_data",
        "SCE",
        "gene_set_enrichment_1vs1_result_tables.rda"
    )
)

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
