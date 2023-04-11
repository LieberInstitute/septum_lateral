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
# 4.31 MB
lobstr::obj_size(modeling_results)
# 11.88 MB
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
dim(modeling_results_broad$pairwise)
# [1] 2635   53
head(names(modeling_results_broad$pairwise), 10)
# [1] "t_stat_LS-Astro"  "p_value_LS-Astro" "fdr_LS-Astro"     "t_stat_LS-Chol"
# [5] "p_value_LS-Chol"  "fdr_LS-Chol"      "t_stat_LS-ChP"    "p_value_LS-ChP"
# [9] "fdr_LS-ChP"       "t_stat_LS-Endo"

class(modeling_results)
# [1] "list"
names(modeling_results)
# [1] "enrichment" "pairwise"
dim(modeling_results$enrichment)
# [1] 4609   32
head(names(modeling_results$enrichment), 10)
# [1] "t_stat_LS_In.C"  "p_value_LS_In.C" "fdr_LS_In.C"     "t_stat_LS_In.D"
# [5] "p_value_LS_In.D" "fdr_LS_In.D"     "t_stat_LS_In.M"  "p_value_LS_In.M"
# [9] "fdr_LS_In.M"     "t_stat_LS_In.N
dim(modeling_results$pairwise)
# [1] 4609  272
head(names(modeling_results$pairwise), 10)
# [1] "t_stat_LS_In.C-LS_In.D"  "p_value_LS_In.C-LS_In.D"
# [3] "fdr_LS_In.C-LS_In.D"     "t_stat_LS_In.C-LS_In.M"
# [5] "p_value_LS_In.C-LS_In.M" "fdr_LS_In.C-LS_In.M"
# [7] "t_stat_LS_In.C-LS_In.N"  "p_value_LS_In.C-LS_In.N"
# [9] "fdr_LS_In.C-LS_In.N"     "t_stat_LS_In.C-LS_In.O"

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

prwiseTab_broad_glFDR05_ctFDR05 <- spatialLIBD::gene_set_enrichment(gene_list = gene_list_FDR05, modeling_results = modeling_results_broad, fdr_cut = 0.05, model_type = "pairwise", reverse = FALSE)
prwiseTab_broad_glFDR05_ctFDR05
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

prwiseTab_broad_glFDR01_ctFDR05 <- spatialLIBD::gene_set_enrichment(gene_list = gene_list_FDR01, modeling_results = modeling_results_broad, fdr_cut = 0.05, model_type = "pairwise", reverse = FALSE)
# Warning message:
# In spatialLIBD::gene_set_enrichment(gene_list = gene_list_FDR01,  :
#   Gene list with n < 25 may have insufficent power for enrichment analysis: all ,positive ,negative
prwiseTab_broad_glFDR01_ctFDR05
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

###############################################################################



################ gene_set_enrichment analysis with 1vsAll data ################

## Running gene_set_enrichment with two sets of DE genes, one with FDR < 0.05 the other with FDR < 0.01
enrichTab_glFDR05_ctFDR05 <- gene_set_enrichment(gene_list = gene_list_FDR05, modeling_results = modeling_results, model_type = "enrichment", fdr_cut = 0.05)
enrichTab_glFDR05_ctFDR05
#           OR         Pval      test NumSig SetSize       ID model_type fdr_cut
# 1  1.8235459 4.825621e-15   LS_In.C    369    1008      all enrichment    0.05
# 2  0.4900993 9.999758e-01   LS_In.C     33     212 positive enrichment    0.05
# 3  2.3671167 1.531071e-25   LS_In.C    336     796 negative enrichment    0.05
# 4  1.7217953 4.121137e-14   LS_In.D    479    1008      all enrichment    0.05
# 5  0.1936702 1.000000e+00   LS_In.D     23     212 positive enrichment    0.05
# 6  2.7040083 1.638271e-36   LS_In.D    456     796 negative enrichment    0.05
# 7  1.4142860 1.049409e-05   LS_In.M    296    1008      all enrichment    0.05
# 8  0.3504751 9.999999e-01   LS_In.M     22     212 positive enrichment    0.05
# 9  1.8575712 3.976873e-13   LS_In.M    274     796 negative enrichment    0.05
# 10 1.9110004 1.433212e-09   LS_In.N    157    1008      all enrichment    0.05
# 11 0.4196753 9.991092e-01   LS_In.N     10     212 positive enrichment    0.05
# 12 2.4140415 4.916727e-15   LS_In.N    147     796 negative enrichment    0.05
# 13 2.3719793 7.922417e-20   LS_In.O    228    1008      all enrichment    0.05
# 14 0.3718389 9.999522e-01   LS_In.O     12     212 positive enrichment    0.05
# 15 3.1155807 2.722198e-30   LS_In.O    216     796 negative enrichment    0.05
# 16 1.5674136 5.018405e-10   LS_In.P    440    1008      all enrichment    0.05
# 17 0.2432223 1.000000e+00   LS_In.P     26     212 positive enrichment    0.05
# 18 2.3113620 2.639740e-26   LS_In.P    414     796 negative enrichment    0.05
# 19 1.9245763 2.086235e-11   LS_In.Q    195    1008      all enrichment    0.05
# 20 0.5394298 9.958132e-01   LS_In.Q     16     212 positive enrichment    0.05
# 21 2.3748882 5.528764e-17   LS_In.Q    179     796 negative enrichment    0.05
# 22 1.4597263 6.356542e-07   LS_In.R    341    1008      all enrichment    0.05
# 23 0.4341617 9.999987e-01   LS_In.R     31     212 positive enrichment    0.05
# 24 1.8822869 1.737899e-14   LS_In.R    310     796 negative enrichment    0.05
# 25 1.7945339 3.422905e-11 Sept_In.G    244    1008      all enrichment    0.05
# 26 0.7913785 8.984085e-01 Sept_In.G     30     212 positive enrichment    0.05
# 27 2.0744899 1.028981e-14 Sept_In.G    214     796 negative enrichment    0.05
# 28 1.7397001 2.955528e-14 Sept_In.I    446    1008      all enrichment    0.05
# 29 0.2123094 1.000000e+00 Sept_In.I     22     212 positive enrichment    0.05
# 30 2.6387371 1.873380e-34 Sept_In.I    424     796 negative enrichment    0.05

enrichTab_glFDR01_ctFDR05 <- gene_set_enrichment(gene_list = gene_list_FDR01, modeling_results = modeling_results, model_type = "enrichment", fdr_cut = 0.05)
enrichTab_glFDR01_ctFDR05
#            OR         Pval      test NumSig SetSize       ID model_type fdr_cut
# 1  1.75890710 4.712605e-03   LS_In.C     40     103      all enrichment    0.05
# 2  0.85258311 7.244407e-01   LS_In.C     10      42 positive enrichment    0.05
# 3  2.68410669 1.427746e-04   LS_In.C     30      61 negative enrichment    0.05
# 4  1.36514463 7.389330e-02   LS_In.D     46     103      all enrichment    0.05
# 5  0.08293878 9.999999e-01   LS_In.D      2      42 positive enrichment    0.05
# 6  4.43370485 2.904790e-08   LS_In.D     44      61 negative enrichment    0.05
# 7  1.36074373 9.795778e-02   LS_In.M     31     103      all enrichment    0.05
# 8  0.23932901 9.991204e-01   LS_In.M      3      42 positive enrichment    0.05
# 9  2.70413107 1.556728e-04   LS_In.M     28      61 negative enrichment    0.05
# 10 1.88037346 1.638193e-02   LS_In.N     18     103      all enrichment    0.05
# 11 0.21114669 9.897404e-01   LS_In.N      1      42 positive enrichment    0.05
# 12 3.45721379 8.803228e-05   LS_In.N     17      61 negative enrichment    0.05
# 13 2.98463183 2.495454e-06   LS_In.O     32     103      all enrichment    0.05
# 14 0.15472600 9.978207e-01   LS_In.O      1      42 positive enrichment    0.05
# 15 6.89968316 2.472182e-12   LS_In.O     31      61 negative enrichment    0.05
# 16 1.26386378 1.464266e-01   LS_In.P     42     103      all enrichment    0.05
# 17 0.30180544 9.994415e-01   LS_In.P      6      42 positive enrichment    0.05
# 18 2.66543211 1.288062e-04   LS_In.P     36      61 negative enrichment    0.05
# 19 1.86761037 1.042370e-02   LS_In.Q     22     103      all enrichment    0.05
# 20 0.16348513 9.970404e-01   LS_In.Q      1      42 positive enrichment    0.05
# 21 3.64086883 1.159161e-05   LS_In.Q     21      61 negative enrichment    0.05
# 22 1.35557369 9.217997e-02   LS_In.R     35     103      all enrichment    0.05
# 23 0.19930023 9.998314e-01   LS_In.R      3      42 positive enrichment    0.05
# 24 2.93314489 3.304163e-05   LS_In.R     32      61 negative enrichment    0.05
# 25 1.83981639 6.382043e-03 Sept_In.G     28     103      all enrichment    0.05
# 26 0.50798536 9.449126e-01 Sept_In.G      4      42 positive enrichment    0.05
# 27 3.21156524 2.781431e-05 Sept_In.G     24      61 negative enrichment    0.05
# 28 1.33646442 9.304785e-02 Sept_In.I     42     103      all enrichment    0.05
# 29 0.09528581 9.999995e-01 Sept_In.I      2      42 positive enrichment    0.05
# 30 3.74110671 4.619979e-07 Sept_In.I     40      61 negative enrichment    0.05

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

save(prwiseTab_broad_glFDR01_ctFDR05,
    prwiseTab_broad_glFDR05_ctFDR05,
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
