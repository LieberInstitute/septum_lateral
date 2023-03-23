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
#           OR         Pval       test NumSig SetSize       ID model_type fdr_cut
# 1  0.7201371 9.955595e-01      Astro    122     254      all enrichment     0.1
# 2  0.4047635 1.000000e+00      Astro     66     191 positive enrichment     0.1
# 3  6.4335193 1.196493e-08      Astro     56      63 negative enrichment     0.1
# 4  0.7323689 9.874757e-01       Chol     66     254      all enrichment     0.1
# 5  0.1857971 1.000000e+00       Chol     16     191 positive enrichment     0.1
# 6  8.3487628 1.114420e-14       Chol     50      63 negative enrichment     0.1
# 7  0.6600557 9.994517e-01        ChP    102     254      all enrichment     0.1
# 8  0.3338661 1.000000e+00        ChP     49     191 positive enrichment     0.1
# 9  5.4023660 1.330319e-08        ChP     53      63 negative enrichment     0.1
# 10 0.6654264 9.993085e-01       Endo    101     254      all enrichment     0.1
# 11 0.3327230 1.000000e+00       Endo     48     191 positive enrichment     0.1
# 12 5.5357139 7.790609e-09       Endo     53      63 negative enrichment     0.1
# 13 0.7038290 9.973599e-01  Ependymal    128     254      all enrichment     0.1
# 14 0.4129035 1.000000e+00  Ependymal     72     191 positive enrichment     0.1
# 15 5.7214936 1.270046e-07  Ependymal     56      63 negative enrichment     0.1
# 16 0.7389478 9.841231e-01        IoC     63     254      all enrichment     0.1
# 17 0.1725121 1.000000e+00        IoC     14     191 positive enrichment     0.1
# 18 8.1608066 1.051887e-14        IoC     49      63 negative enrichment     0.1
# 19 0.7400886 9.912956e-01         LS    141     254      all enrichment     0.1
# 20 0.4691333 9.999999e-01         LS     85     191 positive enrichment     0.1
# 21 4.8798841 2.257994e-06         LS     56      63 negative enrichment     0.1
# 22 0.5973143 9.999574e-01      Micro     84     254      all enrichment     0.1
# 23 0.2388725 1.000000e+00      Micro     32     191 positive enrichment     0.1
# 24 5.9564898 5.700690e-10      Micro     52      63 negative enrichment     0.1
# 25 0.7177497 9.958342e-01         MS    128     254      all enrichment     0.1
# 26 0.4209919 1.000000e+00         MS     72     191 positive enrichment     0.1
# 27 5.8299309 8.823116e-08         MS     56      63 negative enrichment     0.1
# 28 0.5812538 9.999793e-01      Mural     80     254      all enrichment     0.1
# 29 0.2130867 1.000000e+00      Mural     28     191 positive enrichment     0.1
# 30 6.2396232 1.838616e-10      Mural     52      63 negative enrichment     0.1
# 31 0.7038607 9.973630e-01 Neuroblast    113     254      all enrichment     0.1
# 32 0.3766459 1.000000e+00 Neuroblast     58     191 positive enrichment     0.1
# 33 6.2410022 5.368606e-09 Neuroblast     55      63 negative enrichment     0.1
# 34 0.6102602 9.998763e-01      Oligo     74     254      all enrichment     0.1
# 35 0.1891390 1.000000e+00      Oligo     22     191 positive enrichment     0.1
# 36 7.3286987 2.880502e-12      Oligo     52      63 negative enrichment     0.1
# 37 0.6163673 9.998419e-01        OPC     75     254      all enrichment     0.1
# 38 0.1971528 1.000000e+00        OPC     23     191 positive enrichment     0.1
# 39 7.2588604 3.730497e-12        OPC     52      63 negative enrichment     0.1
# 40 0.7345770 9.926059e-01       Sept    140     254      all enrichment     0.1
# 41 0.4631464 9.999999e-01       Sept     84     191 positive enrichment     0.1
# 42 4.9234972 1.941246e-06       Sept     56      63 negative enrichment     0.1
# 43 0.7278418 9.943339e-01        Str    126     254      all enrichment     0.1
# 44 0.4210690 1.000000e+00        Str     70     191 positive enrichment     0.1
# 45 6.0992009 3.595827e-08        Str     56      63 negative enrichment     0.1
# 46 0.7113415 9.966108e-01       Thal    128     254      all enrichment     0.1
# 47 0.4172688 1.000000e+00       Thal     72     191 positive enrichment     0.1
# 48 5.7800185 1.043165e-07       Thal     56      63 negative enrichment     0.1
# 49 0.6076834 9.999120e-01       TNoS     79     254      all enrichment     0.1
# 50 0.2172528 1.000000e+00       TNoS     27     191 positive enrichment     0.1
# 51 6.6360946 3.898949e-11       TNoS     52      63 negative enrichment     0.1
# 52 0.7103653 9.967607e-01   TT.IG.SH    122     254      all enrichment     0.1
# 53 0.3993264 1.000000e+00   TT.IG.SH     66     191 positive enrichment     0.1
# 54 6.3496877 1.574342e-08   TT.IG.SH     56      63 negative enrichment     0.1

enrichTab_FDR05
#           OR         Pval       test NumSig SetSize       ID model_type fdr_cut
# 1  1.1839041 4.342524e-03      Astro    801    1360      all enrichment     0.1
# 2  0.3216382 1.000000e+00      Astro    190     602 positive enrichment     0.1
# 3  3.8836483 9.584513e-54      Astro    611     758 negative enrichment     0.1
# 4  1.6941314 7.399120e-16       Chol    557    1360      all enrichment     0.1
# 5  0.2349213 1.000000e+00       Chol     67     602 positive enrichment     0.1
# 6  5.0234218 6.904242e-89       Chol    490     758 negative enrichment     0.1
# 7  1.2189418 9.205241e-04        ChP    729    1360      all enrichment     0.1
# 8  0.2816549 1.000000e+00        ChP    146     602 positive enrichment     0.1
# 9  4.0002444 1.339747e-60        ChP    583     758 negative enrichment     0.1
# 10 1.2343742 4.575959e-04       Endo    724    1360      all enrichment     0.1
# 11 0.2837443 1.000000e+00       Endo    144     602 positive enrichment     0.1
# 12 4.0145038 2.617540e-61       Endo    580     758 negative enrichment     0.1
# 13 1.1088073 5.703782e-02  Ependymal    823    1360      all enrichment     0.1
# 14 0.3218454 1.000000e+00  Ependymal    206     602 positive enrichment     0.1
# 15 3.5978296 1.154929e-46  Ependymal    617     758 negative enrichment     0.1
# 16 1.6826648 3.375873e-15        IoC    532    1360      all enrichment     0.1
# 17 0.2227243 1.000000e+00        IoC     60     602 positive enrichment     0.1
# 18 4.8719544 1.304711e-85        IoC    472     758 negative enrichment     0.1
# 19 1.0701686 1.567955e-01         LS    865    1360      all enrichment     0.1
# 20 0.3202850 1.000000e+00         LS    228     602 positive enrichment     0.1
# 21 3.6835028 1.283170e-44         LS    637     758 negative enrichment     0.1
# 22 1.2809362 4.897104e-05      Micro    670    1360      all enrichment     0.1
# 23 0.2444891 1.000000e+00      Micro    111     602 positive enrichment     0.1
# 24 4.2379306 1.098409e-68      Micro    559     758 negative enrichment     0.1
# 25 1.1228233 3.764875e-02         MS    820    1360      all enrichment     0.1
# 26 0.3261351 1.000000e+00         MS    205     602 positive enrichment     0.1
# 27 3.6055438 3.634314e-47         MS    615     758 negative enrichment     0.1
# 28 1.3098820 1.101963e-05      Mural    660    1360      all enrichment     0.1
# 29 0.2392201 1.000000e+00      Mural    105     602 positive enrichment     0.1
# 30 4.3429908 2.104907e-71      Mural    555     758 negative enrichment     0.1
# 31 1.2134485 1.213996e-03 Neuroblast    767    1360      all enrichment     0.1
# 32 0.3041092 1.000000e+00 Neuroblast    168     602 positive enrichment     0.1
# 33 4.0052977 2.301670e-58 Neuroblast    599     758 negative enrichment     0.1
# 34 1.3842783 1.964804e-07      Oligo    620    1360      all enrichment     0.1
# 35 0.2276607 1.000000e+00      Oligo     88     602 positive enrichment     0.1
# 36 4.4381449 2.715118e-75      Oligo    532     758 negative enrichment     0.1
# 37 1.3776084 2.869702e-07        OPC    622    1360      all enrichment     0.1
# 38 0.2285101 1.000000e+00        OPC     89     602 positive enrichment     0.1
# 39 4.4189447 8.272987e-75        OPC    533     758 negative enrichment     0.1
# 40 1.0782868 1.303533e-01       Sept    864    1360      all enrichment     0.1
# 41 0.3210299 1.000000e+00       Sept    227     602 positive enrichment     0.1
# 42 3.7203735 2.419678e-45       Sept    637     758 negative enrichment     0.1
# 43 1.1583487 1.152494e-02        Str    813    1360      all enrichment     0.1
# 44 0.3266551 1.000000e+00        Str    199     602 positive enrichment     0.1
# 45 3.7591492 1.379496e-50        Str    614     758 negative enrichment     0.1
# 46 1.1146684 4.808382e-02       Thal    821    1360      all enrichment     0.1
# 47 0.3256236 1.000000e+00       Thal    206     602 positive enrichment     0.1
# 48 3.5706133 1.996935e-46       Thal    615     758 negative enrichment     0.1
# 49 1.3453694 1.631918e-06       TNoS    646    1360      all enrichment     0.1
# 50 0.2341909 1.000000e+00       TNoS     98     602 positive enrichment     0.1
# 51 4.4335286 5.025350e-74       TNoS    548     758 negative enrichment     0.1
# 52 1.1777769 5.535827e-03   TT.IG.SH    804    1360      all enrichment     0.1
# 53 0.3250915 1.000000e+00   TT.IG.SH    193     602 positive enrichment     0.1
# 54 3.8260189 1.594080e-52   TT.IG.SH    611     758 negative enrichment     0.1

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


################################################################
#### Plot enrichment using gene_set_enrichment_plot_complex ####
################################################################

source(
    here(
        "snRNAseq_mouse",
        "code",
        "04_clinical_set_enrichment",
        "gene_set_enrichment_plot_complex.R"
        )
    )


use_gsepc<-function(modeling_results, model_type, gene_list, enrichTab, path_to_plot){
    gene_enrichment_count <- get_gene_enrichment_count(model_results = modeling_results, model_type = model_type, bayes_anno = NULL)
    gene_list_count <- get_gene_list_count(gene_list)

    gse_plot<-gene_set_enrichment_plot_complex(
        enrichment = enrichTab,
        gene_count_col = log10(gene_list_count),
        gene_count_row = gene_enrichment_count,
        anno_title_col = "DE Genes\n(log10)",
        anno_title_row = "Cluster\nGenes"
        )

    pdf(path_to_plot, height = 4, width = 6)
    print(gse_plot)
    dev.off()
}

use_gsepc(
    modeling_results = modeling_results,
    model_type = "enrichment",
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
    model_type = "enrichment",
    gene_list = gene_list_FDR05,
    enrichTab = enrichTab_FDR05,
    path_to_plot = here(
        "snRNAseq_mouse",
        "plots",
        "04_clinical_set_enrichment",
        "Gene_set_enrichment_FDR05.pdf"
    )
)


##############################################################
#### Plot pairwise using gene_set_enrichment_plot_complex ####
##############################################################

use_gsepc(
    modeling_results = modeling_results,
    model_type = "pairwise",
    gene_list = gene_list_FDR01,
    enrichTab = prwiseTab_FDR01,
    path_to_plot = here(
        "snRNAseq_mouse",
        "plots",
        "04_clinical_set_enrichment",
        "Gene_set_pairwise_FDR01.pdf"
    )
)

use_gsepc(
    modeling_results = modeling_results,
    model_type = "pairwise",
    gene_list = gene_list_FDR05,
    enrichTab = prwiseTab_FDR05,
    path_to_plot = here(
        "snRNAseq_mouse",
        "plots",
        "04_clinical_set_enrichment",
        "Gene_set_pairwise_FDR05.pdf"
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
