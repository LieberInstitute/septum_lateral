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
# 2.10 kB
lobstr::obj_size(enrichTab_FDR05)
# 2.10 kB
lobstr::obj_size(prwiseTab_FDR01)
# 5.70 kB
lobstr::obj_size(prwiseTab_FDR05)
# 5.70 kB

enrichTab_FDR01
#          OR         Pval test NumSig SetSize       ID model_type fdr_cut
# 1 0.4608680 8.448794e-32   LS    286    1125      all enrichment     0.1
# 2 0.3646524 5.586260e-46   LS    222    1042 positive enrichment     0.1
# 3 4.6841133 8.061127e-11   LS     64      83 negative enrichment     0.1

enrichTab_FDR05
#          OR          Pval test NumSig SetSize       ID model_type fdr_cut
# 1 0.9144459  1.399671e-02   LS   1413    3533      all enrichment     0.1
# 2 0.3400692 3.221221e-108   LS    474    2282 positive enrichment     0.1
# 3 4.4066653 4.230149e-130   LS    939    1251 negative enrichment     0.1

prwiseTab_FDR01
#           OR          Pval          test NumSig SetSize       ID model_type fdr_cut
# 1  0.8441017  5.641178e-03      LS-Astro    550    1125      all   pairwise     0.1
# 2  0.7452277  3.408002e-06      LS-Astro    478    1042 positive   pairwise     0.1
# 3  5.8306375  1.164104e-10      LS-Astro     72      83 negative   pairwise     0.1
# 4  1.4490954  2.242803e-09       LS-Chol    703    1125      all   pairwise     0.1
# 5  1.3873343  3.507366e-07       LS-Chol    641    1042 positive   pairwise     0.1
# 6  2.5412572  9.977551e-05       LS-Chol     62      83 negative   pairwise     0.1
# 7  0.9890193  8.793224e-01        LS-ChP    540    1125      all   pairwise     0.1
# 8  0.8630534  2.138602e-02        LS-ChP    466    1042 positive   pairwise     0.1
# 9  8.8498671  3.991009e-15        LS-ChP     74      83 negative   pairwise     0.1
# 10 0.5061589  1.189967e-28       LS-Endo    437    1125      all   pairwise     0.1
# 11 0.4319198  4.164560e-39       LS-Endo    367    1042 positive   pairwise     0.1
# 12 4.4072244  2.335464e-08       LS-Endo     70      83 negative   pairwise     0.1
# 13 0.8129481  7.455137e-04  LS-Ependymal    492    1125      all   pairwise     0.1
# 14 0.7034108  3.330382e-08  LS-Ependymal    420    1042 positive   pairwise     0.1
# 15 6.9222103  5.555243e-13  LS-Ependymal     72      83 negative   pairwise     0.1
# 16 1.1020732  1.214186e-01        LS-IoC    701    1125      all   pairwise     0.1
# 17 1.0685494  3.188913e-01        LS-IoC    642    1042 positive   pairwise     0.1
# 18 1.6350796  4.328310e-02        LS-IoC     59      83 negative   pairwise     0.1
# 19 0.2013188 2.750013e-130      LS-Micro    241    1125      all   pairwise     0.1
# 20 0.1435482 3.527064e-161      LS-Micro    170    1042 positive   pairwise     0.1
# 21 4.6136307  2.035377e-08      LS-Micro     71      83 negative   pairwise     0.1
# 22 1.5213887  5.914655e-12         LS-MS    584    1125      all   pairwise     0.1
# 23 1.4843475  4.323549e-10         LS-MS    535    1042 positive   pairwise     0.1
# 24 2.0046850  1.731433e-03         LS-MS     49      83 negative   pairwise     0.1
# 25 0.5164756  4.982618e-27      LS-Mural    443    1125      all   pairwise     0.1
# 26 0.4457561  2.239188e-36      LS-Mural    375    1042 positive   pairwise     0.1
# 27 3.6987841  3.615451e-07      LS-Mural     68      83 negative   pairwise     0.1
# 28 0.6193723  7.808788e-15 LS-Neuroblast    434    1125      all   pairwise     0.1
# 29 0.5367588  5.386146e-22 LS-Neuroblast    368    1042 positive   pairwise     0.1
# 30 3.9034934  3.248638e-08 LS-Neuroblast     66      83 negative   pairwise     0.1
# 31 0.4454148  7.031635e-40      LS-Oligo    454    1125      all   pairwise     0.1
# 32 0.3799063  6.430563e-52      LS-Oligo    382    1042 positive   pairwise     0.1
# 33 4.4486383  1.128896e-07      LS-Oligo     72      83 negative   pairwise     0.1
# 34 0.6891235  1.071237e-09        LS-OPC    499    1125      all   pairwise     0.1
# 35 0.6154754  2.180227e-14        LS-OPC    434    1042 positive   pairwise     0.1
# 36 3.1707107  3.975026e-06        LS-OPC     65      83 negative   pairwise     0.1
# 37 0.9646414  5.772511e-01       LS-Sept    441    1125      all   pairwise     0.1
# 38 0.9118632  1.574422e-01       LS-Sept    395    1042 positive   pairwise     0.1
# 39 1.8654508  4.895547e-03       LS-Sept     46      83 negative   pairwise     0.1
# 40 1.2230371  1.023868e-03        LS-Str    647    1125      all   pairwise     0.1
# 41 1.2769614  1.334139e-04        LS-Str    610    1042 positive   pairwise     0.1
# 42 0.7211399  1.526294e-01        LS-Str     37      83 negative   pairwise     0.1
# 43 1.0034017  9.756247e-01       LS-Thal    498    1125      all   pairwise     0.1
# 44 1.0475984  4.659643e-01       LS-Thal    472    1042 positive   pairwise     0.1
# 45 0.5754190  1.980889e-02       LS-Thal     26      83 negative   pairwise     0.1
# 46 1.1854529  6.052296e-03       LS-TNoS    661    1125      all   pairwise     0.1
# 47 1.1181899  8.186624e-02       LS-TNoS    598    1042 positive   pairwise     0.1
# 48 2.6115074  8.968439e-05       LS-TNoS     63      83 negative   pairwise     0.1
# 49 1.2118497  1.593195e-03   LS-TT.IG.SH    610    1125      all   pairwise     0.1
# 50 1.2406630  6.654299e-04   LS-TT.IG.SH    571    1042 positive   pairwise     0.1
# 51 0.9005673  6.613378e-01   LS-TT.IG.SH     39      83 negative   pairwise     0.1

prwiseTab_FDR05
#           OR          Pval          test NumSig SetSize       ID model_type fdr_cut
# 1  1.0199036  5.920441e-01      LS-Astro   1887    3533      all   pairwise     0.1
# 2  0.5011011  2.139634e-55      LS-Astro    849    2282 positive   pairwise     0.1
# 3  4.5413165 3.182017e-114      LS-Astro   1038    1251 negative   pairwise     0.1
# 4  1.2625743  1.043252e-10       LS-Chol   2081    3533      all   pairwise     0.1
# 5  1.1263701  6.906441e-03       LS-Chol   1290    2282 positive   pairwise     0.1
# 6  1.4998078  6.328508e-12       LS-Chol    791    1251 negative   pairwise     0.1
# 7  1.0914319  1.450608e-02        LS-ChP   1774    3533      all   pairwise     0.1
# 8  0.6390076  3.765369e-24        LS-ChP    869    2282 positive   pairwise     0.1
# 9  2.9150396  1.611027e-69        LS-ChP    905    1251 negative   pairwise     0.1
# 10 0.7941486  1.199010e-10       LS-Endo   1765    3533      all   pairwise     0.1
# 11 0.3857778 7.073521e-102       LS-Endo    766    2282 positive   pairwise     0.1
# 12 3.3681802  4.174628e-78       LS-Endo    999    1251 negative   pairwise     0.1
# 13 1.0382606  2.926893e-01  LS-Ependymal   1750    3533      all   pairwise     0.1
# 14 0.4966355  4.408156e-55  LS-Ependymal    755    2282 positive   pairwise     0.1
# 15 4.3045553 4.848238e-116  LS-Ependymal    995    1251 negative   pairwise     0.1
# 16 1.2470697  2.303866e-09        LS-IoC   2286    3533      all   pairwise     0.1
# 17 1.1309870  5.969170e-03        LS-IoC   1433    2282 positive   pairwise     0.1
# 18 1.4433661  1.743277e-09        LS-IoC    853    1251 negative   pairwise     0.1
# 19 0.6050466  7.512007e-45      LS-Micro   1595    3533      all   pairwise     0.1
# 20 0.2307025 6.567584e-222      LS-Micro    563    2282 positive   pairwise     0.1
# 21 3.8241365  1.909192e-89      LS-Micro   1032    1251 negative   pairwise     0.1
# 22 1.3056766  1.046062e-13         LS-MS   1686    3533      all   pairwise     0.1
# 23 1.3256064  1.067562e-10         LS-MS   1103    2282 positive   pairwise     0.1
# 24 1.2213668  5.607248e-04         LS-MS    583    1251 negative   pairwise     0.1
# 25 0.7895550  4.052953e-11      LS-Mural   1763    3533      all   pairwise     0.1
# 26 0.3616345 6.016172e-115      LS-Mural    737    2282 positive   pairwise     0.1
# 27 3.8763838  8.345538e-93      LS-Mural   1026    1251 negative   pairwise     0.1
# 28 0.9585093  2.392583e-01 LS-Neuroblast   1731    3533      all   pairwise     0.1
# 29 0.4698151  4.804679e-64 LS-Neuroblast    754    2282 positive   pairwise     0.1
# 30 3.7412968  7.619737e-97 LS-Neuroblast    977    1251 negative   pairwise     0.1
# 31 0.8418348  1.912765e-06      LS-Oligo   1974    3533      all   pairwise     0.1
# 32 0.4057361  6.093067e-94      LS-Oligo    890    2282 positive   pairwise     0.1
# 33 4.6023610 5.440185e-101      LS-Oligo   1084    1251 negative   pairwise     0.1
# 34 0.9293898  4.161584e-02        LS-OPC   1826    3533      all   pairwise     0.1
# 35 0.4752699  1.107300e-63        LS-OPC    830    2282 positive   pairwise     0.1
# 36 3.5697189  1.011528e-86        LS-OPC    996    1251 negative   pairwise     0.1
# 37 1.1353371  4.749145e-04       LS-Sept   1511    3533      all   pairwise     0.1
# 38 0.8831922  5.595146e-03       LS-Sept    851    2282 positive   pairwise     0.1
# 39 1.7092135  2.042413e-20       LS-Sept    660    1251 negative   pairwise     0.1
# 40 1.2917757  1.131225e-12        LS-Str   2061    3533      all   pairwise     0.1
# 41 1.3838779  1.527718e-13        LS-Str   1372    2282 positive   pairwise     0.1
# 42 1.1042766  8.841187e-02        LS-Str    689    1251 negative   pairwise     0.1
# 43 1.0635334  8.808364e-02       LS-Thal   1609    3533      all   pairwise     0.1
# 44 0.9344307  1.258586e-01       LS-Thal    973    2282 positive   pairwise     0.1
# 45 1.3206202  1.607132e-06       LS-Thal    636    1251 negative   pairwise     0.1
# 46 1.1947641  8.321928e-07       LS-TNoS   2071    3533      all   pairwise     0.1
# 47 0.9795512  6.468749e-01       LS-TNoS   1238    2282 positive   pairwise     0.1
# 48 1.6807665  4.214663e-18       LS-TNoS    833    1251 negative   pairwise     0.1
# 49 1.2304863  6.685425e-09   LS-TT.IG.SH   1915    3533      all   pairwise     0.1
# 50 1.4081411  4.775489e-15   LS-TT.IG.SH   1312    2282 positive   pairwise     0.1
# 51 0.9435630  3.268131e-01   LS-TT.IG.SH    603    1251 negative   pairwise     0.1


######################################################################
#### Plot using as base gene_set_enrichment_plot from spatialLIBD ####
######################################################################

gene_set_enrichment_plot_mod <- function(
        enrichment,
        cols = colorRamp2(
            c(0, 12),
            c("white", "red")),
    path_to_plot,
    plot_name) {
    enrichment$log10_P_thresh <- round(-log10(enrichment$Pval), 2)
    enrichment$log10_P_thresh[which(enrichment$log10_P_thresh > 12)] <- 12
    enrichment$OR_char <- as.character(round(enrichment$OR, 2))
    enrichment$OR_char[enrichment$log10_P_thresh < 3] <- ""
    make_wide <- function(var = "OR_char") {
        res <- reshape(enrichment,
            idvar = "ID", timevar = "test",
            direction = "wide", drop = colnames(enrichment)[!colnames(enrichment) %in%
                c("ID", "test", var)], sep = "_mypattern_"
        )[, -1, drop = FALSE]
        colnames(res) <- gsub(".*_mypattern_", "", colnames(res))
        rownames(res) <- unique(enrichment$ID)
        res <- res[, levels(as.factor(enrichment$test))]
        t(res)
    }
    wide_or <- make_wide("OR_char")
    wide_p <- make_wide("log10_P_thresh")

    pdf(path_to_plot, height = 4, width = 6)
    plot_gs <- Heatmap(wide_p,
        row_title_side = "left",
        rect_gp = gpar(col = "white", lwd = 2),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(wide_or[i, j], x, y, gp = gpar(fontsize = 10))
        },
        col = cols,
        column_title = plot_name,
        heatmap_legend_param = list(
            title = "-log10(p-val)",
            at = c(0, 2, 4, 6, 8, 10, 12)
        )
    )
    print(plot_gs)
    dev.off()
}

gene_set_enrichment_plot_mod(enrichment = prwiseTab_FDR05, path_to_plot = here("snRNAseq_mouse/", "plots/", "04_clinical_set_enrichment/", "Gene_set_pairwise_FDR05.pdf"), plot_name = "Pairwise FDR < 0.05")

gene_set_enrichment_plot_mod(enrichment = prwiseTab_FDR01, path_to_plot = here("snRNAseq_mouse/", "plots/", "04_clinical_set_enrichment/", "Gene_set_pairwise_FDR01.pdf"), plot_name = "Pairwise FDR < 0.01")

#####################################
#### Reproducibility information ####
#####################################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
