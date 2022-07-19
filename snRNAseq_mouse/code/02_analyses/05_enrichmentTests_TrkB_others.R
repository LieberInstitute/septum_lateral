### mouse LS snRNA-seq analysis
### How do cell type markers (or marker lists) intersect with
# gene sets of particular interest?
### qsub -l bluejay,mf=60G,h_vmem=64G,h_fsize=40G
### Initiated LAR,MNT 03May2022

library(SingleCellExperiment)
library(scater)
library(scran)
library(scry)
library(batchelor)
library(BiocParallel)
library(jaffelab)
library(gridExtra)
library(here)
library(sessioninfo)

here()

### Palette taken from `scater`
tableau10medium <- c(
    "#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
    "#CDCC5D", "#6DCCDA"
)
tableau20 <- c(
    "#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
    "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
    "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
    "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5"
)

# From http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
#      (since there are more than 30 clusters, but <=38)
cbPalette <- c(
    "#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7"
)

# plotExpressionCustom for nicer aesthetics
source("/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/analyses_sn/plotExpressionCustom.R")



### TrkB KO experiment ==========
# original project dir: /dcl01/lieber/ajaffe/Keri/TrkBKO/
# (DE list csv copied over into this dir)

DE.TrkB.KO <- read.csv("/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/processed_data/tables/TrkBKO_DEgenes_FDR05.csv",
    sep = ","
)

head(DE.TrkB.KO)
#                       X             gencodeID Symbol    logFC        t
# 1  ENSMUSG00000060802.8  ENSMUSG00000060802.8    B2m 3.390829 20.05293
# 2  ENSMUSG00000069874.7  ENSMUSG00000069874.7  Irgm2 3.811124 19.80172
# 3  ENSMUSG00000046879.6  ENSMUSG00000046879.6  Irgm1 4.346880 19.65682
# 4 ENSMUSG00000040253.15 ENSMUSG00000040253.15   Gbp7 3.474357 19.52114
# 5  ENSMUSG00000067212.8  ENSMUSG00000067212.8 H2-T23 3.999794 19.48655
# 6  ENSMUSG00000078853.8  ENSMUSG00000078853.8   Igtp 4.561851 19.04770

#        P.Value    adj.P.Val         B          ensemblID Length      gene_type
# 1 1.642178e-08 7.538354e-05 10.375556 ENSMUSG00000060802    860 protein_coding
# 2 1.826103e-08 7.538354e-05 10.302562 ENSMUSG00000069874   3690 protein_coding
# 3 1.942555e-08 7.538354e-05 10.243430 ENSMUSG00000046879   5810 protein_coding
# 4 2.059150e-08 7.538354e-05 10.187026 ENSMUSG00000040253   5614 protein_coding
# 5 2.090102e-08 7.538354e-05 10.170358 ENSMUSG00000067212   2889 protein_coding
# 6 2.531390e-08 7.538354e-05  9.977743 ENSMUSG00000078853   2022 protein_coding

#   EntrezID  AveExpr
# 1    12010 8.238995
# 2    54396 5.594922
# 3    15944 5.908829
# 4   229900 5.396657
# 5    15040 5.503027
# 6    16145 5.691325

# How many downregulated/upregulated?
table(DE.TrkB.KO$logFC > 0)
# FALSE  TRUE
# 1348  2402 (checks out with the volcano plot)




### First: Are any DE genes cell type markers? ========================
load(here("snRNAseq_mouse", "processed_data", "SCE", "sce_updated_LS.rda"), verbose = T)
# sce.ls, annotationTab.ls, cell_colors.ls

load(here(
    "snRNAseq_mouse", "processed_data", "SCE",
    "markers-stats_LS-n4_findMarkers_33cellTypes.rda"
), verbose = T)

length(intersect(DE.TrkB.KO$Symbol, rowData(sce.ls)$gene_name))
# [1] 3524
length(setdiff(DE.TrkB.KO$Symbol, rowData(sce.ls)$gene_name))
# [1] 220

# By Ensembl ID
length(intersect(DE.TrkB.KO$ensemblID, rowData(sce.ls)$gene_id))
# [1] 3536
length(setdiff(DE.TrkB.KO$ensemblID, rowData(sce.ls)$gene_id))
# [1] 214

sort(setdiff(DE.TrkB.KO$Symbol, rowData(sce.ls)$gene_name))
# A lot of GmXXXXXX and XXXXXXXRik genes

pdf(here("snRNAseq_mouse", "plots", "volcanoPlot_intersectingDEgenes-vs-SCE.pdf"))
plot(
    x = DE.TrkB.KO$logFC, y = (-log10(DE.TrkB.KO$P.Value)),
    pch = 16, col = ifelse(!(DE.TrkB.KO$Symbol %in% rowData(sce.ls)$gene_name),
        "red", "black"
    ),
    main = "TrkB-KO significant genes not in SCE (red)"
)
dev.off()



# Create marker list (from PW tests)
markerList.t.pw <- lapply(markers.ls.t.pw, function(x) {
    rownames(x)[x$FDR < 0.05 & x$non0median == TRUE]
})

# Change to gene symbols
markerList.t.pw <- lapply(markerList.t.pw, function(x) {
    rowData(sce.ls)$gene_name[match(x, rowData(sce.ls)$gene_id)]
})

# For reference:
lengths(markerList.t.pw)
# Aqp4.Rbpms      Astro       Endo  Ependymal    Excit_A    Excit_B    Excit_C
#       493         98        328        270         46         25         16
#   Excit_D    Excit_E    Excit_F    Inhib_A    Inhib_B    Inhib_C    Inhib_D
#        69         25         14          2         29         14          0
#   Inhib_E    Inhib_F    Inhib_G    Inhib_H    Inhib_I    Inhib_J    Inhib_K
#        43         26         41         31          0          6          1
#   Inhib_L    Inhib_M    Inhib_N    Inhib_O    Inhib_P    Inhib_Q    Inhib_R
#        24         36         14         17          2          4         20
#     Micro      Mural      Oligo        OPC    OPC_COP
#       180         17        118         50        141

sapply(markerList.t.pw, function(x) {
    round(length(intersect(x, DE.TrkB.KO$Symbol)) / length(x) * 100, 2)
})
# Aqp4.Rbpms      Astro       Endo  Ependymal    Excit_A    Excit_B    Excit_C
#     17.65      33.67      32.01      21.48      10.87      32.00       6.25
#   Excit_D    Excit_E    Excit_F    Inhib_A    Inhib_B    Inhib_C    Inhib_D
#     21.74      12.00      28.57       0.00      10.34      21.43        NaN
#   Inhib_E    Inhib_F    Inhib_G    Inhib_H    Inhib_I    Inhib_J    Inhib_K
#     16.28      30.77      31.71      41.94        NaN      16.67       0.00
#   Inhib_L    Inhib_M    Inhib_N    Inhib_O    Inhib_P    Inhib_Q    Inhib_R
#     16.67       8.33      35.71      11.76       0.00      25.00      15.00
#     Micro      Mural      Oligo        OPC    OPC_COP
#     63.33      52.94      40.68      16.00      29.08


# Split by direction of FC
DE.up <- DE.TrkB.KO[DE.TrkB.KO$logFC > 0, ]
DE.down <- DE.TrkB.KO[DE.TrkB.KO$logFC < 0, ]

# Ratio of up-regulated to down-regulated markers
sapply(markerList.t.pw, function(x) {
    round(length(intersect(x, DE.up$Symbol)) / length(intersect(x, DE.down$Symbol)), 2)
})
## Note that if none are up-regulated, will be 0 (even if some or down-regulated)

# This is interesting...!:
intersect(markerList.t.pw[["Inhib_N"]], DE.down$Symbol)
# [1] "Adra1b"        "Kit"           "A830082K12Rik" "Reln"
# [5] "Slc2a13"


DEmarkersNums <- matrix(c(0, 0, 0, 0), ncol = 1)
for (i in names(markerList.t.pw)) {
    # print(paste0("total DE genes in marker set for ",i,": "))
    nMarkers <- length(markerList.t.pw[[i]])
    total <- length(intersect(markerList.t.pw[[i]], DE.TrkB.KO$Symbol))
    # For up-regulated:
    up <- length(intersect(markerList.t.pw[[i]], DE.TrkB.KO$Symbol[DE.TrkB.KO$logFC > 0]))
    # For down-regulated:
    down <- length(intersect(markerList.t.pw[[i]], DE.TrkB.KO$Symbol[DE.TrkB.KO$logFC < 0]))
    vec <- c(nMarkers, total, up, down)

    DEmarkersNums <- cbind(DEmarkersNums, vec)
}

DEmarkersNums <- DEmarkersNums[, -1]
colnames(DEmarkersNums) <- names(markerList.t.pw)
rownames(DEmarkersNums) <- c("nMarkers", "total.intersecting", "upRegulated", "downRegulated")

DEmarkersNums
#                   Aqp4.Rbpms Astro Endo Ependymal Excit_A Excit_B Excit_C
# nMarkers                  493    98  328       270      46      25      16
# total.intersecting         87    33  105        58       5       8       1
# upRegulated                73    33   91        44       1       0       0
# downRegulated              14     0   14        14       4       8       1

#                    Excit_D Excit_E Excit_F Inhib_A Inhib_B Inhib_C Inhib_D
# nMarkers                69      25      14       2      29      14       0
# total.intersecting      15       3       4       0       3       3       0
# upRegulated              1       1       0       0       0       0       0
# downRegulated           14       2       4       0       3       3       0

#                    Inhib_E Inhib_F Inhib_G Inhib_H Inhib_I Inhib_J Inhib_K
# nMarkers                43      26      41      31       0       6       1
# total.intersecting       7       8      13      13       0       1       0
# upRegulated              1       0       2       1       0       0       0
# downRegulated            6       8      11      12       0       1       0

#                    Inhib_L Inhib_M Inhib_N Inhib_O Inhib_P Inhib_Q Inhib_R
# nMarkers                24      36      14      17       2       4      20
# total.intersecting       4       3       5       2       0       1       3
# upRegulated              1       2       0       0       0       0       2
# downRegulated            3       1       5       2       0       1       1

#                    Micro Mural Oligo OPC OPC_COP
# nMarkers             180    17   118  50     141
# total.intersecting   114     9    48   8      41
# upRegulated          111     9    45   6      31
# downRegulated          3     0     3   2      10


# Intersecting gene IDs (for Lionel to explore) ===
for (i in names(markerList.t.pw)) {
    print(paste0(i, " intersecting (strict-PW) markers:"))
    # For up-regulated:
    cat("Up-regulated:\n")
    print(intersect(markerList.t.pw[[i]], DE.TrkB.KO$Symbol[DE.TrkB.KO$logFC > 0]))
    # For down-regulated:
    cat("\nDown-regulated:\n")
    print(intersect(markerList.t.pw[[i]], DE.TrkB.KO$Symbol[DE.TrkB.KO$logFC < 0]))
    cat("\n\n\n")
}



## Enrichment tests (with the '1vAll' numbers) ===
# (skip for now - may want to run at a later time)

# sapply(markers.ls.t.1vAll, names)
#
# # Take just the enriched set
# markerList.t.1vAll <- lapply(markers.ls.t.1vAll, function(x){
#   rownames(x[[2]])[ x[[2]]$log.FDR < log(0.05) & x[[2]]$non0median==TRUE ]
# }
# )
#
# # Change to gene symbols
# markerList.t.1vAll <- lapply(markerList.t.1vAll, function(x){
#   rowData(sce.ls)$gene_name[match(x, rowData(sce.ls)$gene_id)]
# })
#
# ## From Andrew's code (at /dcl01/lieber/ajaffe/Keri/TrkBKO)
#     # Load the RSE for the DE analyses
#     load("/dcl01/lieber/ajaffe/Keri/TrkBKO/preprocessed_data/rse_gene_TrkB_KO_LS_n8.Rdata",
#          verbose=T)
#         # rse_gene, getRPKM
#
#     ## get phenotype data
#     pd = read.csv("/dcl01/ajaffe/data/Nina/Keri/SunHong_102318/sampleList.txt",as.is=TRUE)
#     rownames(pd) = paste0("Sample.", pd$SampleID)
#     colData(rse_gene) = DataFrame(cbind(pd[colnames(rse_gene),], data.frame(colData(rse_gene))))
#
#     geneExprs = log2(getRPKM(rse_gene,"Length")+1)
#
#     ## check Bdnf
#     bdnf = geneExprs[rownames(rse_gene)[which(rowData(rse_gene)$Symbol == "Bdnf")],]
#     boxplot(bdnf ~ rse_gene$Condition,las=3, ylab = "Bdnf: log2(RPKM+1)")

#         # - more would need to be done (set up and run Fisher's tests)



### Can export previous GO enrichment test gene lists?? ========

# As per /dcl01/lieber/ajaffe/Keri/TrkBKO/de_analysis.R, there's no .rda
#   (or any file) at /dcl01/lieber/ajaffe/Keri/TrkBKO/rdas/

# What about this?
load("/dcl01/lieber/ajaffe/Keri/TrkBKO/gene_sets.rda", verbose = T)

table(goDf$ONTOLOGY)
#  BP   CC   MF
# 2018  240  259
table(goDf$p.adj <= 0.05)
# FALSE  TRUE
#  181  2336
table(goDf$p.adj <= 0.10) # all TRUE - so this is probably the results we want


# Check these are all DE genes
table(unlist(strsplit(goDf$geneID[1], split = "/")) %in% DE.TrkB.KO$Symbol)
# all 103 (these are the intersecting b/tw DE genes) (see GeneRatio)

# B/tw those top to significant terms
table(unlist(strsplit(goDf$geneID[1], split = "/")) %in% unlist(strsplit(goDf$geneID[2], split = "/")))
# Same 103 genes in both terms

# Now check top 100 in the significant GO terms list
unlist(sapply(c(1:100), function(x) {
    table(unlist(strsplit(goDf$geneID[x], split = "/")) %in% DE.TrkB.KO$Symbol)
}))
# Looks like a small number FALSE... this is probably bc the input for clusterProfiler
#     is the Entrez IDs, and the discrepancy/redundancy with gene Symbols, etc


## For future exploration, since this is a lot of terms and genes...
#      can pull out those DE genes in a given significant term of interest as follows:
DEgenesInTerm <- unlist(strsplit(goDf$geneID[goDf$Description == "regulation of synaptic plasticity"],
    split = "/"
))

DEgenesInTerm
# [1] "Rims2"    "Syt4"     "Prkcz"    "Grm5"     "Lrrtm1"   "Shisa9"
# [7] "Ppp1r9a"  "Dgki"     "Kcnb1"    "Rab3a"    "Snca"     "Jph4"
# [13] "Grik2"    "Snap47"   "Mapk1"    "Jph3"     "Ncdn"     "Syp"
# [19] "Camk2b"   "Syngr1"   "Unc13a"   "Grin1"    "Rasgrf1"  "Rims1"
# [25] "Ywhag"    "Egr1"     "Rab3gap1" "Grid2"    "Adgrb1"   "Npas4"
# [31] "Dlg4"     "Rapgef2"  "Rims3"    "Hras"     "Arf1"     "Fgf14"
# [37] "Stau2"    "Kit"      "Stxbp1"   "Crhr2"    "Snap25"   "Shisa7"
# [43] "Nsg1"     "Rims4"    "Vamp2"    "Rasgrf2"  "Syt7"     "Ppfia3"
# [49] "Reln"     "Plk2"     "Slc8a2"   "Nos1"     "Cplx2"    "Pak1"
# [55] "Vgf"      "Nrgn"     "Atp2b2"   "Nsmf"     "Sipa1l1"  "Nlgn3"
# [61] "Mapt"     "Baiap2"   "Gsk3b"


## Save the 'goDf' without the long gene lists ($geneID) as a CSV and commit, so can explore interactively,
#      THEN print out those genes in R
goDf2print <- goDf[, c(1:9)]

write.csv(goDf2print,
    file = here("snRNAseq_mouse", "processed_data", "tables", "GOterms_TrkB-KO-DE_withoutGeneLists.csv"),
    quote = F, row.names = F
)

# jk there are commas in some Descriptions
goDf2print$Description <- gsub(",", " -", goDf2print$Description)


## Reproducibility information ====
print("Reproducibility information:")
Sys.time()
# [1] "2022-05-10 15:48:33 EDT"
proc.time()
#     user    system   elapsed
#    148.449    5.915 5428.779
options(width = 120)
session_info()
# ─ Session info ────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.1.2 Patched (2021-11-04 r81138)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2022-05-10
# pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
#
# ─ Packages ────────────────────────────────────────────────────────────────
# package              * version  date (UTC) lib source
# assertthat             0.2.1    2019-03-21 [2] CRAN (R 4.1.0)
# batchelor            * 1.10.0   2021-10-26 [1] Bioconductor
# beachmat               2.10.0   2021-10-26 [2] Bioconductor
# beeswarm               0.4.0    2021-06-01 [2] CRAN (R 4.1.2)
# Biobase              * 2.54.0   2021-10-26 [2] Bioconductor
# BiocGenerics         * 0.40.0   2021-10-26 [2] Bioconductor
# BiocNeighbors          1.12.0   2021-10-26 [2] Bioconductor
# BiocParallel         * 1.28.3   2021-12-09 [2] Bioconductor
# BiocSingular           1.10.0   2021-10-26 [2] Bioconductor
# bitops                 1.0-7    2021-04-24 [2] CRAN (R 4.1.0)
# bluster                1.4.0    2021-10-26 [2] Bioconductor
# cli                    3.3.0    2022-04-25 [2] CRAN (R 4.1.2)
# cluster                2.1.3    2022-03-28 [3] CRAN (R 4.1.2)
# colorspace             2.0-3    2022-02-21 [2] CRAN (R 4.1.2)
# crayon                 1.5.1    2022-03-26 [2] CRAN (R 4.1.2)
# DBI                    1.1.2    2021-12-20 [2] CRAN (R 4.1.2)
# DelayedArray           0.20.0   2021-10-26 [2] Bioconductor
# DelayedMatrixStats     1.16.0   2021-10-26 [2] Bioconductor
# dplyr                  1.0.9    2022-04-28 [2] CRAN (R 4.1.2)
# dqrng                  0.3.0    2021-05-01 [2] CRAN (R 4.1.2)
# edgeR                  3.36.0   2021-10-26 [2] Bioconductor
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
# fansi                  1.0.3    2022-03-24 [2] CRAN (R 4.1.2)
# fs                     1.5.2    2021-12-08 [2] CRAN (R 4.1.2)
# gargle                 1.2.0    2021-07-02 [2] CRAN (R 4.1.0)
# generics               0.1.2    2022-01-31 [2] CRAN (R 4.1.2)
# GenomeInfoDb         * 1.30.1   2022-01-30 [2] Bioconductor
# GenomeInfoDbData       1.2.7    2021-11-01 [2] Bioconductor
# GenomicRanges        * 1.46.1   2021-11-18 [2] Bioconductor
# ggbeeswarm             0.6.0    2017-08-07 [2] CRAN (R 4.1.2)
# ggplot2              * 3.3.6    2022-05-03 [2] CRAN (R 4.1.2)
# ggrepel                0.9.1    2021-01-15 [2] CRAN (R 4.1.0)
# glue                   1.6.2    2022-02-24 [2] CRAN (R 4.1.2)
# googledrive            2.0.0    2021-07-08 [2] CRAN (R 4.1.0)
# gridExtra            * 2.3      2017-09-09 [2] CRAN (R 4.1.0)
# gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
# here                 * 1.0.1    2020-12-13 [2] CRAN (R 4.1.2)
# igraph                 1.3.1    2022-04-20 [2] CRAN (R 4.1.2)
# IRanges              * 2.28.0   2021-10-26 [2] Bioconductor
# irlba                  2.3.5    2021-12-06 [2] CRAN (R 4.1.2)
# jaffelab             * 0.99.31  2021-12-13 [1] Github (LieberInstitute/jaffelab@2cbd55a)
# lattice                0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
# lifecycle              1.0.1    2021-09-24 [2] CRAN (R 4.1.2)
# limma                  3.50.3   2022-04-07 [2] Bioconductor
# locfit                 1.5-9.5  2022-03-03 [2] CRAN (R 4.1.2)
# magrittr               2.0.3    2022-03-30 [2] CRAN (R 4.1.2)
# Matrix                 1.4-1    2022-03-23 [3] CRAN (R 4.1.2)
# MatrixGenerics       * 1.6.0    2021-10-26 [2] Bioconductor
# matrixStats          * 0.62.0   2022-04-19 [2] CRAN (R 4.1.2)
# metapod                1.2.0    2021-10-26 [2] Bioconductor
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
# pillar                 1.7.0    2022-02-01 [2] CRAN (R 4.1.2)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# purrr                  0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# R6                     2.5.1    2021-08-19 [2] CRAN (R 4.1.2)
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.1.2)
# RColorBrewer           1.1-3    2022-04-03 [2] CRAN (R 4.1.2)
# Rcpp                   1.0.8.3  2022-03-17 [2] CRAN (R 4.1.2)
# RCurl                  1.98-1.6 2022-02-08 [2] CRAN (R 4.1.2)
# ResidualMatrix         1.4.0    2021-10-26 [1] Bioconductor
# rlang                  1.0.2    2022-03-04 [2] CRAN (R 4.1.2)
# rprojroot              2.0.3    2022-04-02 [2] CRAN (R 4.1.2)
# rsvd                   1.0.5    2021-04-16 [2] CRAN (R 4.1.2)
# S4Vectors            * 0.32.4   2022-03-24 [2] Bioconductor
# ScaledMatrix           1.2.0    2021-10-26 [2] Bioconductor
# scales                 1.2.0    2022-04-13 [2] CRAN (R 4.1.2)
# scater               * 1.22.0   2021-10-26 [2] Bioconductor
# scran                * 1.22.1   2021-11-14 [2] Bioconductor
# scry                 * 1.6.0    2021-10-26 [2] Bioconductor
# scuttle              * 1.4.0    2021-10-26 [2] Bioconductor
# segmented              1.3-4    2021-04-22 [1] CRAN (R 4.1.2)
# sessioninfo          * 1.2.2    2021-12-06 [2] CRAN (R 4.1.2)
# SingleCellExperiment * 1.16.0   2021-10-26 [2] Bioconductor
# sparseMatrixStats      1.6.0    2021-10-26 [2] Bioconductor
# statmod                1.4.36   2021-05-10 [2] CRAN (R 4.1.0)
# SummarizedExperiment * 1.24.0   2021-10-26 [2] Bioconductor
# tibble                 3.1.7    2022-05-03 [2] CRAN (R 4.1.2)
# tidyselect             1.1.2    2022-02-21 [2] CRAN (R 4.1.2)
# utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
# vctrs                  0.4.1    2022-04-13 [2] CRAN (R 4.1.2)
# vipor                  0.4.5    2017-03-22 [2] CRAN (R 4.1.2)
# viridis                0.6.2    2021-10-13 [2] CRAN (R 4.1.2)
# viridisLite            0.4.0    2021-04-13 [2] CRAN (R 4.1.0)
# withr                  2.5.0    2022-03-03 [2] CRAN (R 4.1.2)
# XVector                0.34.0   2021-10-26 [2] Bioconductor
# zlibbioc               1.40.0   2021-10-26 [2] Bioconductor
#
# [1] /users/ntranngu/R/4.1.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/R/4.1.x/lib64/R/library
#
# ───────────────────────────────────────────────────────────────────────────
