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
tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                    "#CDCC5D", "#6DCCDA")
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")

# From http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
#      (since there are more than 30 clusters, but <=38)
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7")

# plotExpressionCustom for nicer aesthetics
source("/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/analyses_sn/plotExpressionCustom.R")



### TrkB KO experiment ==========
  # original project dir: /dcl01/lieber/ajaffe/Keri/TrkBKO/
  # (DE list csv copied over into this dir)

DE.TrkB.KO <- read.csv("/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/processed_data/tables/TrkBKO_DEgenes_FDR05.csv",
                       sep=",")

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
    #FALSE  TRUE 
    # 1348  2402 (checks out with the volcano plot)




### First: Are any DE genes cell type markers? ========================
load(here("snRNAseq_mouse", "processed_data","SCE",
          "markers-stats_LS-n4_findMarkers_35cellTypes.rda"), verbose=T)

length(intersect(DE.TrkB.KO$Symbol, rowData(sce.ls)$gene_name))
length(setdiff(DE.TrkB.KO$Symbol, rowData(sce.ls)$gene_name))
    #[1] 220
length(intersect(DE.TrkB.KO$ensemblID, rowData(sce.ls)$gene_id))
length(setdiff(DE.TrkB.KO$ensemblID, rowData(sce.ls)$gene_id))
    #[1] 214

setdiff(DE.TrkB.KO$Symbol, rowData(sce.ls)$gene_name)

pdf(here("snRNAseq_mouse","plots","volcanoPlot_intersectingDEgenes-vs-SCE.pdf"))
plot(x=DE.TrkB.KO$logFC, y=(-log10(DE.TrkB.KO$P.Value)),
     pch=16, col=ifelse(!(DE.TrkB.KO$Symbol %in% rowData(sce.ls)$gene_name),
                        "red", "black"),
     main="TrkB-KO significant genes not in SCE (red)")
dev.off()



# Create marker list (from PW tests)
markerList.t.pw <- lapply(markers.ls.t.pw, function(x){
  rownames(x)[x$FDR < 0.05 & x$non0median == TRUE]
}
)

# Change to gene symbols
markerList.t.pw <- lapply(markerList.t.pw, function(x){
  rowData(sce.ls)$gene_name[match(x, rowData(sce.ls)$gene_id)]
})

# For reference:
lengths(markerList.t.pw)
    # ambig.glial_A  ambig.glial_B   Aqp4.Rbpms_A   Aqp4.Rbpms_B        Astro_A 
    #           142             22              8              0             30 
    #       Astro_B  Astro.OPC_COP   drop.doublet    drop.lowNTx         Endo_A 
    #             0             49              0             33              0 
    #        Endo_B        Excit_A        Excit_B        Excit_C        Excit_D 
    #           511             18             36             44             74 
    #       Excit_E        Excit_F        Excit_G        Inhib_A        Inhib_B 
    #            32             36              2              6             35 
    #       Inhib_C        Inhib_D        Inhib_E        Inhib_F        Inhib_G 
    #            31              1             44             31             42 
    #       Inhib_H          Micro        Mural_A        Mural_B Neuron.mixed_A 
    #            40            148             22            136              2 
    #Neuron.Ppp1r1b        Oligo_A        Oligo_B            OPC        OPC_COP 
    #           195              9            121             45            115 

sapply(markerList.t.pw, function(x){
  round(length(intersect(x, DE.TrkB.KO$Symbol)) / length(x) * 100, 2)
})
    # ambig.glial_A  ambig.glial_B   Aqp4.Rbpms_A   Aqp4.Rbpms_B        Astro_A 
    #          7.75          45.45          12.50            NaN          30.00 
    #       Astro_B  Astro.OPC_COP   drop.doublet    drop.lowNTx         Endo_A 
    #           NaN          30.61            NaN          51.52            NaN 
    #        Endo_B        Excit_A        Excit_B        Excit_C        Excit_D 
    #         36.20           5.56          41.67          15.91          20.27 
    #       Excit_E        Excit_F        Excit_G        Inhib_A        Inhib_B 
    #         18.75          27.78           0.00          16.67          14.29 
    #       Inhib_C        Inhib_D        Inhib_E        Inhib_F        Inhib_G 
    #         19.35           0.00          18.18          41.94          26.19 
    #       Inhib_H          Micro        Mural_A        Mural_B Neuron.mixed_A 
    #         37.50          64.86          40.91          24.26          50.00 
    #Neuron.Ppp1r1b        Oligo_A        Oligo_B            OPC        OPC_COP 
    #         20.51          33.33          33.06          11.11          22.61


# Split by direction of FC
DE.up <- DE.TrkB.KO[DE.TrkB.KO$logFC > 0, ]
DE.down <- DE.TrkB.KO[DE.TrkB.KO$logFC < 0, ]


sapply(markerList.t.pw, function(x){
  length(intersect(x, DE.up$Symbol)) / length(intersect(x, DE.down$Symbol)) 
})

DEmarkersNums <- matrix(c(0,0,0,0), ncol=1)
for(i in names(markerList.t.pw)){
  #print(paste0("total DE genes in marker set for ",i,": "))
  nMarkers <- length(markerList.t.pw[[i]])
  total <- length(intersect(markerList.t.pw[[i]], DE.TrkB.KO$Symbol))
  # For up-regulated:
  up <- length(intersect(markerList.t.pw[[i]], DE.TrkB.KO$Symbol[DE.TrkB.KO$logFC > 0]))
  # For down-regulated:
  down <- length(intersect(markerList.t.pw[[i]], DE.TrkB.KO$Symbol[DE.TrkB.KO$logFC < 0]))
  vec <- c(nMarkers, total, up, down)
  
  DEmarkersNums <- cbind(DEmarkersNums, vec)
  }

DEmarkersNums <- DEmarkersNums[ ,-1]
colnames(DEmarkersNums) <- names(markerList.t.pw)
rownames(DEmarkersNums) <- c("nMarkers","total.intersecting","upRegulated","downRegulated")

DEmarkersNums
    #                    ambig.glial_A ambig.glial_B Aqp4.Rbpms_A Aqp4.Rbpms_B Astro_A Astro_B
    # nMarkers                     142            22            8            0      30       0
    # total.intersecting            11            10            1            0       9       0
    # upRegulated                    9            10            1            0       9       0
    # downRegulated                  2             0            0            0       0       0

    #                    Astro.OPC_COP drop.doublet drop.lowNTx Endo_A Endo_B Excit_A Excit_B Excit_C
    # nMarkers                      49            0          33      0    511      18      36      44
    # total.intersecting            15            0          17      0    185       1      15       7
    # upRegulated                   15            0           5      0    175       0       0       0
    # downRegulated                  0            0          12      0     10       1      15       7

    #                    Excit_D Excit_E Excit_F Excit_G Inhib_A Inhib_B Inhib_C Inhib_D Inhib_E Inhib_F
    # nMarkers                74      32      36       2       6      35      31       1      44      31
    # total.intersecting      15       6      10       0       1       5       6       0       8      13
    # upRegulated              1       2       1       0       0       0       0       0       0       0
    # downRegulated           14       4       9       0       1       5       6       0       8      13

    #                    Inhib_G Inhib_H Micro Mural_A Mural_B Neuron.mixed_A Neuron.Ppp1r1b Oligo_A
    # nMarkers                42      40   148      22     136              2            195       9
    # total.intersecting      11      15    96       9      33              1             40       3
    # upRegulated              0       1    92       8      29              0             27       3
    # downRegulated           11      14     4       1       4              1             13       0

    #                    Oligo_B OPC OPC_COP
    # nMarkers               121  45     115
    # total.intersecting      40   5      26
    # upRegulated             38   1      16
    # downRegulated            2   4      10


# Intersecting gene IDs (for Lionel to explore) ===
for(i in names(markerList.t.pw)){
  print(paste0(i," intersecting (strict-PW) markers:"))
  # For up-regulated:
  cat("Up-regulated:\n")
  print(intersect(markerList.t.pw[[i]], DE.TrkB.KO$Symbol[DE.TrkB.KO$logFC > 0]))
  # For down-regulated:
  cat("\nDown-regulated:\n")
  print(intersect(markerList.t.pw[[i]], DE.TrkB.KO$Symbol[DE.TrkB.KO$logFC < 0]))
  cat("\n\n\n")
}



## Enrichment tests (with the '1vAll' numbers) ===






## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
    #[1] "2022-05-10 15:48:33 EDT"
proc.time()
    #     user    system   elapsed 
    #    148.449    5.915 5428.779 
options(width = 120)
session_info()
    #─ Session info ────────────────────────────────────────────────────────────
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


