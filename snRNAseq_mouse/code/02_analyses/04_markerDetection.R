### mouse LS snRNA-seq analysis
### Feature selection with deviance & clustering
###     qsub -l bluejay,mf=76G,h_vmem=80G,h_fsize=40G
### Initiated LAR,ST,MNT 07Apr2022

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


## ===


load(here("snRNAseq_mouse", "processed_data","SCE", "sce_updated_LS.rda"), verbose=T)
    # sce.ls, annotationTab.ls, cell_colors.ls

table(sce.ls$cellType.final)
    #  Aqp4.Rbpms              Astro       drop.doublet drop.likelyDoublet 
    #         486               5682                150                573 
    # drop.lowNTx               Endo          Ependymal            Excit_A 
    #         180                158                 65                485 
    #     Excit_B            Excit_C            Excit_D            Excit_E 
    #         670                232                 67                 31 
    #     Excit_F            Excit_G            Inhib_A            Inhib_B 
    #         202                728                797                760 
    #     Inhib_C            Inhib_D            Inhib_E            Inhib_F 
    #         191               1075                306               4727 
    #     Inhib_G            Inhib_H            Inhib_I            Inhib_J 
    #        1524                 86                111                402 
    #     Inhib_K            Inhib_L            Inhib_M            Inhib_N 
    #          36                167                 45                 69 
    #     Inhib_O            Inhib_P            Inhib_Q              Micro 
    #         190                 36                 68                166 
    #       Mural       Neuron.mixed              Oligo                OPC 
    #         222                 73               1583                464 
    #     OPC_COP 
    #          53 

## doubletScore distributions / cluster?
cellClust.idx <- splitit(sce.ls$cellType.final)
sapply(cellClust.idx, function(x){round(quantile(sce.ls$doubletScore[x]), 2)})
    #      Aqp4.Rbpms Astro drop.doublet drop.likelyDoublet drop.lowNTx Endo
    # 0%         0.00  0.00         0.06               0.00        0.00 0.00
    # 25%        0.04  0.06         3.78               0.62        0.01 0.04
    # 50%        0.10  0.18         4.83               2.01        0.03 0.06
    # 75%        0.46  0.48         6.72               2.94        0.06 0.17
    # 100%      10.39 14.84         9.39              11.38        7.93 5.98
    #      Ependymal Excit_A Excit_B Excit_C Excit_D Excit_E Excit_F Excit_G Inhib_A
    # 0%        0.02    0.00    0.02    0.03    0.04    0.12    0.02    0.04    0.01
    # 25%       0.36    0.04    0.15    0.14    0.38    0.61    0.13    0.17    0.13
    # 50%       0.49    0.07    0.39    0.46    1.28    0.78    0.23    0.30    0.25
    # 75%       0.53    0.27    1.31    1.30    1.39    1.21    0.52    1.03    0.62
    # 100%      7.10   12.53   10.41    6.82    9.93    7.56    7.62    9.25   10.44
    #      Inhib_B Inhib_C Inhib_D Inhib_E Inhib_F Inhib_G Inhib_H Inhib_I Inhib_J
    # 0%      0.00    0.01    0.01    0.00    0.00    0.00    0.09    0.00    0.00
    # 25%     0.03    0.23    0.17    0.02    0.15    0.13    0.30    0.12    0.13
    # 50%     0.07    0.43    0.36    0.06    0.29    0.27    0.52    0.23    0.28
    # 75%     0.78    1.08    0.93    0.52    0.90    0.74    0.91    0.87    0.67
    # 100%   11.46    5.98   11.49    6.60   10.09   11.42    8.40    5.92   11.43
    #      Inhib_K Inhib_L Inhib_M Inhib_N Inhib_O Inhib_P Inhib_Q Micro Mural
    # 0%      0.01    0.01    0.06    0.13    0.02    0.04    0.07  0.00  0.00
    # 25%     0.26    0.11    0.32    0.28    0.17    0.19    0.16  0.06  0.04
    # 50%     0.39    0.18    0.42    0.40    0.25    0.34    0.23  0.09  0.05
    # 75%     0.94    0.30    0.74    0.79    0.42    0.84    0.34  0.17  0.07
    # 100%    9.41    7.01    2.50    8.65    7.82    6.07    3.02  3.62  9.43
    #      Neuron.mixed Oligo  OPC OPC_COP
    # 0%           0.04  0.00 0.00    0.05
    # 25%          0.17  0.03 0.01    0.24
    # 50%          0.33  0.07 0.04    0.66
    # 75%          0.69  0.19 0.91    0.93
    # 100%         5.91  7.20 9.30    8.48

sapply(cellClust.idx, function(x){quantile(sce.ls$sum[x])})
    #      Aqp4.Rbpms    Astro drop.doublet drop.likelyDoublet drop.lowNTx     Endo
    # 0%         1475   770.00      2524.00                995       716.0   912.00
    # 25%        5142  2780.00      5443.75               4256       959.0  2506.75
    # 50%        7382  3679.50      7017.00               7501      1279.0  5411.00
    # 75%       11051  4991.75      8601.25              12432      1872.5  9237.75
    # 100%     181764 32378.00     19299.00              56422     29767.0 48895.00

    #      Ependymal Excit_A Excit_B  Excit_C Excit_D Excit_E  Excit_F  Excit_G
    # 0%        3270     980  1377.0  3157.00  3191.0    5919  3194.00  1972.00
    # 25%       5620    4952  8945.0  9379.25  5481.5    9051 10119.75 10323.25
    # 50%       7739    6694 12002.5 13230.50  8429.0   10714 14389.50 14305.00
    # 75%      12165    9044 16862.0 18958.50 10386.5   14075 19935.75 20099.00
    # 100%     38139   43894 64847.0 58550.00 23709.0   35325 60231.00 54045.00

    #      Inhib_A  Inhib_B Inhib_C Inhib_D  Inhib_E  Inhib_F  Inhib_G  Inhib_H
    # 0%      1947   933.00  3595.0  1069.0  1707.00    850.0  1768.00  2817.00
    # 25%     5937  2680.00  9057.5  6909.5  3605.50   7688.5  6224.00  7266.25
    # 50%     8099  4112.00 12421.0  9056.0  5034.00  10395.0  8654.00  9130.50
    # 75%    11043  7178.75 17005.0 12044.5  7344.25  13950.0 11837.75 13190.75
    # 100%   53149 83174.00 33512.0 45681.0 40784.00 134792.0 51958.00 50031.00

    #      Inhib_I  Inhib_J  Inhib_K Inhib_L Inhib_M Inhib_N Inhib_O  Inhib_P
    # 0%    5982.0  2179.00  4684.00    3068    4600    3005  6673.0  5602.00
    # 25%  11439.5  7938.00  8697.00    7046    8307    6080 11536.5  8687.25
    # 50%  16729.0 10361.00 11063.50    9079   11863    7613 13859.5 10348.50
    # 75%  21484.5 13901.75 12831.25   11844   16078   10506 17912.5 13666.25
    # 100% 57135.0 44825.00 36869.00   36314   40109   21836 57038.0 26839.00

    #       Inhib_Q    Micro    Mural Neuron.mixed   Oligo      OPC OPC_COP
    # 0%    5720.00   881.00   737.00         3855   847.0  1037.00    2805
    # 25%  11226.75  2009.25  2244.75         8877  2216.5  3564.25    5800
    # 50%  13534.50  2908.50  3384.00        12019  2947.0  5289.00    7404
    # 75%  17733.25  4202.00  5403.00        16701  4100.5  7997.25   11560
    # 100% 51899.00 36079.00 36689.00        47248 52916.0 41110.00   36554




# First drop any flagged clusters for dropping    - not doing this time (until we have final annotations)
sce.ls <- sce.ls[ ,-grep("drop.",sce.ls$cellType.final)]
sce.ls <- sce.ls[ ,-grep("Neuron.mixed", sce.ls$cellType.final)]
sce.ls$cellType.final <- droplevels(sce.ls$cellType.final)

# And corresponding cell type colors
cell_colors.ls <- cell_colors.ls[-grep("drop", names(cell_colors.ls))]
cell_colors.ls <- cell_colors.ls[-grep("Neuron.mixed", names(cell_colors.ls))]


# Remove 0 genes across all nuclei
sce.ls <- sce.ls[!rowSums(assay(sce.ls, "counts"))==0, ]  #


## Re-create 'logcounts' (don't want to use 'multiBatchNorm's down-scaling across donor 'batches')
    # First 'hold' the MBN 'logcounts' for printing
sce.hold <- sce.ls

assay(sce.ls, "logcounts") <- NULL
sizeFactors(sce.ls) <- NULL
sce.ls <- logNormCounts(sce.ls)


### First make a list of Boolean param / cell subtype ===
  # Will use this to assess more 'valid', non-noise-driving markers
cellClust.idx <- splitit(sce.ls$cellType.final)
medianNon0.ls <- lapply(cellClust.idx, function(x){
  apply(as.matrix(assay(sce.ls, "logcounts")), 1, function(y){
    median(y[x]) > 0
  })
})

sapply(medianNon0.ls, table)
    #      Aqp4.Rbpms Astro  Endo Ependymal Excit_A Excit_B Excit_C Excit_D Excit_E
    # FALSE      25636 27052 26439     25184   26193   24722   24492   25803   24898
    # TRUE        2115   699  1312      2567    1558    3029    3259    1948    2853

    #       Excit_F Excit_G Inhib_A Inhib_B Inhib_C Inhib_D Inhib_E Inhib_F Inhib_G
    # FALSE   24414   24422   25734   26804   24630   25312   26549   25221   25505
    # TRUE     3337    3329    2017     947    3121    2439    1202    2530    2246

    #       Inhib_H Inhib_I Inhib_J Inhib_K Inhib_L Inhib_M Inhib_N Inhib_O Inhib_P
    # FALSE   25242   23913   25004   24636   25474   24957   25698   24134   24702
    # TRUE     2509    3838    2747    3115    2277    2794    2053    3617    3049

    #       Inhib_Q Micro Mural Oligo   OPC OPC_COP
    # FALSE   24050 27078 27209 27217 26439   25887
    # TRUE     3701   673   542   534  1312    1864


# # Just for interactive exploration of some of these
# plotExpressionCustom(sce = sce.ls,
#                      exprs_values = "logcounts",
#                      #
#                      features = head(Oligo_A_markers,4),
#                      features_name = "custom-selected",
#                      anno_name = "cellType.final",
#                      ncol=2, point_alpha=0.4, point_size=0.9,
#                      scales="free_y", swap_rownames="gene_name") +
#     ggtitle(label=paste0("Lionel's custom markers of interest")) +
#     theme(plot.title = element_text(size = 12),
#           axis.text.x = element_text(size=7))


## Traditional t-test, pairwise ===
mod <- with(colData(sce.ls), model.matrix(~ Sample))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

# Run pairwise t-tests
markers.ls.t.pw <- findMarkers(sce.ls, groups=sce.ls$cellType.final,
                               assay.type="logcounts", design=mod, test="t",
                               direction="up", pval.type="all", full.stats=T)

sapply(markers.ls.t.pw, function(x){table(x$FDR<0.05)})
    #




# Add respective 'non0median' column to the stats for each set of markers
for(i in names(markers.ls.t.pw)){
  markers.ls.t.pw[[i]] <- cbind(markers.ls.t.pw[[i]],
                                medianNon0.ls[[i]][match(rownames(markers.ls.t.pw[[i]]),
                                                         names(medianNon0.ls[[i]]))])
  colnames(markers.ls.t.pw[[i]])[36] <- "non0median"
}

sapply(markers.ls.t.pw, function(x){table(x$FDR<0.05 & x$non0median == TRUE)["TRUE"]})
    #Aqp4.Rbpms.TRUE      Astro.TRUE       Endo.TRUE  Ependymal.TRUE    Excit_A.TRUE 
    #            496              98             328             272              47 
    #   Excit_B.TRUE    Excit_C.TRUE    Excit_D.TRUE    Excit_E.TRUE    Excit_F.TRUE 
    #             26              15              69              24              24 
    #   Excit_G.TRUE    Inhib_A.TRUE    Inhib_B.TRUE    Inhib_C.TRUE      Inhib_D.NA 
    #              2               3              34              14              NA 
    #   Inhib_E.TRUE    Inhib_F.TRUE      Inhib_G.NA    Inhib_H.TRUE    Inhib_I.TRUE 
    #             41              13              NA              36               6 
    #   Inhib_J.TRUE    Inhib_K.TRUE    Inhib_L.TRUE    Inhib_M.TRUE    Inhib_N.TRUE 
    #              1              25              36              14              15 
    #   Inhib_O.TRUE    Inhib_P.TRUE    Inhib_Q.TRUE      Micro.TRUE      Mural.TRUE 
    #              2               5              21             180              17 
    #     Oligo.TRUE        OPC.TRUE    OPC_COP.TRUE 
    #            119              50             141 

## Save these
save(markers.ls.t.pw, medianNon0.ls,
     file=here("snRNAseq_mouse", "processed_data","SCE",
               "markers-stats_LS-n4_findMarkers_33cellTypes.rda"))

    # # As needed
    # load(here("snRNAseq_mouse", "processed_data","SCE",
    #           "markers-stats_LS-n4_findMarkers_33cellTypes.rda"), verbose=T)


# Print these to pngs
markerList.t.pw <- lapply(markers.ls.t.pw, function(x){
  rownames(x)[x$FDR < 0.05 & x$non0median == TRUE]
  }
)

# Change to gene symbols
markerList.t.pw <- lapply(markerList.t.pw, function(x){
  rowData(sce.ls)$gene_name[match(x, rowData(sce.ls)$gene_id)]
})

genes.top40.t <- lapply(markerList.t.pw, function(x){head(x, n=40)})

smaller.set <- names(genes.top40.t)[lengths(genes.top40.t) <= 20]
left.set <- setdiff(names(genes.top40.t), smaller.set)

# Now remove those with no significant markers
smaller.set <- setdiff(smaller.set,
                       names(genes.top40.t)[lengths(genes.top40.t) == 0])

# Smaller graphical window
dir.create(here("snRNAseq_mouse","plots","markers"))
for(i in smaller.set){
  png(here("snRNAseq_mouse","plots","markers",
           paste0("LS_t_pairwise_top40markers-", i, "_logExprs.png")), height=950, width=1200)
  print(
    plotExpressionCustom(sce = sce.hold,
                         features = genes.top40.t[[i]],
                         features_name = i,
                         anno_name = "cellType.final",
                         ncol=5, point_alpha=0.4,
                         scales="free_y", swap_rownames="gene_name") +
      scale_color_manual(values=cell_colors.ls) +
      ggtitle(label=paste0("LS ", i, " top markers: single-nucleus-level p.w. t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}

# 20-40 markers
for(i in left.set){
  png(here("snRNAseq_mouse","plots","markers",
           paste0("LS_t_pairwise_top40markers-", i, "_logExprs.png")), height=1900, width=1200)
  print(
    plotExpressionCustom(sce = sce.hold,
                         features = genes.top40.t[[i]],
                         features_name = i,
                         anno_name = "cellType.final",
                         ncol=5, point_alpha=0.4,
                         scales="free_y", swap_rownames="gene_name") +
      scale_color_manual(values=cell_colors.ls) +
      ggtitle(label=paste0("LS ", i, " top markers: single-nucleus-level p.w. t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}



### Cluster-vs-all-others single-nucleus-level iteration ========

# (Pre-process as above, if needed)
mod <- with(colData(sce.ls), model.matrix(~ Sample))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

markers.ls.t.1vAll <- list()
for(i in levels(sce.ls$cellType.final)){
  # Make temporary contrast
  sce.ls$contrast <- ifelse(sce.ls$cellType.final==i, 1, 0)
  # Test cluster vs. all others
  markers.ls.t.1vAll[[i]] <- findMarkers(sce.ls, groups=sce.ls$contrast,
                                         assay.type="logcounts", design=mod, test="t",
                                         std.lfc=TRUE,
                                         direction="up", pval.type="all", full.stats=T)
}


## Save these
save(markers.ls.t.pw, markers.ls.t.1vAll, medianNon0.ls,
     file=here("snRNAseq_mouse", "processed_data","SCE",
               "markers-stats_LS-n4_findMarkers_33cellTypes.rda"))


class(markers.ls.t.1vAll[[1]])
    # a SimpleList of length 2, named "0" and "1" (from the temporary 'contrast')
    # -> we want the second entry, named "1"
    #    (for other purposes, might be interesting to look into that "0" entry, which
    #     is basically what genes are depleted in the cell type of interest)


sapply(markers.ls.t.1vAll, function(x){
  table(x[["1"]]$stats.0$log.FDR < log(.001))
})
    #

# Do some reorganizing
markers.ls.t.1vAll <- lapply(markers.ls.t.1vAll, function(x){
  # Basically take the 'stats.[1 or 0]' since is redundant with the 'summary'-level stats
  lapply(x, function(y){ y[ ,4] })
})

# Re-name std.lfc column and the entries; add non-0-median info
for(i in names(markers.ls.t.1vAll)){
  colnames(markers.ls.t.1vAll[[i]][["0"]])[1] <- "std.logFC"
  colnames(markers.ls.t.1vAll[[i]][["1"]])[1] <- "std.logFC"
  # Add non0median Boolean - might be informative for both sets of stats
  markers.ls.t.1vAll[[i]][["0"]] <- cbind(markers.ls.t.1vAll[[i]][["0"]],
                                          medianNon0.ls[[i]][match(rownames(markers.ls.t.1vAll[[i]][["0"]]),
                                                                   names(medianNon0.ls[[i]]))])
  colnames(markers.ls.t.1vAll[[i]][["0"]])[4] <- "non0median"

  # "1" aka 'enriched'
  markers.ls.t.1vAll[[i]][["1"]] <- cbind(markers.ls.t.1vAll[[i]][["1"]],
                                          medianNon0.ls[[i]][match(rownames(markers.ls.t.1vAll[[i]][["1"]]),
                                                                   names(medianNon0.ls[[i]]))])
  colnames(markers.ls.t.1vAll[[i]][["1"]])[4] <- "non0median"

  # Then re-name the entries to more interpretable, because we'll keeping both contrasts
  names(markers.ls.t.1vAll[[i]]) <- paste0(i,c("_depleted", "_enriched"))
}



## Marker numbers with the non-0-median filter
sapply(markers.ls.t.1vAll, function(x){
  table(x[[2]]$log.FDR < log(.001) & x[[2]]$non0median == TRUE)
})
    #       Aqp4.Rbpms Astro  Endo Ependymal Excit_A Excit_B Excit_C Excit_D Excit_E
    # FALSE      26399 27270 26938     26714   27036   26120   26740   27243   27395
    # TRUE        1352   481   813      1037     715    1631    1011     508     356

    #       Excit_F Excit_G Inhib_A Inhib_B Inhib_C Inhib_D Inhib_E Inhib_F Inhib_G
    # FALSE   26935   26309   26845   27193   26882   26237   27222   26237   26350
    # TRUE      816    1442     906     558     869    1514     529    1514    1401

    #       Inhib_H Inhib_I Inhib_J Inhib_K Inhib_L Inhib_M Inhib_N Inhib_O Inhib_P
    # FALSE   27248   27028   26528   27446   26903   27451   27328   26648   27452
    # TRUE      503     723    1223     305     848     300     423    1103     299

    #       Inhib_Q Micro Mural Oligo   OPC OPC_COP
    # FALSE   26921 27338 27529 27367 27213   27300
    # TRUE      830   413   222   384   538     451


## Print these to pngs
markerList.t.1vAll <- lapply(markers.ls.t.1vAll, function(x){
  rownames(x[[2]])[ x[[2]]$log.FDR < log(0.05) & x[[2]]$non0median==TRUE ]
 }
)

# Change to gene symbols
markerList.t.1vAll <- lapply(markerList.t.1vAll, function(x){
  rowData(sce.ls)$gene_name[match(x, rowData(sce.ls)$gene_id)]
})

genes.top40.t <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t)){
  png(here("snRNAseq_mouse","plots","markers",
           paste0("LS_t_1vALL_top40markers-",i,"_logExprs.png")), height=1900, width=1200)
  print(
    plotExpressionCustom(sce = sce.hold,
                         features = genes.top40.t[[i]], 
                         features_name = i,
                         anno_name = "cellType.final",
                         ncol=5, point_alpha=0.4,
                         scales="free_y", swap_rownames="gene_name") +
      scale_color_manual(values=cell_colors.ls) +  
      ggtitle(label=paste0("LS ", i, " top markers: 'cluster-vs-all-others' t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}


## Write these top 40 lists to a csv ===
names(markerList.t.pw) <- paste0(names(markerList.t.pw),"_pw")
names(markerList.t.1vAll) <- paste0(names(markerList.t.1vAll),"_1vAll")

# Many of the PW results don't have 40 markers:
extend.idx <- names(which(lengths(markerList.t.pw) < 40))
for(i in extend.idx){
  markerList.t.pw[[i]] <- c(markerList.t.pw[[i]], rep("", 40-length(markerList.t.pw[[i]])))
}

top40genes <- cbind(sapply(markerList.t.pw, function(x) head(x, n=40)),
                    sapply(markerList.t.1vAll, function(y) head(y, n=40)))
top40genes <- top40genes[ ,sort(colnames(top40genes))]

write.csv(top40genes, file=here("snRNAseq_mouse", "processed_data","tables",
                                "top40genesLists_LS-n4_33finalCellTypes.csv"),
          row.names=FALSE)





## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
#[1] "2022-06-24 16:42:42 EDT"

proc.time()
#     user    system   elapsed 
#11692.065   150.947 11849.534
options(width = 120)
session_info()
#─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.1.2 Patched (2021-11-04 r81138)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2022-06-24
# pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
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
# cowplot                1.1.1    2020-12-30 [2] CRAN (R 4.1.2)
# crayon                 1.5.1    2022-03-26 [2] CRAN (R 4.1.2)
# DBI                    1.1.3    2022-06-18 [2] CRAN (R 4.1.2)
# DelayedArray           0.20.0   2021-10-26 [2] Bioconductor
# DelayedMatrixStats     1.16.0   2021-10-26 [2] Bioconductor
# digest                 0.6.29   2021-12-01 [2] CRAN (R 4.1.2)
# dplyr                  1.0.9    2022-04-28 [2] CRAN (R 4.1.2)
# dqrng                  0.3.0    2021-05-01 [2] CRAN (R 4.1.2)
# edgeR                  3.36.0   2021-10-26 [2] Bioconductor
# ellipsis               0.3.2    2021-04-29 [2] CRAN (R 4.1.0)
# fansi                  1.0.3    2022-03-24 [2] CRAN (R 4.1.2)
# farver                 2.1.0    2021-02-28 [2] CRAN (R 4.1.0)
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
# igraph                 1.3.2    2022-06-13 [2] CRAN (R 4.1.2)
# IRanges              * 2.28.0   2021-10-26 [2] Bioconductor
# irlba                  2.3.5    2021-12-06 [2] CRAN (R 4.1.2)
# jaffelab             * 0.99.31  2021-12-13 [1] Github (LieberInstitute/jaffelab@2cbd55a)
# labeling               0.4.2    2020-10-20 [2] CRAN (R 4.1.0)
# lattice                0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
# lifecycle              1.0.1    2021-09-24 [2] CRAN (R 4.1.2)
# limma                  3.50.3   2022-04-07 [2] Bioconductor
# locfit                 1.5-9.5  2022-03-03 [2] CRAN (R 4.1.2)
# magrittr               2.0.3    2022-03-30 [2] CRAN (R 4.1.2)
# MASS                   7.3-56   2022-03-23 [3] CRAN (R 4.1.2)
# Matrix                 1.4-1    2022-03-23 [3] CRAN (R 4.1.2)
# MatrixGenerics       * 1.6.0    2021-10-26 [2] Bioconductor
# matrixStats          * 0.62.0   2022-04-19 [2] CRAN (R 4.1.2)
# metapod                1.2.0    2021-10-26 [2] Bioconductor
# munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
# nlme                   3.1-157  2022-03-25 [3] CRAN (R 4.1.2)
# pillar                 1.7.0    2022-02-01 [2] CRAN (R 4.1.2)
# pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
# purrr                  0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# R6                     2.5.1    2021-08-19 [2] CRAN (R 4.1.2)
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.1.2)
# RColorBrewer           1.1-3    2022-04-03 [2] CRAN (R 4.1.2)
# Rcpp                   1.0.8.3  2022-03-17 [2] CRAN (R 4.1.2)
# RCurl                  1.98-1.7 2022-06-09 [2] CRAN (R 4.1.2)
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
# segmented              1.6-0    2022-05-31 [1] CRAN (R 4.1.2)
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
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────

