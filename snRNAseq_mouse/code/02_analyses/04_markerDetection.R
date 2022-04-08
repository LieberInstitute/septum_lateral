### mouse LS snRNA-seq analysis
### Feature selection with deviance & clustering
###     qrsh -l bluejay,mf=60G,h_vmem=64G,h_fsize=40G
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

# table(sce.ls$cellType)
#     # ambig.glial_A  ambig.glial_B   Aqp4.Rbpms_A   Aqp4.Rbpms_B        Astro_A 
#     #            63             76            106            146           4153 
#     #       Astro_B  Astro.OPC_COP   drop.doublet    drop.lowNTx         Endo_A 
#     #           171           1529            150            180            118 
#     #        Endo_B        Excit_A        Excit_B        Excit_C        Excit_D 
#     #            40            451            670            232             67 
#     #       Excit_E        Excit_F        Excit_G        Inhib_A        Inhib_B 
#     #            34            202            728            797            760 
#     #       Inhib_C        Inhib_D        Inhib_E        Inhib_F        Inhib_G 
#     #           191           4327            306           4364            363 
#     #       Inhib_H          Micro        Mural_A        Mural_B Neuron.mixed_A 
#     #            86            166            112             34             73 
#     #Neuron.Ppp1r1b        Oligo_A        Oligo_B            OPC        OPC_COP 
#     #            65           1478            105            464             53 
# 
# ## doubletScore distributions / cluster?
# cellClust.idx <- splitit(sce.ls$cellType)
# sapply(cellClust.idx, function(x){round(quantile(sce.ls$sum[x]), 2)})
# 
# 
# sapply(cellSubtype.idx, function(x){quantile(sce.ls$sum[x])})
#     #      ambig.glial_A ambig.glial_B Aqp4.Rbpms_A Aqp4.Rbpms_B Astro_A Astro_B
#     # 0%          9117.0        784.00      5539.00         3802     875    1475
#     # 25%        14281.5       1690.75      7880.50         5266    2820    3790
#     # 50%        17203.0       2100.00      9336.50         6232    3752    5009
#     # 75%        26055.0       2555.50     10952.25         7650    5090    8110
#     # 100%      181764.0      10268.00     57022.00        32300   32378   32256
# 
#     #      Astro.OPC_COP drop.doublet drop.lowNTx   Endo_A   Endo_B Excit_A Excit_B
#     # 0%             770      2524.00       716.0   912.00  3874.00   980.0  1377.0
#     # 25%           2692      5443.75       959.0  2279.00  8155.75  4883.5  8945.0
#     # 50%           3506      7017.00      1279.0  3459.00 11137.50  6601.0 12002.5
#     # 75%           4605      8601.25      1872.5  5996.25 20012.00  9014.0 16862.0
#     # 100%         28529     19299.00     29767.0 17195.00 48895.00 43894.0 64847.0
# 
#     #       Excit_C Excit_D Excit_E  Excit_F  Excit_G Inhib_A  Inhib_B Inhib_C
#     # 0%    3157.00  3191.0    3851  3194.00  1972.00    1947   933.00  3595.0
#     # 25%   9379.25  5481.5    5900 10119.75 10323.25    5937  2680.00  9057.5
#     # 50%  13230.50  8429.0    7018 14389.50 14305.00    8099  4112.00 12421.0
#     # 75%  18958.50 10386.5    9027 19935.75 20099.00   11043  7178.75 17005.0
#     # 100% 58550.00 23709.0   32275 60231.00 54045.00   53149 83174.00 33512.0
# 
#     #      Inhib_D  Inhib_E   Inhib_F Inhib_G  Inhib_H    Micro  Mural_A  Mural_B
#     # 0%     995.0  1707.00    850.00  2053.0  2817.00   881.00   737.00  3047.00
#     # 25%   6726.5  3605.50   7810.75  6602.0  7266.25  2009.25  2773.00  5336.25
#     # 50%   9336.0  5034.00  10514.50  8630.0  9130.50  2908.50  3796.00  6998.00
#     # 75%  12936.5  7344.25  13983.25 13227.5 13190.75  4202.00  5435.75 11711.75
#     # 100% 57135.0 40784.00 134792.00 47249.0 50031.00 36079.00 28953.00 36689.00
# 
#     #      Neuron.mixed_A Neuron.Ppp1r1b Oligo_A Oligo_B      OPC OPC_COP
#     # 0%             3855           3270   847.0    5196  1037.00    2805
#     # 25%            8877           5620  2183.5    6891  3564.25    5800
#     # 50%           12019           7739  2843.5    8311  5289.00    7404
#     # 75%           16701          12165  3784.0   11659  7997.25   11560
#     # 100%          47248          38139 32078.0   52916 41110.00   36554


# # First drop any flagged clusters for dropping    - not doing this time (until we have final annotations)
# sce.ls <- sce.ls[ ,-grep("drop.",sce.ls$cellType)]
# sce.ls$cellType <- droplevels(sce.ls$cellType)

# Remove 0 genes across all nuclei
sce.ls <- sce.ls[!rowSums(assay(sce.ls, "counts"))==0, ]  #


## Re-create 'logcounts' (don't want to use 'multiBatchNorm's down-scaling across donor 'batches')
    # First 'hold' the MBN 'logcounts' for printing
sce.hold <- sce.ls

assay(sce.ls, "logcounts") <- NULL
sizeFactors(sce.ls) <- NULL
sce.ls <- logNormCounts(sce.ls)


# ### First make a list of Boolean param / cell subtype ===
#   # Will use this to assess more 'valid', non-noise-driving markers
# cellSubtype.idx <- splitit(sce.ls$cellType)
# medianNon0.ls <- lapply(cellSubtype.idx, function(x){
#   apply(as.matrix(assay(sce.ls, "logcounts")), 1, function(y){
#     median(y[x]) > 0
#   })
# })
# 
# sapply(medianNon0.ls, table)
#     #      ambig.glial_A ambig.glial_B Aqp4.Rbpms_A Aqp4.Rbpms_B Astro_A Astro_B
#     # FALSE         22941         27523        25205        26064   27084   26362
#     # TRUE           4856           274         2592         1733     713    1435
#     #       Astro.OPC_COP drop.doublet drop.lowNTx Endo_A Endo_B Excit_A Excit_B
#     # FALSE         27132        26062       27619  26978  24033   26253   24768
#     # TRUE            665         1735         178    819   3764    1544    3029
#     #       Excit_C Excit_D Excit_E Excit_F Excit_G Inhib_A Inhib_B Inhib_C Inhib_D
#     # FALSE   24538   25849   25958   24460   24468   25780   26850   24676   25403
#     # TRUE     3259    1948    1839    3337    3329    2017     947    3121    2394
#     #       Inhib_E Inhib_F Inhib_G Inhib_H Micro Mural_A Mural_B Neuron.mixed_A
#     # FALSE   26595   25263   25526   25288 27124   27012   25577          24928
#     # TRUE     1202    2534    2271    2509   673     785    2220           2869
#     #       Neuron.Ppp1r1b Oligo_A Oligo_B   OPC OPC_COP
#     # FALSE          25230   27311   25609 26485   25933
#     # TRUE            2567     486    2188  1312    1864
# 
# 
# # Just for interactive exploration of some of these
# plotExpressionCustom(sce = sce.ls,
#                      exprs_values = "logcounts",
#                      # 
#                      features = head(Oligo_A_markers,4),
#                      features_name = "custom-selected",
#                      anno_name = "cellType",
#                      ncol=2, point_alpha=0.4, point_size=0.9,
#                      scales="free_y", swap_rownames="gene_name") +
#     ggtitle(label=paste0("Lionel's custom markers of interest")) +
#     theme(plot.title = element_text(size = 12),
#           axis.text.x = element_text(size=7))
# 
# 
# ## Traditional t-test, pairwise ===
# mod <- with(colData(sce.ls), model.matrix(~ Sample))
# mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`
# 
# # Run pairwise t-tests
# markers.ls.t.pw <- findMarkers(sce.ls, groups=sce.ls$cellType,
#                                assay.type="logcounts", design=mod, test="t",
#                                direction="up", pval.type="all", full.stats=T)
# 
# sapply(markers.ls.t.pw, function(x){table(x$FDR<0.05)})
#     #      ambig.glial_A ambig.glial_B Aqp4.Rbpms_A Aqp4.Rbpms_B Astro_A Astro_B
#     # FALSE         27511         27687        27697        27746   27745   27758
#     # TRUE            286           110          100           51      52      39
#     #       Astro.OPC_COP drop.doublet drop.lowNTx Endo_A Endo_B Excit_A Excit_B
#     # FALSE         27713        27789       27729  27758  27006   27760   27716
#     # TRUE             84            8          68     39    791      37      81
#     #       Excit_C Excit_D Excit_E Excit_F Excit_G Inhib_A Inhib_B Inhib_C Inhib_D
#     # FALSE   27691   27572   27690   27705   27790   27780   27626   27708   27796
#     # TRUE      106     225     107      92       7      17     171      89       1
#     #       Inhib_E Inhib_F Inhib_G Inhib_H Micro Mural_A Mural_B Neuron.mixed_A
#     # FALSE   27715   27753   27729   27731 27146   27584   27460          27761
#     # TRUE       82      44      68      66   651     213     337             36
#     #       Neuron.Ppp1r1b Oligo_A Oligo_B   OPC OPC_COP
#     # FALSE          27300   27783   27613 27714   27587
#     # TRUE             497      14     184    83     210
# 
# 
# 
# 
# # Add respective 'non0median' column to the stats for each set of markers
# for(i in names(markers.ls.t.pw)){
#   markers.ls.t.pw[[i]] <- cbind(markers.ls.t.pw[[i]],
#                                 medianNon0.ls[[i]][match(rownames(markers.ls.t.pw[[i]]),
#                                                          names(medianNon0.ls[[i]]))])
#   colnames(markers.ls.t.pw[[i]])[38] <- "non0median"
# }
# 
# sapply(markers.ls.t.pw, function(x){table(x$FDR<0.05 & x$non0median == TRUE)["TRUE"]})
#     # ambig.glial_A.TRUE  ambig.glial_B.TRUE   Aqp4.Rbpms_A.TRUE     Aqp4.Rbpms_B.NA 
#     #                142                  22                   8                  NA 
#     # Astro_A.TRUE          Astro_B.NA  Astro.OPC_COP.TRUE     drop.doublet.NA 
#     #           30                  NA                  49                  NA 
#     # drop.lowNTx.TRUE           Endo_A.NA         Endo_B.TRUE        Excit_A.TRUE 
#     #           33                  NA                 511                  18 
#     # Excit_B.TRUE        Excit_C.TRUE        Excit_D.TRUE        Excit_E.TRUE 
#     #           36                  44                  74                  32 
#     # Excit_F.TRUE        Excit_G.TRUE        Inhib_A.TRUE        Inhib_B.TRUE 
#     #           36                   2                   6                  35 
#     # Inhib_C.TRUE        Inhib_D.TRUE        Inhib_E.TRUE        Inhib_F.TRUE 
#     #           31                   1                  44                  31 
#     # Inhib_G.TRUE        Inhib_H.TRUE          Micro.TRUE        Mural_A.TRUE 
#     #           42                  40                 148                  22 
#     # Mural_B.TRUE Neuron.mixed_A.TRUE Neuron.Ppp1r1b.TRUE        Oligo_A.TRUE 
#     #          136                   2                 195                   9 
#     # Oligo_B.TRUE            OPC.TRUE        OPC_COP.TRUE 
#     #          121                  45                 115 
# 
# ## Save these
# save(markers.ls.t.pw, medianNon0.ls,
#      file=here("snRNAseq_mouse", "processed_data","SCE",
#                "markers-stats_LS-n4_findMarkers_35cellTypes.rda"))

    # As needed
    load(here("snRNAseq_mouse", "processed_data","SCE",
              "markers-stats_LS-n4_findMarkers_35cellTypes.rda"), verbose=T)


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
#dir.create(here("snRNAseq_mouse","plots","markers"))
for(i in smaller.set){
  png(here("snRNAseq_mouse","plots","markers",
           paste0("LS_t_pairwise_top40markers-", i, "_logExprs.png")), height=950, width=1200)
  print(
    plotExpressionCustom(sce = sce.hold,
                         features = genes.top40.t[[i]], 
                         features_name = i,
                         anno_name = "cellType",
                         ncol=5, point_alpha=0.4) +
      scale_color_manual(values = cell_colors.ls) +  
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
                         anno_name = "cellType",
                         ncol=5, point_alpha=0.4) +
      scale_color_manual(values = cell_colors.ls) +  
      ggtitle(label=paste0("LS ", i, " top markers: single-nucleus-level p.w. t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}



### Cluster-vs-all-others single-nucleus-level iteration ========

# # (Pre-process as above, if needed)
# mod <- with(colData(sce.ls), model.matrix(~ Sample))
# mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`
# 
# markers.ls.t.1vAll <- list()
# for(i in levels(sce.ls$cellType)){
#   # Make temporary contrast
#   sce.ls$contrast <- ifelse(sce.ls$cellType==i, 1, 0)
#   # Test cluster vs. all others
#   markers.ls.t.1vAll[[i]] <- findMarkers(sce.ls, groups=sce.ls$contrast,
#                                          assay.type="logcounts", design=mod, test="t",
#                                          std.lfc=TRUE,
#                                          direction="up", pval.type="all", full.stats=T)
# }
# 
# 
# ## Save these
# save(markers.ls.t.pw, markers.ls.t.1vAll, medianNon0.ls,
#      file=here("snRNAseq_mouse", "processed_data","SCE",
#                "markers-stats_LS-n4_findMarkers_35cellTypes.rda"))
# 
# 
# class(markers.ls.t.1vAll[[1]])
#     # a SimpleList of length 2, named "0" and "1" (from the temporary 'contrast')
#     # -> we want the second entry, named "1"
#     #    (for other purposes, might be interesting to look into that "0" entry, which
#     #     is basically what genes are depleted in the cell type of interest)
# 
# 
# sapply(markers.ls.t.1vAll, function(x){
#   table(x[["1"]]$stats.0$log.FDR < log(.001))
# })
#     #      ambig.glial_A ambig.glial_B Aqp4.Rbpms_A Aqp4.Rbpms_B Astro_A Astro_B
#     # FALSE         24628         26869        24484        24607   23731   24970
#     # TRUE           3169           928         3313         3190    4066    2827
#     #       Astro.OPC_COP drop.doublet drop.lowNTx Endo_A Endo_B Excit_A Excit_B
#     # FALSE         25358        26213       26725  25611  25350   25947   23953
#     # TRUE           2439         1584        1072   2186   2447    1850    3844
#     #       Excit_C Excit_D Excit_E Excit_F Excit_G Inhib_A Inhib_B Inhib_C Inhib_D
#     # FALSE   25948   26705   27314   26445   24746   25739   24724   26264   20760
#     # TRUE     1849    1092     483    1352    3051    2058    3073    1533    7037
#     #       Inhib_E Inhib_F Inhib_G Inhib_H Micro Mural_A Mural_B Neuron.mixed_A
#     # FALSE   26556   23649   26609   27084 25406   26004   26488          26904
#     # TRUE     1241    4148    1188     713  2391    1793    1309            893
#     #       Neuron.Ppp1r1b Oligo_A Oligo_B   OPC OPC_COP
#     # FALSE          25365   25805   26341 26083   26941
#     # TRUE            2432    1992    1456  1714     856
# 
# # Do some reorganizing
# markers.ls.t.1vAll <- lapply(markers.ls.t.1vAll, function(x){
#   # Basically take the 'stats.[1 or 0]' since is redundant with the 'summary'-level stats
#   lapply(x, function(y){ y[ ,4] }) 
# })
# 
# # Re-name std.lfc column and the entries; add non-0-median info
# for(i in names(markers.ls.t.1vAll)){
#   colnames(markers.ls.t.1vAll[[i]][["0"]])[1] <- "std.logFC"
#   colnames(markers.ls.t.1vAll[[i]][["1"]])[1] <- "std.logFC"
#   # Add non0median Boolean - might be informative for both sets of stats
#   markers.ls.t.1vAll[[i]][["0"]] <- cbind(markers.ls.t.1vAll[[i]][["0"]],
#                                           medianNon0.ls[[i]][match(rownames(markers.ls.t.1vAll[[i]][["0"]]),
#                                                                    names(medianNon0.ls[[i]]))])
#   colnames(markers.ls.t.1vAll[[i]][["0"]])[4] <- "non0median"
#   
#   # "1" aka 'enriched'
#   markers.ls.t.1vAll[[i]][["1"]] <- cbind(markers.ls.t.1vAll[[i]][["1"]],
#                                           medianNon0.ls[[i]][match(rownames(markers.ls.t.1vAll[[i]][["1"]]),
#                                                                    names(medianNon0.ls[[i]]))])
#   colnames(markers.ls.t.1vAll[[i]][["1"]])[4] <- "non0median"
#   
#   # Then re-name the entries to more interpretable, because we'll keeping both contrasts
#   names(markers.ls.t.1vAll[[i]]) <- paste0(i,c("_depleted", "_enriched"))
# }
# 
# 
# 
# ## Marker numbers with the non-0-median filter
# sapply(markers.ls.t.1vAll, function(x){
#   table(x[[2]]$log.FDR < log(.001) & x[[2]]$non0median == TRUE)
# })
#     #       ambig.glial_A ambig.glial_B Aqp4.Rbpms_A Aqp4.Rbpms_B Astro_A Astro_B
#     # FALSE         25817         27636        26385        26767   27312   27063
#     # TRUE           1980           161         1412         1030     485     734
# 
#     #       Astro.OPC_COP drop.doublet drop.lowNTx Endo_A Endo_B Excit_A Excit_B
#     # FALSE         27343        27109       27688  27288  26278   27088   26166
#     # TRUE            454          688         109    509   1519     709    1631
# 
#     #       Excit_C Excit_D Excit_E Excit_F Excit_G Inhib_A Inhib_B Inhib_C Inhib_D
#     # FALSE   26776   27284   27543   26979   26348   26889   27243   26921   26077
#     # TRUE     1021     513     254     818    1449     908     554     876    1720
# 
#     #       Inhib_E Inhib_F Inhib_G Inhib_H Micro Mural_A Mural_B Neuron.mixed_A
#     # FALSE   27264   26268   27157   27289 27387   27505   27179          27294
#     # TRUE      533    1529     640     508   410     292     618            503
# 
#     #       Neuron.Ppp1r1b Oligo_A Oligo_B   OPC OPC_COP
#     # FALSE          26756   27434   26998 27259   27350
#     # TRUE            1041     363     799   538     447


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
  png(here("snRNA-seq","plots","markers",
           paste0("LS_t_1vALL_top40markers-",i,"_logExprs.png")), height=1900, width=1200)
  print(
    plotExpressionCustom(sce = sce.hold,
                         features = genes.top40.t[[i]], 
                         features_name = i,
                         anno_name = "cellType",
                         ncol=5, point_alpha=0.4) +
      scale_color_manual(values = cell_colors.ls) +  
      ggtitle(label=paste0("LS ", i, " top markers: 'cluster-vs-all-others' t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}


# ## Write these top 40 lists to a csv ===
# names(markerList.t.pw) <- paste0(names(markerList.t.pw),"_pw")
# names(markerList.t.1vAll) <- paste0(names(markerList.t.1vAll),"_1vAll")
# 
# # Many of the PW results don't have 40 markers:
# extend.idx <- names(which(lengths(markerList.t.pw) < 40))
# for(i in extend.idx){
#   markerList.t.pw[[i]] <- c(markerList.t.pw[[i]], rep("", 40-length(markerList.t.pw[[i]])))
# }
# 
# top40genes <- cbind(sapply(markerList.t.pw, function(x) head(x, n=40)),
#                     sapply(markerList.t.1vAll, function(y) head(y, n=40)))
# top40genes <- top40genes[ ,sort(colnames(top40genes))]
# #dir.create(here("snRNAseq_mouse", "processed_data","tables"))
# write.csv(top40genes, file=here("snRNAseq_mouse", "processed_data","tables",
#                                 "top40genesLists_LS-n4_35cellTypes.csv"),
#           row.names=FALSE)



