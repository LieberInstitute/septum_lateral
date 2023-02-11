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
library(pheatmap)

here()

#If you need to mamually set working directory
setwd("/dcs04/lieber/marmaypag/pilotLS_LIBD1070/")
getwd() #to check it's in the right place


# plotExpressionCustom for nicer aesthetics
source("/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/analyses_sn/plotExpressionCustom.R")


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

## ===


load("/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/processed_data/SCE/sce_updated_LS.rda", verbose = T)
# sce.ls, annotationTab.ls, cell_colors.ls

table(sce.ls$cellType.final)
#             Astro          Chol_Ex.D                ChP       drop.doublet
#              5682                 67                 65                150
# drop.likelyDoublet        drop.lowNTx               Endo          Ependymal
#               573                180                158                486
#          IoC_In.E            LS_In.C            LS_In.D            LS_In.M
#               306                191               1075                167
#           LS_In.N            LS_In.O            LS_In.P            LS_In.Q
#                45                 69                190                 36
#           LS_In.R              Micro            MS_In.J            MS_In.K
#                68                166                111                402
#             Mural       Neuron.mixed              Oligo                OPC
#               222                 73               1583                464
#           OPC_COP          Sept_In.G          Sept_In.I           Str_In.A
#                53                363               1524                797
#          Str_In.F           Str_In.H           Str_In.L          Thal_Ex.B
#              4364                 86                 36                670
#         TNoS_Ex.A      TT.IG.SH_Ex.C      TT.IG.SH_Ex.E      TT.IG.SH_Ex.F
#               485                232                 31                930
#        Ventr_In.B
#               760

## doubletScore distributions / cluster?
cellClust.idx <- splitit(sce.ls$cellType.final)
sapply(cellClust.idx, function(x) {
  round(quantile(sce.ls$doubletScore[x]), 2)
})
#      Astro Chol_Ex.D  ChP drop.doublet drop.likelyDoublet drop.lowNTx Endo
# 0%    0.00      0.04 0.02         0.06               0.00        0.00 0.00
# 25%   0.06      0.38 0.36         3.78               0.62        0.01 0.04
# 50%   0.18      1.28 0.49         4.83               2.01        0.03 0.06
# 75%   0.48      1.39 0.53         6.72               2.94        0.06 0.17
# 100% 14.84      9.93 7.10         9.39              11.38        7.93 5.98

#      Ependymal IoC_In.E LS_In.C LS_In.D LS_In.M LS_In.N LS_In.O LS_In.P LS_In.Q
# 0%        0.00     0.00    0.01    0.01    0.01    0.06    0.13    0.02    0.04
# 25%       0.04     0.02    0.23    0.17    0.11    0.32    0.28    0.17    0.19
# 50%       0.10     0.06    0.43    0.36    0.18    0.42    0.40    0.25    0.34
# 75%       0.46     0.52    1.08    0.93    0.30    0.74    0.79    0.42    0.84
# 100%     10.39     6.60    5.98   11.49    7.01    2.50    8.65    7.82    6.07

#      LS_In.R Micro MS_In.J MS_In.K Mural Neuron.mixed Oligo  OPC OPC_COP
# 0%      0.07  0.00    0.00    0.00  0.00         0.04  0.00 0.00    0.05
# 25%     0.16  0.06    0.12    0.13  0.04         0.17  0.03 0.01    0.24
# 50%     0.23  0.09    0.23    0.28  0.05         0.33  0.07 0.04    0.66
# 75%     0.34  0.17    0.87    0.67  0.07         0.69  0.19 0.91    0.93
# 100%    3.02  3.62    5.92   11.43  9.43         5.91  7.20 9.30    8.48

#      Sept_In.G Sept_In.I Str_In.A Str_In.F Str_In.H Str_In.L Thal_Ex.B
# 0%        0.00      0.00     0.01     0.00     0.09     0.01      0.02
# 25%       0.04      0.13     0.13     0.16     0.30     0.26      0.15
# 50%       0.11      0.27     0.25     0.29     0.52     0.39      0.39
# 75%       0.84      0.74     0.62     0.90     0.91     0.94      1.31
# 100%      9.50     11.42    10.44    10.09     8.40     9.41     10.41

#      TNoS_Ex.A TT.IG.SH_Ex.C TT.IG.SH_Ex.E TT.IG.SH_Ex.F Ventr_In.B
# 0%        0.00          0.03          0.12          0.02       0.00
# 25%       0.04          0.14          0.61          0.16       0.03
# 50%       0.07          0.46          0.78          0.29       0.07
# 75%       0.27          1.30          1.21          0.85       0.78
# 100%     12.53          6.82          7.56          9.25      11.46

sapply(cellClust.idx, function(x) {
  quantile(sce.ls$sum[x])
})
#         Astro Chol_Ex.D   ChP drop.doublet drop.likelyDoublet drop.lowNTx
# 0%     770.00    3191.0  3270      2524.00                995       716.0
# 25%   2780.00    5481.5  5620      5443.75               4256       959.0
# 50%   3679.50    8429.0  7739      7017.00               7501      1279.0
# 75%   4991.75   10386.5 12165      8601.25              12432      1872.5
# 100% 32378.00   23709.0 38139     19299.00              56422     29767.0

#          Endo Ependymal IoC_In.E LS_In.C LS_In.D LS_In.M LS_In.N LS_In.O
# 0%     912.00      1475  1707.00  3595.0  1069.0    3068    4600    3005
# 25%   2506.75      5142  3605.50  9057.5  6909.5    7046    8307    6080
# 50%   5411.00      7382  5034.00 12421.0  9056.0    9079   11863    7613
# 75%   9237.75     11051  7344.25 17005.0 12044.5   11844   16078   10506
# 100% 48895.00    181764 40784.00 33512.0 45681.0   36314   40109   21836

#      LS_In.P  LS_In.Q  LS_In.R    Micro MS_In.J  MS_In.K    Mural Neuron.mixed
# 0%    6673.0  5602.00  5720.00   881.00  5982.0  2179.00   737.00         3855
# 25%  11536.5  8687.25 11226.75  2009.25 11439.5  7938.00  2244.75         8877
# 50%  13859.5 10348.50 13534.50  2908.50 16729.0 10361.00  3384.00        12019
# 75%  17912.5 13666.25 17733.25  4202.00 21484.5 13901.75  5403.00        16701
# 100% 57038.0 26839.00 51899.00 36079.00 57135.0 44825.00 36689.00        47248

#        Oligo      OPC OPC_COP Sept_In.G Sept_In.I Str_In.A  Str_In.F Str_In.H
# 0%     847.0  1037.00    2805    2053.0   1768.00     1947    850.00  2817.00
# 25%   2216.5  3564.25    5800    6602.0   6224.00     5937   7810.75  7266.25
# 50%   2947.0  5289.00    7404    8630.0   8654.00     8099  10514.50  9130.50
# 75%   4100.5  7997.25   11560   13227.5  11837.75    11043  13983.25 13190.75
# 100% 52916.0 41110.00   36554   47249.0  51958.00    53149 134792.00 50031.00

#      Str_In.L Thal_Ex.B TNoS_Ex.A TT.IG.SH_Ex.C TT.IG.SH_Ex.E TT.IG.SH_Ex.F
# 0%    4684.00    1377.0       980       3157.00          5919        1972.0
# 25%   8697.00    8945.0      4952       9379.25          9051       10227.0
# 50%  11063.50   12002.5      6694      13230.50         10714       14337.0
# 75%  12831.25   16862.0      9044      18958.50         14075       20052.5
# 100% 36869.00   64847.0     43894      58550.00         35325       60231.0

#      Ventr_In.B
# 0%       933.00
# 25%     2680.00
# 50%     4112.00
# 75%     7178.75
# 100%   83174.00

# First drop any flagged clusters for dropping
sce.ls <- sce.ls[, -grep("drop.", sce.ls$cellType.final)]
sce.ls <- sce.ls[, -grep("Neuron.mixed", sce.ls$cellType.final)]
sce.ls$cellType.final <- droplevels(sce.ls$cellType.final)

# And corresponding cell type colors
cell_colors.ls <- cell_colors.ls[-grep("drop", names(cell_colors.ls))]
cell_colors.ls <- cell_colors.ls[-grep("Neuron.mixed", names(cell_colors.ls))]


# Remove 0 genes across all nuclei
sce.ls <- sce.ls[!rowSums(assay(sce.ls, "counts")) == 0, ] #


## Re-create 'logcounts' (don't want to use 'multiBatchNorm's down-scaling across donor 'batches')
# First 'hold' the MBN 'logcounts' for printing
sce.hold <- sce.ls

assay(sce.ls, "logcounts") <- NULL
sizeFactors(sce.ls) <- NULL
sce.ls <- logNormCounts(sce.ls)


### First make a list of Boolean param / cell subtype ===
# Will use this to assess more 'valid', non-noise-driving markers
cellClust.idx <- splitit(sce.ls$cellType.final)
medianNon0.ls <- lapply(cellClust.idx, function(x) {
  apply(as.matrix(assay(sce.ls, "logcounts")), 1, function(y) {
    median(y[x]) > 0
  })
})

sapply(medianNon0.ls, table)
#


## Traditional t-test, pairwise ===
mod <- with(colData(sce.ls), model.matrix(~Sample))
mod <- mod[, -1, drop = F] # intercept otherwise automatically dropped by `findMarkers()`

# Run pairwise t-tests
markers.ls.t.pw <- findMarkers(sce.ls,
                               groups = sce.ls$cellType.final,
                               assay.type = "logcounts", design = mod, test = "t",
                               direction = "up", pval.type = "all", full.stats = T
)

sapply(markers.ls.t.pw, function(x) {
  table(x$FDR < 0.05)
})
#




# Add respective 'non0median' column to the stats for each set of markers
for (i in names(markers.ls.t.pw)) {
  markers.ls.t.pw[[i]] <- cbind(
    markers.ls.t.pw[[i]],
    medianNon0.ls[[i]][match(
      rownames(markers.ls.t.pw[[i]]),
      names(medianNon0.ls[[i]])
    )]
  )
  colnames(markers.ls.t.pw[[i]])[36] <- "non0median"
}

sapply(markers.ls.t.pw, function(x) {
  table(x$FDR < 0.05 & x$non0median == TRUE)["TRUE"]
})
#     Astro.TRUE     Chol_Ex.D.TRUE           ChP.TRUE          Endo.TRUE
#             98                 69                270                328
# Ependymal.TRUE      IoC_In.E.TRUE       LS_In.C.TRUE         LS_In.D.NA
#            493                 43                 14                 NA
#   LS_In.M.TRUE       LS_In.N.TRUE       LS_In.O.TRUE       LS_In.P.TRUE
#             36                 14                 17                  2
#   LS_In.Q.TRUE       LS_In.R.TRUE         Micro.TRUE       MS_In.J.TRUE
#              4                 20                180                  6
#   MS_In.K.TRUE         Mural.TRUE         Oligo.TRUE           OPC.TRUE
#              1                 17                118                 50
#   OPC_COP.TRUE     Sept_In.G.TRUE       Sept_In.I.NA      Str_In.A.TRUE
#            141                 41                 NA                  2
#  Str_In.F.TRUE      Str_In.H.TRUE      Str_In.L.TRUE     Thal_Ex.B.TRUE
#             26                 31                 24                 25
# TNoS_Ex.A.TRUE TT.IG.SH_Ex.C.TRUE TT.IG.SH_Ex.E.TRUE TT.IG.SH_Ex.F.TRUE
#             46                 16                 25                 14
# Ventr_In.B.TRUE
#             29

## Save these
save(markers.ls.t.pw, medianNon0.ls,
     file = here(
       "snRNAseq_mouse", "processed_data", "SCE",
       "markers-stats_LS-n4_findMarkers_33cellTypes.rda"
     )
)

# # As needed
# load(here("snRNAseq_mouse", "processed_data","SCE",
#           "markers-stats_LS-n4_findMarkers_33cellTypes.rda"), verbose=T)


# Print these to pngs
markerList.t.pw <- lapply(markers.ls.t.pw, function(x) {
  rownames(x)[x$FDR < 0.05 & x$non0median == TRUE]
})

# Change to gene symbols
markerList.t.pw <- lapply(markerList.t.pw, function(x) {
  rowData(sce.ls)$gene_name[match(x, rowData(sce.ls)$gene_id)]
})

genes.top40.t <- lapply(markerList.t.pw, function(x) {
  head(x, n = 40)
})

smaller.set <- names(genes.top40.t)[lengths(genes.top40.t) <= 20]
left.set <- setdiff(names(genes.top40.t), smaller.set)

# Now remove those with no significant markers
smaller.set <- setdiff(
  smaller.set,
  names(genes.top40.t)[lengths(genes.top40.t) == 0]
)

# Smaller graphical window
dir.create(here("snRNAseq_mouse", "plots", "markers"))
for (i in smaller.set) {
  png(here(
    "snRNAseq_mouse", "plots", "markers",
    paste0("LS_t_pairwise_top40markers-", i, "_logExprs.png")
  ), height = 950, width = 1200)
  print(
    plotExpressionCustom(
      sce = sce.hold,
      features = genes.top40.t[[i]],
      features_name = i,
      anno_name = "cellType.final",
      ncol = 5, point_alpha = 0.4,
      scales = "free_y", swap_rownames = "gene_name"
    ) +
      scale_color_manual(values = cell_colors.ls) +
      ggtitle(label = paste0("LS ", i, " top markers: single-nucleus-level p.w. t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}

# 20-40 markers
for (i in left.set) {
  png(here(
    "snRNAseq_mouse", "plots", "markers",
    paste0("LS_t_pairwise_top40markers-", i, "_logExprs.png")
  ), height = 1900, width = 1200)
  print(
    plotExpressionCustom(
      sce = sce.hold,
      features = genes.top40.t[[i]],
      features_name = i,
      anno_name = "cellType.final",
      ncol = 5, point_alpha = 0.4,
      scales = "free_y", swap_rownames = "gene_name"
    ) +
      scale_color_manual(values = cell_colors.ls) +
      ggtitle(label = paste0("LS ", i, " top markers: single-nucleus-level p.w. t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}



### Cluster-vs-all-others single-nucleus-level iteration ========

# (Pre-process as above, if needed)
mod <- with(colData(sce.ls), model.matrix(~Sample))
mod <- mod[, -1, drop = F] # intercept otherwise automatically dropped by `findMarkers()`

markers.ls.t.1vAll <- list()
for (i in levels(sce.ls$cellType.final)) {
  # Make temporary contrast
  sce.ls$contrast <- ifelse(sce.ls$cellType.final == i, 1, 0)
  # Test cluster vs. all others
  markers.ls.t.1vAll[[i]] <- findMarkers(sce.ls,
                                         groups = sce.ls$contrast,
                                         assay.type = "logcounts", design = mod, test = "t",
                                         std.lfc = TRUE,
                                         direction = "up", pval.type = "all", full.stats = T
  )
}


## Save these
save(markers.ls.t.pw, markers.ls.t.1vAll, medianNon0.ls,
     file = here(
       "snRNAseq_mouse", "processed_data", "SCE",
       "markers-stats_LS-n4_findMarkers_33cellTypes.rda"
     )
)


class(markers.ls.t.1vAll[[1]])
# a SimpleList of length 2, named "0" and "1" (from the temporary 'contrast')
# -> we want the second entry, named "1"
#    (for other purposes, might be interesting to look into that "0" entry, which
#     is basically what genes are depleted in the cell type of interest)


sapply(markers.ls.t.1vAll, function(x) {
  table(x[["1"]]$stats.0$log.FDR < log(.001))
})
#

# Do some reorganizing
markers.ls.t.1vAll <- lapply(markers.ls.t.1vAll, function(x) {
  # Basically take the 'stats.[1 or 0]' since is redundant with the 'summary'-level stats
  lapply(x, function(y) {
    y[, 4]
  })
})

# Re-name std.lfc column and the entries; add non-0-median info
for (i in names(markers.ls.t.1vAll)) {
  colnames(markers.ls.t.1vAll[[i]][["0"]])[1] <- "std.logFC"
  colnames(markers.ls.t.1vAll[[i]][["1"]])[1] <- "std.logFC"
  # Add non0median Boolean - might be informative for both sets of stats
  markers.ls.t.1vAll[[i]][["0"]] <- cbind(
    markers.ls.t.1vAll[[i]][["0"]],
    medianNon0.ls[[i]][match(
      rownames(markers.ls.t.1vAll[[i]][["0"]]),
      names(medianNon0.ls[[i]])
    )]
  )
  colnames(markers.ls.t.1vAll[[i]][["0"]])[4] <- "non0median"
  
  # "1" aka 'enriched'
  markers.ls.t.1vAll[[i]][["1"]] <- cbind(
    markers.ls.t.1vAll[[i]][["1"]],
    medianNon0.ls[[i]][match(
      rownames(markers.ls.t.1vAll[[i]][["1"]]),
      names(medianNon0.ls[[i]])
    )]
  )
  colnames(markers.ls.t.1vAll[[i]][["1"]])[4] <- "non0median"
  
  # Then re-name the entries to more interpretable, because we'll keeping both contrasts
  names(markers.ls.t.1vAll[[i]]) <- paste0(i, c("_depleted", "_enriched"))
}



## Marker numbers with the non-0-median filter
sapply(markers.ls.t.1vAll, function(x) {
  table(x[[2]]$log.FDR < log(.001) & x[[2]]$non0median == TRUE)
})
#       Astro Chol_Ex.D   ChP  Endo Ependymal IoC_In.E LS_In.C LS_In.D LS_In.M
# FALSE 27270     27243 26714 26938     26399    27222   26882   26237   26903
# TRUE    481       508  1037   813      1352      529     869    1514     848

#       LS_In.N LS_In.O LS_In.P LS_In.Q LS_In.R Micro MS_In.J MS_In.K Mural Oligo
# FALSE   27451   27328   26648   27452   26921 27338   27028   26528 27529 27367
# TRUE      300     423    1103     299     830   413     723    1223   222   384

#         OPC OPC_COP Sept_In.G Sept_In.I Str_In.A Str_In.F Str_In.H Str_In.L
# FALSE 27213   27300     27109     26350    26845    26227    27248    27446
# TRUE    538     451       642      1401      906     1524      503      305

#       Thal_Ex.B TNoS_Ex.A TT.IG.SH_Ex.C TT.IG.SH_Ex.E TT.IG.SH_Ex.F Ventr_In.B
# FALSE     26120     27036         26740         27395         26282      27193
# TRUE       1631       715          1011           356          1469        558

## Print these to pngs
markerList.t.1vAll <- lapply(markers.ls.t.1vAll, function(x) {
  rownames(x[[2]])[x[[2]]$log.FDR < log(0.05) & x[[2]]$non0median == TRUE]
})

# Change to gene symbols
markerList.t.1vAll <- lapply(markerList.t.1vAll, function(x) {
  rowData(sce.ls)$gene_name[match(x, rowData(sce.ls)$gene_id)]
})

genes.top40.t <- lapply(markerList.t.1vAll, function(x) {
  head(x, n = 40)
})

for (i in names(genes.top40.t)) {
  png(here(
    "snRNAseq_mouse", "plots", "markers",
    paste0("LS_t_1vALL_top40markers-", i, "_logExprs.png")
  ), height = 1900, width = 1200)
  print(
    plotExpressionCustom(
      sce = sce.hold,
      features = genes.top40.t[[i]],
      features_name = i,
      anno_name = "cellType.final",
      ncol = 5, point_alpha = 0.4,
      scales = "free_y", swap_rownames = "gene_name"
    ) +
      scale_color_manual(values = cell_colors.ls) +
      ggtitle(label = paste0("LS ", i, " top markers: 'cluster-vs-all-others' t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}


## Write these top 40 lists to a csv ===
names(markerList.t.pw) <- paste0(names(markerList.t.pw), "_pw")
names(markerList.t.1vAll) <- paste0(names(markerList.t.1vAll), "_1vAll")

# Many of the PW results don't have 40 markers:
extend.idx <- names(which(lengths(markerList.t.pw) < 40))
for (i in extend.idx) {
  markerList.t.pw[[i]] <- c(markerList.t.pw[[i]], rep("", 40 - length(markerList.t.pw[[i]])))
}

top40genes <- cbind(
  sapply(markerList.t.pw, function(x) head(x, n = 40)),
  sapply(markerList.t.1vAll, function(y) head(y, n = 40))
)
top40genes <- top40genes[, sort(colnames(top40genes))]

write.csv(top40genes,
          file = here(
            "snRNAseq_mouse", "processed_data", "tables",
            "top40genesLists_LS-n4_33finalCellTypes.csv"
          ),
          row.names = FALSE
)






## One more iteration Jan2023 - binned by region-specific(-ish) annotation ===================
# (load SCE; this annotation is in)

# First drop any flagged clusters for dropping
sce.ls <- sce.ls[, -grep("drop.", sce.ls$cellType.broad)]
sce.ls <- sce.ls[, -grep("Neuron.mixed", sce.ls$cellType.broad)]
sce.ls$cellType.broad <- droplevels(sce.ls$cellType.broad)




### First make a list of Boolean param / cell subtype ===
# Will use this to assess more 'valid', non-noise-driving markers
cellClust.idx <- splitit(sce.ls$cellType.broad)
medianNon0.ls.broad <- lapply(cellClust.idx, function(x) {
  apply(as.matrix(assay(sce.ls, "logcounts")), 1, function(y) {
    median(y[x]) > 0
  })
})

sapply(medianNon0.ls.broad, table)
#



## Traditional t-test, pairwise ===
mod <- with(colData(sce.ls), model.matrix(~Sample))
mod <- mod[, -1, drop = F] # intercept otherwise automatically dropped by `findMarkers()`

# Run pairwise t-tests
markers.ls.t.pw.broad <- findMarkers(sce.ls,
                                     groups = sce.ls$cellType.broad,
                                     assay.type = "logcounts", design = mod, test = "t",
                                     direction = "up", pval.type = "all", full.stats = T
)

sapply(markers.ls.t.pw.broad, function(x) {
  table(x$FDR < 0.05)
})
#       Astro  Chol   ChP  Endo Ependymal   IoC    LS Micro    MS Mural Neuroblast Oligo   OPC  Sept
# FALSE 32010 31922 31390 31268     30909 32084 32229 31285 32218 31896      31914 31887 32079 32275
# TRUE    275   363   895  1017      1376   201    56  1000    67   389        371   398   206    10

#         Str  Thal  TNoS TT.IG.SH
# FALSE 32089 32061 32100    32099
# TRUE    196   224   185      186



# Add respective 'non0median' column to the stats for each set of markers
for (i in names(markers.ls.t.pw.broad)) {
  markers.ls.t.pw.broad[[i]] <- cbind(
    markers.ls.t.pw.broad[[i]],
    medianNon0.ls.broad[[i]][match(
      rownames(markers.ls.t.pw.broad[[i]]),
      names(medianNon0.ls.broad[[i]])
    )]
  )
  colnames(markers.ls.t.pw.broad[[i]])[21] <- "non0median"
}

sapply(markers.ls.t.pw.broad, function(x) {
  table(x$FDR < 0.05 & x$non0median == TRUE)["TRUE"]
})
# Astro.TRUE       Chol.TRUE        ChP.TRUE       Endo.TRUE  Ependymal.TRUE        IoC.TRUE 
#        119             126             357             360             564              98 
#    LS.TRUE      Micro.TRUE         MS.TRUE      Mural.TRUE Neuroblast.TRUE      Oligo.TRUE 
#         24             207              17              22              76             167 
#   OPC.TRUE       Sept.TRUE        Str.TRUE       Thal.TRUE       TNoS.TRUE   TT.IG.SH.TRUE 
#        104               2             116             107              86             100




markerList.t.pw.broad <- lapply(markers.ls.t.pw.broad, function(x) {
  rownames(x)[x$FDR < 0.05& x$non0median == TRUE]
})


# Change to gene symbols
markerList.t.pw.broad <- lapply(markerList.t.pw.broad, function(x) {
  rowData(sce.ls)$gene_name[match(x, rowData(sce.ls)$gene_id)]
})

genes.top40.t <- lapply(markerList.t.pw.broad, function(x) {
  head(x, n = 40)
})

smaller.set <- names(genes.top40.t)[lengths(genes.top40.t) <= 20]
left.set <- setdiff(names(genes.top40.t), smaller.set)

# Now remove those with no significant markers
smaller.set <- setdiff(
  smaller.set,
  names(genes.top40.t)[lengths(genes.top40.t) == 0]
)

# Smaller graphical window
for (i in smaller.set) {
  png(
    paste0("/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/plots/markers/LS_t_pairwise_top40markers_cellType.broad-", i, ".png"),
    height = 950, width = 1200)
  print(
    plotExpressionCustom(
      sce = sce.hold,
      features = genes.top40.t[[i]],
      features_name = i,
      anno_name = "cellType.broad",
      ncol = 5, point_alpha = 0.4,
      scales = "free_y", swap_rownames = "gene_name"
    ) +
      #scale_color_manual(values = tableau20) +
      ggtitle(label = paste0("LS ", i, " top markers: single-nucleus-level p.w. t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}

# 20-40 markers
for (i in left.set) {
  png(
    paste0("/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/plots/markers/LS_t_pairwise_top40markers_cellType.broad-", i, ".png"),
    height = 1900, width = 1200)
  print(
    plotExpressionCustom(
      sce = sce.hold,
      features = genes.top40.t[[i]],
      features_name = i,
      anno_name = "cellType.broad",
      ncol = 5, point_alpha = 0.4,
      scales = "free_y", swap_rownames = "gene_name"
    ) +
      
      #scale_color_manual(values = tableau20) +
      ggtitle(label = paste0("LS ", i, " top markers: single-nucleus-level p.w. t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}




## 1-v-rest ===
markers.ls.t.1vAll.broad <- list()
for (i in levels(sce.ls$cellType.broad)) {
  # Make temporary contrast
  sce.ls$contrast <- ifelse(sce.ls$cellType.broad == i, 1, 0)
  # Test cluster vs. all others
  markers.ls.t.1vAll.broad[[i]] <- findMarkers(sce.ls,
                                               groups = sce.ls$contrast,
                                               assay.type = "logcounts", design = mod, test = "t",
                                               std.lfc = TRUE,
                                               direction = "up", pval.type = "all", full.stats = T
  )
}

# Do some reorganizing
markers.ls.t.1vAll.broad <- lapply(markers.ls.t.1vAll.broad, function(x) {
  # Basically take the 'stats.[1 or 0]' since is redundant with the 'summary'-level stats
  lapply(x, function(y) {
    y[, 4]
  })
})

# Re-name std.lfc column and the entries; add non-0-median info
for (i in names(markers.ls.t.1vAll.broad)) {
  colnames(markers.ls.t.1vAll.broad[[i]][["0"]])[1] <- "std.logFC"
  colnames(markers.ls.t.1vAll.broad[[i]][["1"]])[1] <- "std.logFC"
  # Add non0median Boolean - might be informative for both sets of stats
  markers.ls.t.1vAll.broad[[i]][["0"]] <- cbind(
    markers.ls.t.1vAll.broad[[i]][["0"]],
    medianNon0.ls.broad[[i]][match(
      rownames(markers.ls.t.1vAll.broad[[i]][["0"]]),
      names(medianNon0.ls.broad[[i]])
    )]
  )
  colnames(markers.ls.t.1vAll.broad[[i]][["0"]])[4] <- "non0median"
  
  # "1" aka 'enriched'
  markers.ls.t.1vAll.broad[[i]][["1"]] <- cbind(
    markers.ls.t.1vAll.broad[[i]][["1"]],
    medianNon0.ls.broad[[i]][match(
      rownames(markers.ls.t.1vAll.broad[[i]][["1"]]),
      names(medianNon0.ls.broad[[i]])
    )]
  )
  colnames(markers.ls.t.1vAll.broad[[i]][["1"]])[4] <- "non0median"
  
  # Then re-name the entries to more interpretable, because we'll keeping both contrasts
  names(markers.ls.t.1vAll.broad[[i]]) <- paste0(i, c("_depleted", "_enriched"))
}


## Print these to pngs
markerList.t.1vAll.broad <- lapply(markers.ls.t.1vAll.broad, function(x) {
  rownames(x[[2]])[x[[2]]$log.FDR < log(0.05) & x[[2]]$non0median == TRUE]
})


lengths(markerList.t.1vAll.broad)
# Astro       Chol        ChP       Endo  Ependymal        IoC         LS      Micro         MS 
#   488        700       1392        890       1439        607       1925        453       1757 
# Mural Neuroblast      Oligo        OPC       Sept        Str       Thal       TNoS   TT.IG.SH 
#   253        602        388        633       1503       1592       1900        841       1790 


# Change to gene symbols
markerList.t.1vAll.broad <- lapply(markerList.t.1vAll.broad, function(x) {
  rowData(sce.ls)$gene_name[match(x, rowData(sce.ls)$gene_id)]
})

genes.top40.t <- lapply(markerList.t.1vAll.broad, function(x) {
  head(x, n = 40)
})

for (i in names(genes.top40.t)) {
  png(paste0("/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/plots/markers/LS_t_1vALL_top40markers_cellType.broad", i, ".png"),
      height = 1900, width = 1200)
  print(
    plotExpressionCustom(
      sce = sce.hold,
      features = genes.top40.t[[i]],
      features_name = i,
      anno_name = "cellType.broad",
      ncol = 5, point_alpha = 0.4,
      scales = "free_y", swap_rownames = "gene_name"
    ) +
      #scale_color_manual(values = cell_colors.ls) +
      ggtitle(label = paste0("LS ", i, "(broad) top markers: 'cluster-vs-all-others' t-tests (FDR<0.05)")) +
      theme(plot.title = element_text(size = 20))
  )
  dev.off()
}






save(markers.ls.t.pw.broad, markers.ls.t.1vAll.broad, medianNon0.ls.broad,
     file="/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/processed_data/SCE/markers-stats_LS-n4_findMarkers_33cellTypes.broad.rda")



## Write these top 40 lists to a csv ===
names(markerList.t.pw.broad) <- paste0(names(markerList.t.pw.broad), "_pw")
names(markerList.t.1vAll.broad) <- paste0(names(markerList.t.1vAll.broad), "_1vAll")

# Many of the PW results don't have 40 markers:
extend.idx <- names(which(lengths(markerList.t.pw.broad) < 40))
for (i in extend.idx) {
  markerList.t.pw.broad[[i]] <- c(markerList.t.pw.broad[[i]], rep("", 40 - length(markerList.t.pw.broad[[i]])))
}

top40genes <- cbind(
  sapply(markerList.t.pw.broad, function(x) head(x, n = 40)),
  sapply(markerList.t.1vAll.broad, function(y) head(y, n = 40))
)
top40genes <- top40genes[, sort(colnames(top40genes))]

write.csv(top40genes,
          file = "/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/processed_data/tables/top40genesLists_LS-n4_finalCellTypes_broad.csv",
          row.names = FALSE
)









rm(list = ls())

## Reproducibility information ====
print("Reproducibility information:")
Sys.time()
# [1] "2022-06-30 13:22:13 EDT"
proc.time()
#     user    system   elapsed
#   11601.10   158.64 11841.45
options(width = 120)
session_info()
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.1.2 Patched (2021-11-04 r81138)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2022-06-30
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
# rlang                  1.0.3    2022-06-27 [2] CRAN (R 4.1.2)
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

#Additions from Lionel below:

#Broad Marker Expression Heatmap for Fig 1
sce.ls <- sce.ls[, -grep("drop", sce.ls$cellType.final)]
sce.ls <- sce.ls[, -grep("mixed", sce.ls$cellType.final)]
sce.ls$cellType.final <- droplevels(sce.ls$cellType.final)

#to check if it dropped
levels(sce.ls$cellType.final)

dat <- assay(sce.ls, "logcounts")
rownames(dat) <- rowData(sce.ls)$gene_name

genes <- c("Snap25", "Slc17a7", "Slc17a6", "Gad1", "Gad2", "Mbp", "Mobp", "Pdgfra",
           "Vcan", "Bcan", "Sox4", "Slc1a2", "Aqp4", "Flt1", "Cldn5",
           "Cx3cr1", "Csf1r", "Col1a2", "Rbpms")


dat <- dat[genes, ]

cellClust.idx <- splitit(sce.ls$cellType.final)

current_dat <- do.call(cbind, lapply(cellClust.idx, function(ii) {rowMeans(dat[ , ii])}))

pdf('/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/plots/heatmap_broadMarkers_LAH.pdf', useDingbats=TRUE, height=6, width=10)

pheatmap(current_dat, cluster_rows = FALSE, cluster_cols = FALSE, breaks = seq(0.02, 3, length.out = 101),
         color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "OrRd"))(100),
         fontsize_row = 14.5, fontsize_col = 14.5)
grid::grid.text(label="log2-\nExprs", x=0.96, y=0.6, gp=grid::gpar(fontsize=10))

dev.off()

#Broad Marker Expression Heatmap for Fig 1
neuron.sce.ls <- sce.ls
neuron.sce.ls <- neuron.sce.ls[, -grep("Astro", neuron.sce.ls$cellType.final)]
neuron.sce.ls <- neuron.sce.ls[, -grep("ChP", neuron.sce.ls$cellType.final)]
neuron.sce.ls <- neuron.sce.ls[, -grep("Endo", neuron.sce.ls$cellType.final)]
neuron.sce.ls <- neuron.sce.ls[, -grep("Ependymal", neuron.sce.ls$cellType.final)]
neuron.sce.ls <- neuron.sce.ls[, -grep("Mural", neuron.sce.ls$cellType.final)]
neuron.sce.ls <- neuron.sce.ls[, -grep("Neuroblast", neuron.sce.ls$cellType.final)]
neuron.sce.ls <- neuron.sce.ls[, -grep("Oligo", neuron.sce.ls$cellType.final)]
neuron.sce.ls <- neuron.sce.ls[, -grep("OPC", neuron.sce.ls$cellType.final)]
neuron.sce.ls <- neuron.sce.ls[, -grep("Micro", neuron.sce.ls$cellType.final)]

neuron.sce.ls$cellType.final <- droplevels(neuron.sce.ls$cellType.final)

levels(neuron.sce.ls$cellType.final)

#Generate heatmap for only neurons
dat <- assay(neuron.sce.ls, "logcounts")
rownames(dat) <- rowData(neuron.sce.ls)$gene_name

genes <- c("Sv2b", "Pcsk5", "Satb2", "Samd3", "Sema3a", "Synpo2", "Bcl11b", "Rarb",
           "Tac1", "Penk", "Trpc5", "Elavl2", "Trpc4","Homer2", "Ptpn3", "Trhde", "Cpne7",
           "Nrp1", "Pkib", "Drd3", "Chat")

dat <- dat[genes, ]

cellClust.idx <- splitit(neuron.sce.ls$cellType.final)

current_dat <- do.call(cbind, lapply(cellClust.idx, function(ii) {rowMeans(dat[ , ii])}))

pdf('/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/plots/heatmap_NeuronalRegionMarkers 2.pdf', useDingbats=TRUE, height=6, width=16)

pheatmap(current_dat, cluster_rows = FALSE, cluster_cols = FALSE, breaks = seq(0.02, 3, length.out = 101),
         color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "OrRd"))(100),
         fontsize_row = 14.5, fontsize_col = 14.5)
grid::grid.text(label="", x=0.96, y=0.6, gp=grid::gpar(fontsize=10))

dev.off()

#Generate heatmap for only LS clusters
ls.sce.ls <- neuron.sce.ls
ls.sce.ls <- ls.sce.ls[, -grep("MS", ls.sce.ls$cellType.final)]
ls.sce.ls <- ls.sce.ls[, -grep("Chol", ls.sce.ls$cellType.final)]
ls.sce.ls <- ls.sce.ls[, -grep("IoC", ls.sce.ls$cellType.final)]
ls.sce.ls <- ls.sce.ls[, -grep("Str", ls.sce.ls$cellType.final)]
ls.sce.ls <- ls.sce.ls[, -grep("Thal", ls.sce.ls$cellType.final)]
ls.sce.ls <- ls.sce.ls[, -grep("TNoS", ls.sce.ls$cellType.final)]
ls.sce.ls <- ls.sce.ls[, -grep("TT.IG.SH", ls.sce.ls$cellType.final)]

ls.sce.ls$cellType.final <- droplevels(ls.sce.ls$cellType.final)

levels(ls.sce.ls$cellType.final)


#Generate heatmap for only LS clusters
dat <- assay(ls.sce.ls, "logcounts")
rownames(dat) <- rowData(ls.sce.ls)$gene_name

genes <- c("Ar", "Pgr", "Esr1", "Esr2", "Crhr1", "Crhr2", "Ghr", "Nr3c1", "Nr3c2")

dat <- dat[genes, ]

cellClust.idx <- splitit(ls.sce.ls$cellType.final)

current_dat <- do.call(cbind, lapply(cellClust.idx, function(ii) {rowMeans(dat[ , ii])}))

pdf('/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/plots/heatmap_LSClusters_Hormone Receptors.pdf', useDingbats=TRUE, height=6, width=7)

pheatmap(current_dat, cluster_rows = FALSE, cluster_cols = FALSE, breaks = seq(0.02, 2.5, length.out = 101),
         color = colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "OrRd"))(100),
         fontsize_row = 14.5, fontsize_col = 14.5)
grid::grid.text(label="log2-\nExprs", x=0.96, y=0.6, gp=grid::gpar(fontsize=10))

dev.off()



#Monoamine Receptors: c("Htr1a", "Htr1b", "Htr1d", "Htr1f", "Htr2a", "Htr2b", "Htr2c", "Htr4", "Htr5a", "Htr5b", "Htr6", "Htr7",
#"Adra1a", "Adra1b", "Adra1d","Adra2a", "Adra2b", "Adra2c", "Adrb1", "Adrb2", "Adrb3", "Drd1", "Drd2", "Drd3")

#Cluster Markers: c("Kcnmb2", "Gpr176", "Col19a1", "Tac1", "Nos1", "Npy", "Sst", "Stac2", "Grid2ip", "Trpc6", "Trpc3", "Npas1", "Cnr1", "Kcnq4", "Reln", 
#"Pax6", "Gpc3", "Gabrg1", "Rnf207", "St8sia6", "Crhr2", "Tafa1", "Htr1b", "Slc18a2", "Col15a1", "Ano2", "Ntf3", "Vipr2", "Cpa6", "Lgr5")

#GABAergic Markers: c("Slc17a6", "Slc17a7", "Gad1", "Gad2", "Sst", "Nts", "Pvalb", "Calb1", "Calb2", "Npy", "Nos1", "Penk", "Tac1", "Pdyn")

