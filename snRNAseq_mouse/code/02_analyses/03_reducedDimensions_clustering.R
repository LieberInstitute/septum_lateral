### mouse LS snRNA-seq analysis
### Feature selection with deviance & clustering
###     qrsh -l bluejay,mf=92G,h_vmem=96G,h_fsize=40G
### Initiated LAR,ST,MNT 22Mar2022


library(SingleCellExperiment)
library(DropletUtils)
library(scater)
library(scran)
library(scry)
library(batchelor)
library(bluster)
library(BiocParallel)
library(jaffelab)
library(gridExtra)
library(here)
library(sessioninfo)


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
here()
# [1] "/dcs04/lieber/marmaypag/pilotLS_LIBD1070"

source("/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/analyses_sn/plotExpressionCustom.R")


### Feature selection with deviance residuals =====

# Load working SCE
load(here("snRNAseq_mouse", "processed_data", "SCE", "sce_working_LS.rda"), verbose = T)
# sce.ls, sce.ls.unfiltered

sce.ls
# class: SingleCellExperiment
# dim: 32285 22860
# metadata(1): Samples
# assays(1): counts
# rownames(32285): ENSMUSG00000051951 ENSMUSG00000089699 ...
#   ENSMUSG00000095019 ENSMUSG00000095041
# rowData names(6): source type ... gene_name gene_type
# colnames(22860): 1_AAACCCAAGGTACATA-1 1_AAACCCACATCCGAGC-1 ...
#   4_TTTGTTGCATACAGCT-1 4_TTTGTTGGTCAAACGG-1
# colData names(11): Sample Barcode ... high.mito doubletScore
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):

sce.ls <- devianceFeatureSelection(sce.ls,
                                   assay = "counts", fam = "binomial", sorted = F,
                                   # these are default params btw
                                   batch = as.factor(sce.ls$Sample)
)


# Btw:
table(is.na(rowData(sce.ls)$binomial_deviance))
# FALSE  TRUE
# 29556  7045

# Observe:
pdf(here("snRNAseq_mouse", "plots", "featureSelxn_binomialDeviance-byGene.pdf"), height = 5)
plot(sort(rowData(sce.ls)$binomial_deviance, decreasing = T),
     type = "l", xlab = "ranked genes",
     ylab = "binomial deviance", main = "Feature Selection with Deviance: LS (n=4)"
)
abline(v = 2000, lty = 2, col = "red")
dev.off()


## Approximate GLM-PCA with PCA on null deviance residuals
Sys.time()
# [1] "2022-03-23 10:37:02 EDT"
sce.ls <- nullResiduals(sce.ls,
                        assay = "counts", fam = "binomial", # default params
                        # type="deviance")#, batch=as.factor(sce.ls$Sample))
                        type = "pearson"
)
# MNT comment: previously in other projects, 'batch=' threw an error
#             (and we perform MNN-batch correction anyway, if relevant)

#           - Using Pearson residuals bc deviance took nearly an hour
#             with up to 96G RAM...
Sys.time()
# [1] "2022-03-23 11:52:14 EDT"

# *The above adds a new computed assay called 'binomial_deviance_residuals':
sce.ls
# class: SingleCellExperiment
# dim: 32285 22860
# metadata(1): Samples
# assays(2): counts binomial_pearson_residuals
# rownames(32285): ENSMUSG00000051951 ENSMUSG00000089699 ...
#   ENSMUSG00000095019 ENSMUSG00000095041
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(22860): 1_AAACCCAAGGTACATA-1 1_AAACCCACATCCGAGC-1 ...
#   4_TTTGTTGCATACAGCT-1 4_TTTGTTGGTCAAACGG-1
# colData names(11): Sample Barcode ... high.mito doubletScore
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):


## Save into a new file, call 'updated' as opposed to 'working':
save(sce.ls,
     file = here("snRNAseq_mouse", "processed_data", "SCE", "sce_updated_LS.rda")
)



## Feature selection: Take top 2000 highly-deviant genes ('HDG's) for PCA
#      (based on total binomial deviance, computed above)
hdgs.ls <- rownames(sce.ls)[order(rowData(sce.ls)$binomial_deviance, decreasing = T)][1:2000]
hdgs.symbols <- rowData(sce.ls)$gene_name[match(hdgs.ls, rowData(sce.ls)$gene_id)]
# Out of curiosity
c("Snap25", "Mbp", "Gad1") %in% hdgs.symbols # all TRUE

# Run PCA:
Sys.time()
# [1] "2022-03-23 12:44:41 EDT"
set.seed(109)
sce.ls <- runPCA(sce.ls,
                 # exprs_values="binomial_deviance_residuals",
                 exprs_values = "binomial_pearson_residuals",
                 subset_row = hdgs.ls, ncomponents = 100,
                 name = "GLMPCA_approx",
                 BSPARAM = BiocSingular::IrlbaParam()
)
Sys.time()
# [1] "2022-03-23 14:22:15 EDT"      - this took ~1hr40min to compute...



## **Optional: batch correction with batchelor::reducedMNN ===
# Sys.time()
#         #
# glmpca.mnn <- reducedMNN(reducedDim(sce.ls, "GLMPCA_approx"),
#                          batch=as.factor(sce.ls$Sample),
#                          merge.order=c("Br8079_LC", "Br2701_LC", "Br6522_LC")
# )
# Sys.time()
#         #
#
# # Store this
# reducedDim(sce.ls, "GLMPCA_MNN") <- glmpca.mnn$corrected


## Make some other reduced dims for visualization ===

# # Take top 50 PCs to quicken computation
# reducedDim(sce.ls, "glmpca_mnn_50") <- reducedDim(sce.ls, "GLMPCA_MNN")[ ,1:50]
# # oh didn't have to use this - can just use 'n_dimred'
#
# UMAP
set.seed(109)
sce.ls <- runUMAP(sce.ls,
                  dimred = "GLMPCA_approx",
                  n_dimred = 50, name = "UMAP"
)
# t-SNE
set.seed(109)
sce.ls <- runTSNE(sce.ls,
                  dimred = "GLMPCA_approx",
                  n_dimred = 50, name = "TSNE"
)

# Visualize top PCs and these 2D embeddings:
pdf(here("snRNAseq_mouse", "plots", "reducedDims_mouseLS-n4_noBatchCorrxn.pdf"))
plotReducedDim(sce.ls,
               dimred = "GLMPCA_approx", colour_by = "Sample",
               ncomponents = 4, point_alpha = 0.3, point_size = 1.5
)
# UMAPs, along with some other metrics
plotReducedDim(sce.ls,
               dimred = "UMAP", colour_by = "Sample",
               point_alpha = 0.3, point_size = 1.5
)
plotReducedDim(sce.ls,
               dimred = "UMAP", colour_by = "sum",
               point_alpha = 0.3, point_size = 1.5
)
plotReducedDim(sce.ls,
               dimred = "UMAP", colour_by = "doubletScore",
               point_alpha = 0.3, point_size = 1.5
)
# TSNE
plotReducedDim(sce.ls,
               dimred = "TSNE", colour_by = "Sample",
               point_alpha = 0.3, point_size = 1.5
)
dev.off()

## STOP: Are there batch effects that we can diagnose by eye?
##     - If so, proceed in the batch corrected space:
# # UMAP
# set.seed(109)
# sce.ls <- runUMAP(sce.ls, dimred="glmpca_mnn_50",
#                   n_dimred=50, name="UMAP_mnn")
# # t-SNE
# set.seed(109)
# sce.ls <- runTSNE(sce.ls, dimred="glmpca_mnn_50",
#                   n_dimred=50, name="TSNE_mnn")



## Save for now
save(sce.ls,
     file = here("snRNAseq_mouse", "processed_data", "SCE", "sce_updated_LS.rda")
)



### Clustering ====================
# Perform graph-based clustering, as in Tran-Maynard, et al. Neuron 2021

reducedDim(sce.ls, "GLMPCA_50") <- reducedDim(sce.ls, "GLMPCA_approx")[, 1:50]

snn.gr.glmpca <- buildSNNGraph(sce.ls, k = 20, use.dimred = "GLMPCA_50")
clusters.glmpca <- igraph::cluster_walktrap(snn.gr.glmpca)
table(clusters.glmpca$membership)
#   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16
# 797  760  191  451   63 4153 1478 4327  166  464  670  306 1529  112 4364  118
#  17   18   19   20   21   22   23   24   25   26   27   28   29   30   31   32
# 106  105   65  232   67  171  363  146  202  728   40  180   34   86   73  150
#  33   34   35
#  53   34   76

# Store this info
sce.ls$clusters.glmpca <- factor(clusters.glmpca$membership)

# Also save the graph & community info
save(snn.gr.glmpca, clusters.glmpca,
     file = here("snRNAseq_mouse", "processed_data", "SCE", "graph_clusters_glmpca_LS-n3.rda")
)

table(sce.ls$clusters.glmpca, sce.ls$Sample)
# pretty even distribution across samples

## doubletScore distributions / cluster?
cellClust.idx <- splitit(sce.ls$clusters.glmpca)
sapply(cellClust.idx, function(x) {
  round(quantile(sce.ls$doubletScore[x]), 2)
})
#          1     2    3     4    5     6    7     8    9   10    11   12    13
# 0%    0.01  0.00 0.01  0.00 0.01  0.00 0.00  0.00 0.00 0.00  0.02 0.00  0.00
# 25%   0.13  0.03 0.23  0.04 0.08  0.07 0.02  0.15 0.06 0.01  0.15 0.02  0.05
# 50%   0.25  0.07 0.43  0.07 0.12  0.20 0.07  0.34 0.09 0.04  0.39 0.06  0.15
# 75%   0.62  0.78 1.08  0.26 0.41  0.52 0.18  1.16 0.17 0.91  1.31 0.52  0.41
# 100% 10.44 11.46 5.98 12.53 8.94 14.84 7.20 11.49 3.62 9.30 10.41 6.60 14.29

#        14    15   16   17   18   19   20   21    22   23   24   25   26   27
# 0%   0.00  0.00 0.00 0.01 0.04 0.02 0.03 0.04  0.01 0.00 0.00 0.02 0.04 0.03
# 25%  0.04  0.16 0.04 0.04 0.17 0.36 0.14 0.38  0.06 0.04 0.04 0.13 0.17 0.07
# 50%  0.05  0.29 0.05 0.08 0.24 0.49 0.46 1.28  0.22 0.11 0.07 0.23 0.30 0.10
# 75%  0.11  0.90 0.16 0.19 0.43 0.53 1.30 1.39  1.99 0.84 0.22 0.52 1.03 0.19
# 100% 8.53 10.09 5.98 8.16 6.46 7.10 6.82 9.93 10.39 9.50 5.37 7.62 9.25 3.39

#        28   29   30   31  *32   33    34   35
# 0%   0.00 0.02 0.09 0.04 0.06 0.05  0.00 0.00
# 25%  0.01 0.05 0.30 0.17 3.78 0.24  0.06 0.01
# 50%  0.03 0.06 0.52 0.33 4.83 0.66  0.17 0.04   Cluster 32 (150 nuclei) to keep
# 75%  0.06 0.10 0.91 0.69 6.72 0.93  0.30 0.05   an eye on...
# 100% 7.93 9.43 8.40 5.91 9.39 8.48 10.09 6.53

# Compute logcounts to visualize expression in the traditional way
sce.ls <- multiBatchNorm(sce.ls, batch = sce.ls$Sample)

# # For sake of filtering, compute some median expression info across graph-based clusters
# medianNon0.ls <- lapply(cellClust.idx, function(x){
#     apply(as.matrix(assay(sce.ls, "logcounts")), 1, function(y){
#         median(y[x]) > 0
#     })
# })
#
# head(medianNon0.ls[[1]])

# Save
save(sce.ls,
     file = here("snRNAseq_mouse", "processed_data", "SCE", "sce_updated_LS.rda")
)


## Broad markers of interest:
markers.mathys.tran <- list(
  "neuron" = c("SYT1", "SNAP25", "GRIN1"),
  "excit_neuron" = c("SLC17A7", "SLC17A6", "SLC17A8"),
  "inhib_neuron" = c("GAD1", "GAD2", "SLC32A1"),
  # Norepinephrine & serotonergic markers
  "neuron.NE" = c("TH", "DBH", "SLC6A2", "SLC18A2", "GCH1", "DDC"), # SLC6A3 - saw no DAT
  "neuron.5HT" = c("SLC6A4", "TPH1", "TPH2", "DDC"),
  # SERT, serotonin T (aka 5-HTT);
  "monoamine.metab" = c("COMT", "MAOA", "MAOB"),
  # MSN markers
  "MSNs.pan" = c("PPP1R1B", "BCL11B"), # "CTIP2")
  "MSNs.D1" = c("DRD1", "PDYN", "TAC1"),
  "MSNs.D2" = c("DRD2", "PENK"),
  ## Non-neuronal:
  "oligodendrocyte" = c("MBP", "MOBP", "PLP1"),
  "oligo_precursor" = c("PDGFRA", "VCAN", "CSPG4"),
  "microglia" = c("CD74", "CSF1R", "C3"),
  "astrocyte" = c("GFAP", "TNC", "AQP4", "SLC1A2"),
  "endothelial" = c("CLDN5", "FLT1", "VTN"),
  # Post-hoc from Tran-Maynard, et al. Neuron 2021
  "differn_committed_OPC" = c("SOX4", "BCAN", "GPR17", "TNS3"),
  "Tcell" = c("SKAP1", "ITK", "CD247"),
  "Mural" = c("COL1A2", "TBX18", "RBPMS"),
  "Macro" = c("CD163", "SIGLEC1", "F13A1")
)


# Will have to 'make these mouse'
broadMarkers <- markers.mathys.tran
for (i in 1:length(broadMarkers)) {
  broadMarkers[[i]] <- paste0(
    substr(broadMarkers[[i]], 1, 1),
    tolower(substr(broadMarkers[[i]], 2, nchar(broadMarkers[[i]])))
  )
}

table(unname(unlist(broadMarkers)) %in% rowData(sce.ls)$gene_name) # all good


## Expected region/sub-region markers / other curated sets ===
markers.curated <- list(
  "Accumbens" = c(
    "Ppp1r1b", "Adora2a", "Rgs9", "Syndig1l", "Ric8b",
    "Rgs4", "Rgs7bp", "Dgki"
  ),
  "Pan.LS" = c(
    "Dgkg", "Dgkh", "Prkcd", "Glp1r", "Trpc4", "Zic1",
    "Homer2", "Ptpn3"
  ),
  "MS.specific" = c(
    "Nacc2", "Sgpp2", "Kcnab3", "Ngfr", "Lgi2", "Nrip3",
    "Lrrc55", "Trpc5", "Tshz3", "Scn1a"
  ),
  "Septal.Hip" = c("Grid2ip", "Slc17a7", "Matn2", "Bok", "Egr4"),
  "IG" = c("Cabp7", "Them6", "Kcnk2"),
  "TT.IG" = c("Sv2b", "Mtmr12", "Zbtb16", "Cck"),
  "Ventricle.ependymal" = c("Lrrc74b"),
  "TNoS" = c("Eomes", "Col9a1", "Eps8", "Sln", "Adgrd1"),
  "BNST" = c("Tac2", "Glra3", "Tmem145", "Fxyd7", "Rcn1"),
  "Diag.band" = c(
    "Ntrk1", "Coro6", "Itih3", "Arhgap12", "Chat", "Kcnc2",
    "Elavl2", "Lgi2", "Nacc2", "Syt17", "Rnf227", "Ngfr", "Slc18a3"
  )
)

markers.curated.2 <- list(
  "GABA_sub" = c("Sst", "Pvalb", "Cck", "Npy", "Calb1", "Calb2", "Vip", "Nts"),
  "Peptides" = c("Penk", "Pdyn", "Tac1", "Ghrh"),
  "PeptideR" = c(
    "Oxtr", "Crhr1", "Crhr2", "Avpr1a", "Avpr1b", "Avpr2",
    "Galr1", "Galr2", "Galr3", "Ghrhr", "Nr3c2"
  ),
  "SexHormoneR" = c("Ar", "Esr1", "Esr2", "Pgr"),
  "MetaboSerotoninR" = c(
    "Htr1a", "Htr1b", "Htr1d", "Htr1f",
    "Htr2a", "Htr2b", "Htr2c", "Htr4", "Htr5a",
    "Htr5b", "Htr6", "Htr7"
  ),
  "IonoSerotoninR" = c("Htr3a", "Htr3b"),
  "AdrenergicR" = c(
    "Adra1a", "Adra1b", "Adra1d", "Adra2a",
    "Adra2b", "Adra2c", "Adrb1", "Adrb2", "Adrb3"
  ),
  "DopamineR" = c("Drd1", "Drd2", "Drd3", "Drd4", "Drd5"),
  "Other" = c("Ntrk2", "Bdnf")
)


markers.besnard <- list(
  "Rostrocaudal.Deep" = c("Cbln2", "Cbln4", "Pax6", "Sema3a"),
  "Rostrocaudal.Superficial" = c("Col15a1", "Matn2", "Wfs1"),
  "Rostral.Deep" = c("Drd3", "Galr1", "Igfbp5", "Pou6f2"),
  "Rostral.Superficial" = c(
    "Foxp2", "Ndst4", "Gdpd2", "Cdh7",
    "Crhr2", "Dach2", "Fst", "Lhx2"
  ),
  "Caudal.Deep" = c("Asb4", "Otx2", "Pde1c", "Tacr1", "Trpc6"),
  "Caudal.Superficial" = c("Cd24a", "Igfbp4")
)



# pdf(here("snRNAseq_mouse","plots",paste0("LS-n4_expression_broadMarkers_GLMPCA-graphClusters.pdf")),
# pdf(here("snRNAseq_mouse","plots",paste0("LS-n4_expression_broadMarkers_GLMPCA-graphClusters_annotated.pdf")),
pdf(here("snRNAseq_mouse", "plots", paste0("LS-n4_expression_broadMarkers_GLMPCA-graphClusters_finalAnnotations.pdf")), height = 6, width = 14)

# pdf(here("snRNAseq_mouse","plots",paste0("LS-n4_expression_curatedMarkers_GLMPCA-graphClusters_finalAnnotations.pdf")),
for (i in 1:length(broadMarkers)) {
  
  # pdf(here("snRNAseq_mouse","plots","supplemental", paste0("LS-n4_expression_", names(markers.curated)[i],"_Markers_finalAnnotations.pdf")),
  #       height=6, width=14)
  print(
    plotExpressionCustom(
      sce = sce.ls,
      exprs_values = "logcounts",
      features = broadMarkers[[i]],
      features_name = names(broadMarkers)[i],
      # anno_name = "clusters.glmpca",
      anno_name = "cellType.final",
      ncol = 4,
      point_alpha = 0.4, point_size = 0.9,
      scales = "free_y", swap_rownames = "gene_name"
    ) +
      ggtitle(label = paste0(
        "mouse LS (n4) clusters: ",
        names(broadMarkers)[i], " markers from Mathys, et al. (2019) & Tran, Maynard, et al. (2021)"
      )) +
      theme(
        plot.title = element_text(size = 12),
        axis.text.x = element_text(size = 7)
      ) +
      scale_color_manual(values = cell_colors.ls)
  )
}
dev.off()


# Add'ly temp/for Lionel:
pdf(here("snRNAseq_mouse", "plots", paste0("temp_interactiveMarkerExplore.pdf")), height = 6, width = 14)
plotExpressionCustom(
  sce = sce.ls,
  exprs_values = "logcounts",
  # Ependymal cell markers
  features = c("Rarres2", "Ccdc153", "Tmem212", "S100b", "Acta2", "Foxj1"),
  features_name = "custom-selected",
  anno_name = "cellType",
  ncol = 3, point_alpha = 0.4, point_size = 0.9,
  scales = "free_y", swap_rownames = "gene_name"
) +
  scale_color_manual(values = c(cbPalette, tableau20, tableau10medium)) +
  ggtitle(label = paste0("Literature-based selection of ependymal markers")) +
  theme(
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7)
  )
#
plotExpressionCustom(
  sce = sce.ls,
  exprs_values = "logcounts",
  # Ependymal cell markers
  features = c("Drd3", "Htr4"),
  features_name = "custom-selected",
  anno_name = "cellType",
  ncol = 3, point_alpha = 0.4, point_size = 0.9,
  scales = "free_y", swap_rownames = "gene_name"
) +
  scale_color_manual(values = c(cbPalette, tableau20, tableau10medium)) +
  ggtitle(label = paste0("Lionel's custom markers of interest")) +
  theme(
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7)
  )
dev.off()



## Annotate clusters ===
annotationTab.ls <- data.frame(cluster = c(1:35))
annotationTab.ls$cellType <- NA
annotationTab.ls$cellType[c(4, 11, 20, 21, 34, 25, 26)] <- paste0("Excit_", c("A", "B", "C", "D", "E", "F", "G"))
annotationTab.ls$cellType[c(1, 2, 3, 8, 12, 15, 23, 30)] <-
  paste0("Inhib_", c("A", "B", "C", "D", "E", "F", "G", "H"))
annotationTab.ls$cellType[c(31)] <- "Neuron.mixed_A"

annotationTab.ls$cellType[c(6, 13)] <- paste0("Astro_", c("A", "B"))
annotationTab.ls$cellType[c(17, 24, 5, 22)] <- paste0("Aqp4.Rbpms_", c("A", "B", "C", "D"))
annotationTab.ls$cellType[c(16, 27)] <- paste0("Endo_", c("A", "B"))
annotationTab.ls$cellType[c(14, 29, 35)] <- paste0("Mural_", c("A", "B", "C"))
annotationTab.ls$cellType[c(9)] <- "Micro"
annotationTab.ls$cellType[c(7, 18)] <- paste0("Oligo_", c("A", "B"))
annotationTab.ls$cellType[10] <- c("OPC")
annotationTab.ls$cellType[33] <- c("OPC_COP")
annotationTab.ls$cellType[c(19)] <- "Ependymal"

annotationTab.ls$cellType[28] <- c("drop.lowNTx")
annotationTab.ls$cellType[32] <- c("drop.doublet")


sce.ls$cellType <- annotationTab.ls$cellType[match(
  sce.ls$clusters.glmpca,
  annotationTab.ls$cluster
)]
sce.ls$cellType <- factor(sce.ls$cellType)

table(sce.ls$cellType)
#  Aqp4.Rbpms_A   Aqp4.Rbpms_B   Aqp4.Rbpms_C   Aqp4.Rbpms_D        Astro_A        Astro_B
#           106            146             63            171           4153           1529
#  drop.doublet    drop.lowNTx         Endo_A         Endo_B      Ependymal        Excit_A
#           150            180            118             40             65            451
#       Excit_B        Excit_C        Excit_D        Excit_E        Excit_F        Excit_G
#           670            232             67             34            202            728
#       Inhib_A        Inhib_B        Inhib_C        Inhib_D        Inhib_E        Inhib_F
#           797            760            191           4327            306           4364
#       Inhib_G        Inhib_H          Micro        Mural_A        Mural_B        Mural_C
#           363             86            166            112             34             76
# Neuron.mixed_A        Oligo_A        Oligo_B            OPC        OPC_COP
#            73           1478            105            464             53


## Save
save(sce.ls, annotationTab.ls,
     file = here("snRNAseq_mouse", "processed_data", "SCE", "sce_updated_LS.rda")
)


# 21-22Apr (now 24Jun): And now since we're happy with the decided drop. clusters, drop those
#                       and re-print the reducedDims & marker expression plots
sce.ls <- sce.ls[, -grep("drop.", sce.ls$cellType.final)]
sce.ls <- sce.ls[, -grep("Neuron.mixed", sce.ls$cellType.final)]
sce.ls$cellType.final <- droplevels(sce.ls$cellType.final)

# 24Jun: Drop those from the cell_colors.ls too
#     (based on observations at https://github.com/lmweber/locus-c/blob/f2509bb1fc74c01155c9f3a74faf222bf2f1f1a1/code/analyses_sn/04_clusterAgglomeration.R#L348)
cell_colors.ls <- cell_colors.ls[-grep("drop", names(cell_colors.ls))]
cell_colors.ls <- cell_colors.ls[-grep("Neuron.mixed", names(cell_colors.ls))]



## re-print reducedDims with these new annotations ===
##As needed, remove the dropped clusters

#pdf(here("snRNAseq_mouse", "plots", "reducedDims_mouseLS-n4_graph-basedClusters_annotated.pdf"))
pdf("/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/plots/reducedDims_mouseLS-n4_graph-basedClusters_FINALannot_Nov2022.pdf")

# UMAPs, along with some other metrics
plotReducedDim(sce.ls,
               dimred = "UMAP", colour_by = "Sample",
               point_alpha = 0.3, point_size = 1.5
) +
  ggtitle("UMAP of LS (n=4)")
plotReducedDim(sce.ls,
               dimred = "UMAP", colour_by = "Sex",
               point_alpha = 0.3, point_size = 1.5
) +
  ggtitle("UMAP of LS (n=4)")
plotReducedDim(sce.ls,
               dimred = "UMAP", colour_by = "sum",
               point_alpha = 0.3, point_size = 1.5
) +
  ggtitle("UMAP of LS (n=4)")
plotReducedDim(sce.ls,
               dimred = "UMAP", colour_by = "doubletScore",
               point_alpha = 0.3, point_size = 1.5
) +
  ggtitle("UMAP of LS (n=4)")
# With text
plotReducedDim(sce.ls,
               dimred = "UMAP", colour_by = "cellType.final",
               text_by = "cellType.final",
               text_size = 3,
               point_alpha = 0.3, point_size = 1.5
) +
  scale_color_manual(
    values = cell_colors.ls,
    labels = paste0(
      levels(sce.ls$cellType.final), " (",
      table(sce.ls$cellType.final), ")"
    )
  ) +
  guides(color = guide_legend(ncol = 1)) +
  ggtitle("UMAP of LS (n=4), colored by annotated graph-based clusters") +
  labs(colour = "Cell type")
# Without text
plotReducedDim(sce.ls,
               dimred = "UMAP", colour_by = "cellType.final",
               text_size = 3,
               point_alpha = 0.3, point_size = 1.5
) +
  scale_color_manual(
    values = cell_colors.ls,
    labels = paste0(
      levels(sce.ls$cellType.final), " (",
      table(sce.ls$cellType.final), ")"
    )
  ) +
  guides(color = guide_legend(ncol = 1)) +
  ggtitle("UMAP of LS (n=4), colored by annotated graph-based clusters") +
  labs(colour = "Cell type")

# TSNE
plotReducedDim(sce.ls,
               dimred = "TSNE", colour_by = "cellType.final",
               text_by = "cellType.final", text_size = 3,
               point_alpha = 0.3, point_size = 1.5
) +
  scale_color_manual(
    values = cell_colors.ls,
    labels = paste0(
      levels(sce.ls$cellType.final), " (",
      table(sce.ls$cellType.final), ")"
    )
  ) +
  guides(color = guide_legend(ncol = 1)) +
  ggtitle("TSNE of LS (n=4), colored by annotated graph-based clusters") +
  labs(colour = "Cell type")

plotReducedDim(sce.ls,
               dimred = "TSNE", colour_by = "cellType.final",
               text_size = 3,
               point_alpha = 0.3, point_size = 1.5
) +
  scale_color_manual(
    values = cell_colors.ls,
    labels = paste0(
      levels(sce.ls$cellType.final), " (",
      table(sce.ls$cellType.final), ")"
    )
  ) +
  guides(color = guide_legend(ncol = 1)) +
  ggtitle("TSNE of LS (n=4), colored by annotated graph-based clusters") +
  labs(colour = "Cell type")
dev.off()


## Look at some distributions with these annotations
cellClust.idx <- splitit(sce.ls$cellType)
sapply(cellClust.idx, function(x) {
  round(quantile(sce.ls$sum[x]), 2)
})

sapply(cellClust.idx, function(x) {
  round(quantile(sce.ls$doubletScore[x]), 2)
})



### Inhib_D: huge (>4300 nuclei) cluster with little PW markers =======
# Sub-clustering likely necessary
library(pheatmap)

# Compute cluster modularity ratio (a measure of cluster separation)
mod.ratio.merged.HC <- pairwiseModularity(
  graph = snn.gr.glmpca,
  clusters = sce.ls$cellType.final,
  as.ratio = TRUE
)

# Heatmap
# pdf(here("snRNAseq_mouse","plots","clusterModRatio_35clusters_LS-n4.pdf"))
pdf(here("snRNAseq_mouse", "plots", "clusterModRatio_37clusters_LS-n4.pdf"))
pheatmap(log2(mod.ratio.merged.HC + 1),
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("white", "blue"))(100),
         main = "Modularity ratio for 37 graph-based clusters in LS (n=4)",
         fontsize_row = 7.5, fontsize_col = 7.5, angle_col = 90,
         display_numbers = T, number_format = "%.1f", fontsize_number = 5.5,
         na_col = "darkgrey"
)
grid::grid.text(label = "log2(ratio)", x = 0.97, y = 0.64, gp = grid::gpar(fontsize = 7))
dev.off()
# Inhib_D for sure should be re-clustered; maybe Inhib_F for the low
# within-cluster modularity ratio, and that it is an equally large cluster



## Subclustering Inhib_D ===

# Try scran's `clusterCells` with defaults (NNGraphParam input)
subclust.Inhib_D <- clusterCells(sce.ls[, sce.ls$cellType == "Inhib_D"], use.dimred = "GLMPCA_50")
# 23 populations with default k=10 neighbors & walktrap method

subclust.Inhib_D <- clusterCells(sce.ls[, sce.ls$cellType == "Inhib_D"],
                                 use.dimred = "GLMPCA_50",
                                 BLUSPARAM = NNGraphParam(k = 20)
)
# 13 populations

sce.inhibD <- sce.ls[, sce.ls$cellType == "Inhib_D"]
sce.inhibD$cellType <- droplevels(sce.inhibD$cellType)
sce.inhibD$cellType.sub <- clusterCells(sce.inhibD,
                                        use.dimred = "GLMPCA_50",
                                        BLUSPARAM = NNGraphParam(k = 20)
)

pdf(here("snRNAseq_mouse", "plots", paste0("LS-n4_expression_broadMarkers_Inhib_D-subClusters.pdf")),
    height = 6, width = 14
)
for (i in 1:length(broadMarkers)) {
  print(
    plotExpressionCustom(
      sce = sce.inhibD,
      exprs_values = "logcounts",
      features = broadMarkers[[i]],
      features_name = names(broadMarkers)[[i]],
      # anno_name = "clusters.glmpca",
      anno_name = "cellType.sub",
      ncol = 4,
      point_alpha = 0.4, point_size = 0.9,
      scales = "free_y", swap_rownames = "gene_name"
    ) +
      ggtitle(label = paste0(
        "mouse LS (n4) subclusters of 'Inhib_D': ",
        names(broadMarkers)[[i]], " markers"
      )) +
      theme(
        plot.title = element_text(size = 12),
        axis.text.x = element_text(size = 7)
      )
  )
}
dev.off()

## These look pretty good - will store these and merge non-neuronal populations

sce.ls$cellType.final <- as.character(sce.ls$cellType)
sce.ls$cellType.final[grep("Aqp4.Rbpms", sce.ls$cellType.final)] <- "Aqp4.Rbpms"
sce.ls$cellType.final[grep("Mural", sce.ls$cellType.final)] <- "Mural"
sce.ls$cellType.final[grep("Oligo", sce.ls$cellType.final)] <- "Oligo"
sce.ls$cellType.final[grep("Astro", sce.ls$cellType.final)] <- "Astro"
sce.ls$cellType.final[grep("Endo", sce.ls$cellType.final)] <- "Endo"


## What if just didn't subset the SCE?
sce.ls$cellType.final[sce.ls$cellType.final == "Inhib_D"] <-
  as.character(clusterCells(sce.ls[, sce.ls$cellType.final == "Inhib_D"],
                            use.dimred = "GLMPCA_50",
                            BLUSPARAM = NNGraphParam(k = 20)
  ))

# Then check this
sce.inhibD <- sce.ls[, sce.ls$cellType == "Inhib_D"]

plotExpressionCustom(
  sce = sce.inhibD,
  exprs_values = "logcounts",
  features = broadMarkers[["excit_neuron"]],
  features_name = "",
  anno_name = "cellType.final",
  ncol = 4,
  point_alpha = 0.4, point_size = 0.9,
  scales = "free_y", swap_rownames = "gene_name"
)
# this looks right
#   -> go ahead and re-assign these to new annotations (mostly inhib)

# First - Inhib_D.1 looks very much like oligos....
plotReducedDim(sce.inhibD, "UMAP", colour_by = "cellType.final", point_size = 2.5)
# it's quite dispersed
cell.sub.idx <- splitit(sce.inhibD$cellType.final)
sapply(cell.sub.idx, function(x) {
  quantile(sce.inhibD$sum[x])
})
# total UMI distribution is normal for oligos
sapply(cell.sub.idx, function(x) {
  round(quantile(sce.inhibD$doubletScore[x]), 3)
})
#           1    10    11    12    13      2      3     4      5     6     7     8     9
# 0%    0.000 0.134 0.022 0.036 0.073  0.007  0.000 0.000  0.000 0.015 0.015 0.061 0.121
# 25%   0.621 0.279 0.170 0.190 0.157  0.169  0.126 0.122  0.134 0.261 0.112 0.316 0.614
# 50%   2.008 0.400 0.254 0.339 0.234  0.365  0.267 0.231  0.281 0.394 0.182 0.425 0.779
# 75%   2.937 0.790 0.425 0.835 0.340  0.928  0.741 0.872  0.673 0.937 0.300 0.740 1.206
# 100% 11.377 8.651 7.822 6.075 3.018 11.487 11.420 5.916 11.432 9.406 7.009 2.497 7.562
# A higher median score than typical; but also not _as_ high as we often see
#     for more convincing doublets... maybe just flag for now


## First make some merges based on deep dives with markers ===
# Excit_A + Excit_E (co-express most of their top markers)
sce.ls$cellType.final[sce.ls$cellType.final %in% c("Excit_A", "Excit_E")] <- "Excit_A"
# Same with Excit_F + Excit_G
sce.ls$cellType.final[sce.ls$cellType.final %in% c("Excit_F", "Excit_G")] <- "Excit_F"


# Re-assignment of Inhib_D subclusters:
# - Keep '2' as Inhib_D, re-assign the rest
sce.ls$cellType.final[sce.ls$cellType.final == "1"] <- "drop.likelyDoublet"
sce.ls$cellType.final[sce.ls$cellType.final == "2"] <- "Inhib_D"
# and the rest - already have an _E : _H -->
sce.ls$cellType.final[sce.ls$cellType.final == "3"] <- "Inhib_I"
sce.ls$cellType.final[sce.ls$cellType.final == "4"] <- "Inhib_J"
sce.ls$cellType.final[sce.ls$cellType.final == "5"] <- "Inhib_K"
sce.ls$cellType.final[sce.ls$cellType.final == "6"] <- "Inhib_L"
sce.ls$cellType.final[sce.ls$cellType.final == "7"] <- "Inhib_M"
sce.ls$cellType.final[sce.ls$cellType.final == "8"] <- "Inhib_N"
sce.ls$cellType.final[sce.ls$cellType.final == "9"] <- "Excit_E"
# A 'hidden' excit - Re-assign as 'Excit_E' bc merged that with _A, above
sce.ls$cellType.final[sce.ls$cellType.final == "10"] <- "Inhib_O"
sce.ls$cellType.final[sce.ls$cellType.final == "11"] <- "Inhib_P"
sce.ls$cellType.final[sce.ls$cellType.final == "12"] <- "Inhib_Q"
sce.ls$cellType.final[sce.ls$cellType.final == "13"] <- "Inhib_R"

# Re-annotate 'Neuron.mixed_A' - only 'mixed' Neuronal population
sce.ls$cellType.final[sce.ls$cellType.final == "Neuron.mixed_A"] <- "Neuron.mixed"



# Re-factor
sce.ls$cellType.final <- factor(sce.ls$cellType.final)

# Check
plotReducedDim(sce.ls, "UMAP",
               colour_by = "cellType.final", point_size = 1.5,
               text_by = "cellType.final", text_size = 3
) +
  scale_color_manual(
    values = cell_colors.ls,
    labels = paste0(
      levels(sce.ls$cellType.final), " (",
      table(sce.ls$cellType.final), ")"
    )
  )


## Can go ahead and assign final cell class color (37 levels)
cell_colors.ls <- c(tableau20, tableau10medium, cbPalette[-6])
names(cell_colors.ls) <- levels(sce.ls$cellType.final)

# Save
save(sce.ls, annotationTab.ls, cell_colors.ls,
     file = here("snRNAseq_mouse", "processed_data", "SCE", "sce_updated_LS.rda")
)





### FINAL addition: (sub)region-specific annotations; re-annotation of 'Ependymal' ===========
#   - New insights from (sub)region localization of top marker genes to add to $cellType.final
#   - Additionally, Aqp4.Rbmps top markers --> ependymal, whereas the formerly-
#     called 'Ependymal' is actually choroid plexus

load(here("snRNAseq_mouse", "processed_data", "SCE", "sce_updated_LS.rda"), verbose = T)
# sce.ls, annotationTab.ls, cell_colors.ls

# Call the $cellType.final --> $cellType.exp[anded]
sce.ls$cellType.exp <- sce.ls$cellType.final
sce.ls$cellType.final <- NULL

## Expand annotationTab.ls to capture summary of data analysis
annotation.exp <- data.frame(matrix(unlist(rep(
  annotationTab.ls[which(annotationTab.ls$cellType == "Inhib_D"), ],
  12
)),
ncol = 2, byrow = T
))
colnames(annotation.exp) <- colnames(annotationTab.ls)

annotation.exp$cellType.exp <- c(
  "drop.likelyDoublet", "Excit_E",
  paste0("Inhib_", c("I", "J", "K", "L", "M", "N", "O", "P", "Q", "R"))
)

## And as before:
annotationTab.ls$cellType.exp <- annotationTab.ls$cellType
annotationTab.ls$cellType.exp[grep("Aqp4.Rbpms", annotationTab.ls$cellType.exp)] <- "Aqp4.Rbpms"
annotationTab.ls$cellType.exp[grep("Mural", annotationTab.ls$cellType.exp)] <- "Mural"
annotationTab.ls$cellType.exp[grep("Oligo", annotationTab.ls$cellType.exp)] <- "Oligo"
annotationTab.ls$cellType.exp[grep("Astro", annotationTab.ls$cellType.exp)] <- "Astro"
annotationTab.ls$cellType.exp[grep("Endo", annotationTab.ls$cellType.exp)] <- "Endo"

# Excit_A + Excit_E (co-express most of their top markers)
annotationTab.ls$cellType.exp[annotationTab.ls$cellType.exp %in% c("Excit_A", "Excit_E")] <- "Excit_A"
# Same with Excit_F + Excit_G
annotationTab.ls$cellType.exp[annotationTab.ls$cellType.exp %in% c("Excit_F", "Excit_G")] <- "Excit_F"

annotationTab.ls$cellType.exp[annotationTab.ls$cellType.exp == "Neuron.mixed_A"] <- "Neuron.mixed"

# Add expansion
annotationTab.ls <- rbind(annotationTab.ls, annotation.exp)


## Now finally add the new insights & (sub)region-specific prefixes, based on the above ===
annotationTab.ls$cellType.final <- annotationTab.ls$cellType.exp
annotationTab.ls$cellType.final[grep("Ependymal", annotationTab.ls$cellType.final)] <- "ChP"
annotationTab.ls$cellType.final[grep("Aqp4.Rbpms", annotationTab.ls$cellType.final)] <- "Ependymal"

annotationTab.ls$cellType.final[grep("Excit_A", annotationTab.ls$cellType.final)] <- "TNoS_Ex.A"
annotationTab.ls$cellType.final[grep("Excit_B", annotationTab.ls$cellType.final)] <- "Thal_Ex.B"
annotationTab.ls$cellType.final[grep("Excit_C", annotationTab.ls$cellType.final)] <- "TT.IG.SH_Ex.C"
annotationTab.ls$cellType.final[grep("Excit_D", annotationTab.ls$cellType.final)] <- "Chol_Ex.D"
annotationTab.ls$cellType.final[grep("Excit_E", annotationTab.ls$cellType.final)] <- "TT.IG.SH_Ex.E"
annotationTab.ls$cellType.final[grep("Excit_F", annotationTab.ls$cellType.final)] <- "TT.IG.SH_Ex.F"

annotationTab.ls$cellType.final[grep("Inhib_A", annotationTab.ls$cellType.final)] <- "Str_In.A"
annotationTab.ls$cellType.final[grep("Inhib_B", annotationTab.ls$cellType.final)] <- "Ventr_In.B"
annotationTab.ls$cellType.final[grep("Inhib_C", annotationTab.ls$cellType.final)] <- "LS_In.C"
annotationTab.ls$cellType.final[grep("Inhib_D", annotationTab.ls$cellType.final)] <- "LS_In.D"
annotationTab.ls$cellType.final[grep("Inhib_E", annotationTab.ls$cellType.final)] <- "IoC_In.E"
annotationTab.ls$cellType.final[grep("Inhib_F", annotationTab.ls$cellType.final)] <- "Str_In.F"
annotationTab.ls$cellType.final[grep("Inhib_G", annotationTab.ls$cellType.final)] <- "Sept_In.G"
annotationTab.ls$cellType.final[grep("Inhib_H", annotationTab.ls$cellType.final)] <- "Str_In.H"
annotationTab.ls$cellType.final[grep("Inhib_I", annotationTab.ls$cellType.final)] <- "Sept_In.I"

annotationTab.ls$cellType.final[grep("Inhib_J", annotationTab.ls$cellType.final)] <- "MS_In.J"
annotationTab.ls$cellType.final[grep("Inhib_K", annotationTab.ls$cellType.final)] <- "MS_In.K"
annotationTab.ls$cellType.final[grep("Inhib_L", annotationTab.ls$cellType.final)] <- "Str_In.L"
annotationTab.ls$cellType.final[grep("Inhib_M", annotationTab.ls$cellType.final)] <- "LS_In.M"
annotationTab.ls$cellType.final[grep("Inhib_N", annotationTab.ls$cellType.final)] <- "LS_In.N"
annotationTab.ls$cellType.final[grep("Inhib_O", annotationTab.ls$cellType.final)] <- "LS_In.O"
annotationTab.ls$cellType.final[grep("Inhib_P", annotationTab.ls$cellType.final)] <- "LS_In.P"
annotationTab.ls$cellType.final[grep("Inhib_Q", annotationTab.ls$cellType.final)] <- "LS_In.Q"
annotationTab.ls$cellType.final[grep("Inhib_R", annotationTab.ls$cellType.final)] <- "LS_In.R"

# Then re-annotate
sce.ls$cellType.final <- annotationTab.ls$cellType.final[match(
  sce.ls$cellType.exp,
  annotationTab.ls$cellType.exp
)]
sce.ls$cellType.final <- factor(sce.ls$cellType.final)
table(sce.ls$cellType.final)

## Re-assign final cell class color (37 levels)
cell_colors.ls <- c(tableau20, tableau10medium, cbPalette[-6])
names(cell_colors.ls) <- levels(sce.ls$cellType.final)


# Save
save(sce.ls, annotationTab.ls, cell_colors.ls,
     file = here("snRNAseq_mouse", "processed_data", "SCE", "sce_updated_LS.rda")
)


## 17Nov2022 final final re-annotations:
sce.ls$cellType.final <- as.character(sce.ls$cellType.final)
sce.ls$cellType.final[sce.ls$cellType.final=="Ventr_In.B"] <- "Neuroblast"
sce.ls$cellType.final[sce.ls$cellType.final=="Str_In.L"] <- "LS_In.L"
sce.ls$cellType.final <- factor(sce.ls$cellType.final)

# Re-assign the colors too
names(cell_colors.ls)[names(cell_colors.ls)=="Ventr_In.B"] <- "Neuroblast"
names(cell_colors.ls)[names(cell_colors.ls)=="Str_In.L"] <- "LS_In.L"

save(sce.ls, annotationTab.ls, cell_colors.ls,
     file="/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/processed_data/SCE/sce_updated_LS.rda")


# Re-compute cluster modularity 
load("/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/processed_data/SCE/graph_clusters_glmpca_LS-n3.rda", verbose=T)
mod.ratio.merged.HC <- pairwiseModularity(
  graph = snn.gr.glmpca,
  clusters = sce.ls$cellType.final,
  as.ratio = TRUE
)

# Heatmap
pdf("/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/plots/clusterModRatio_33clusters_Nov2022.pdf")
pheatmap(log2(mod.ratio.merged.HC + 1),
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("white", "blue"))(100),
         main = "Modularity ratio for 37 (including drops) graph-based clusters in LS (n=4)",
         fontsize_row = 7.5, fontsize_col = 7.5, angle_col = 90,
         display_numbers = T, number_format = "%.1f", fontsize_number = 5.5,
         na_col = "darkgrey"
)
grid::grid.text(label = "log2(ratio)", x = 0.97, y = 0.64, gp = grid::gpar(fontsize = 7))
dev.off()




## One more iteration Jan2023 - binned by region-specific(-ish) annotation ===
# (load SCE containing .rda)

sce.ls$cellType.broad <- factor(ss(as.character(sce.ls$cellType.final), "_", 1))

# Add this to the annotationTab 'summary'
annotationTab.ls$cellType.broad <- ss(annotationTab.ls$cellType.final, "_", 1)


save(sce.ls, annotationTab.ls, cell_colors.ls,
     file="/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/processed_data/SCE/sce_updated_LS.rda")














### For iSEE shiny app =============
load(here("snRNAseq_mouse", "processed_data", "SCE", "sce_updated_LS.rda"), verbose = T)

# As-is size:
lobstr::obj_size(sce.ls) / 1024^3 # 7.17 GB

# Remove the clusters we're not reporting
sce.ls <- sce.ls[, -grep("drop.", sce.ls$cellType.final)]
sce.ls <- sce.ls[, -grep("Neuron.mixed", sce.ls$cellType.final)]
sce.ls$cellType.final <- droplevels(sce.ls$cellType.final)

# Corresponding hexadecimal colors
cell_colors.ls <- cell_colors.ls[-grep("drop", names(cell_colors.ls))]
cell_colors.ls <- cell_colors.ls[-grep("Neuron.mixed", names(cell_colors.ls))]

sce.hold <- sce.ls

# Remove 'counts' & 'binomial_pearson_residuals' (just keep 'logcounts')
assay(sce.ls, "counts") <- NULL
assay(sce.ls, "binomial_pearson_residuals") <- NULL
lobstr::obj_size(sce.ls) / 1024^3 # 0.82 GB

# 'GLMPCA_50' is just the top 50 PCs from the 100 computed in 'GLMPCA_approx'
reducedDim(sce.ls, "GLMPCA_approx") <- NULL
lobstr::obj_size(sce.ls) / 1024^3 # 0.80 GB

# Because rownames are currently in Ensembl ID, re-assign & 'uniquify' to gene symbols
rownames(sce.ls) <- uniquifyFeatureNames(rowData(sce.ls)$gene_id, rowData(sce.ls)$gene_name)

lobstr::obj_size(sce.ls) / 1024^3

## Some add'l metrics for rowData ===
#    As created/implemented in https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/master/shiny_apps/00_clean_functions.R
rowData(sce.ls)$propNucleiExprs <- apply(
  assay(sce.hold, "counts"),
  1,
  function(x) {
    round(mean(x != 0), 3)
  }
)

# The above, by cell type ===
cellType.idx <- splitit(sce.ls$cellType.final)
rowdat.sce <- rowData(sce.ls)
for (i in names(cellType.idx)) {
  message(Sys.time(), " computing propNucleiExprs for ", i)
  rowdat.sce[, paste0("propExprsIn.", i)] <- apply(
    assay(sce.hold, "counts")[, cellType.idx[[i]]],
    1,
    function(x) {
      round(mean(x != 0), 3)
    }
  )
}
rowData(sce.ls) <- rowdat.sce

## Clean up unnecessary info
keepCols <- setdiff(
  colnames(colData(sce.ls)),
  c("high.mito", "sizeFactor", "clusters.glmpca", "cellType", "cellType.exp")
)
colData(sce.ls) <- colData(sce.ls)[, keepCols]

lobstr::obj_size(sce.ls) / 1024^3 # 0.82 GB

# Re-name these
sce.ls.small <- sce.ls
cell_cols.clean <- cell_colors.ls

## Save this
save(sce.ls.small, cell_cols.clean,
     file = here("snRNAseq_mouse", "processed_data", "SCE", "sce_for_iSEE_LS.rda")
)


##
pdf("/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/plots/FIX_THIS.pdf",height = 7, width = 10)

plotExpressionCustom(
  sce = sce.ls,
  exprs_values = "logcounts",
  # Ependymal cell markers
  features = c("Drd3", "Htr4"), #edit this with markers
  features_name = "custom-selected", 
  anno_name = "cellType.final",
  ncol = 3, point_alpha = 0.4, point_size = 0.9,
  scales = "free_y", swap_rownames = "gene_name"
) +
  scale_color_manual(values = c(cbPalette, tableau20, tableau10medium)) +
  ggtitle(label = paste0("Lionel's custom markers of interest")) +
  theme(
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7)
  )
dev.off()


##modifying script above for data subsetted for neurons or LS only clusters
pdf("/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/plots/Monoamine Receptors Example Expression in LS Clusters.pdf",height = 3, width = 6)

plotExpressionCustom(
  sce = ls.sce.ls,
  exprs_values = "logcounts",
  # 
  features = c("Htr1f", "Htr2c", "Htr4", "Htr7", "Adra1a", "Adra1b"), #edit this with markers
  features_name = "Region Markers", 
  anno_name = "cellType.final",
  ncol = 3, point_alpha = 0.4, point_size = 0.9,
  scales = "free_y", swap_rownames = "gene_name"
) +
  scale_color_manual(values = c(cbPalette, tableau20, tableau10medium)) +
  ggtitle(label = paste0("Region Markers")) +
  theme(
    plot.title = element_text(size = 12),
    axis.text.x = element_text(size = 7)
  )
dev.off()


## Reproducibility information ====
print("Reproducibility information:")
Sys.time()
# [1] "2022-07-01 01:05:06 EDT"
proc.time()
#     user    system   elapsed
#   1497.611  135.700 8502.373
options(width = 120)
session_info()
# ─ Session info ───────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.1.2 Patched (2021-11-04 r81138)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2022-07-01
# pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────
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
# bluster              * 1.4.0    2021-10-26 [2] Bioconductor
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
# DropletUtils         * 1.14.2   2022-01-09 [2] Bioconductor
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
# HDF5Array              1.22.1   2021-11-14 [2] Bioconductor
# here                 * 1.0.1    2020-12-13 [2] CRAN (R 4.1.2)
# igraph                 1.3.2    2022-06-13 [2] CRAN (R 4.1.2)
# IRanges              * 2.28.0   2021-10-26 [2] Bioconductor
# irlba                  2.3.5    2021-12-06 [2] CRAN (R 4.1.2)
# jaffelab             * 0.99.31  2021-12-13 [1] Github (LieberInstitute/jaffelab@2cbd55a)
# labeling               0.4.2    2020-10-20 [2] CRAN (R 4.1.0)
# lattice                0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
# lifecycle              1.0.1    2021-09-24 [2] CRAN (R 4.1.2)
# limma                  3.50.3   2022-04-07 [2] Bioconductor
# lobstr                 1.1.2    2022-06-22 [2] CRAN (R 4.1.2)
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
# prettyunits            1.1.1    2020-01-24 [2] CRAN (R 4.1.0)
# purrr                  0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
# R.methodsS3            1.8.2    2022-06-13 [2] CRAN (R 4.1.2)
# R.oo                   1.25.0   2022-06-12 [2] CRAN (R 4.1.2)
# R.utils                2.12.0   2022-06-28 [2] CRAN (R 4.1.2)
# R6                     2.5.1    2021-08-19 [2] CRAN (R 4.1.2)
# rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.1.2)
# RColorBrewer           1.1-3    2022-04-03 [2] CRAN (R 4.1.2)
# Rcpp                   1.0.8.3  2022-03-17 [2] CRAN (R 4.1.2)
# RCurl                  1.98-1.7 2022-06-09 [2] CRAN (R 4.1.2)
# ResidualMatrix         1.4.0    2021-10-26 [1] Bioconductor
# rhdf5                  2.38.1   2022-03-10 [2] Bioconductor
# rhdf5filters           1.6.0    2021-10-26 [2] Bioconductor
# Rhdf5lib               1.16.0   2021-10-26 [2] Bioconductor
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
# ──────────────────────────────────────────────────────────────────────────────────────────



