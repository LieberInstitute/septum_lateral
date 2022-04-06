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

here()
    #[1] "/dcs04/lieber/marmaypag/pilotLS_LIBD1070"

source("/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/analyses_sn/plotExpressionCustom.R")


### Feature selection with deviance residuals =====

# Load working SCE
load(here("snRNAseq_mouse", "processed_data","SCE", "sce_working_LS.rda"), verbose=T)
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
                                   assay="counts", fam="binomial", sorted=F,
                                   # these are default params btw
                                   batch=as.factor(sce.ls$Sample))


# Btw:
table(is.na(rowData(sce.ls)$binomial_deviance))
    # FALSE  TRUE 
    # 29556  7045

# Observe:
pdf(here("snRNAseq_mouse","plots","featureSelxn_binomialDeviance-byGene.pdf"), height=5)
plot(sort(rowData(sce.ls)$binomial_deviance, decreasing=T),
     type="l", xlab="ranked genes",
     ylab="binomial deviance", main="Feature Selection with Deviance: LS (n=4)")
abline(v=2000, lty=2, col="red")
dev.off()


## Approximate GLM-PCA with PCA on null deviance residuals
Sys.time()
    #[1] "2022-03-23 10:37:02 EDT"
sce.ls <- nullResiduals(sce.ls, assay="counts", fam="binomial",  # default params
                        #type="deviance")#, batch=as.factor(sce.ls$Sample))
                        type="pearson")
            # MNT comment: previously in other projects, 'batch=' threw an error
            #             (and we perform MNN-batch correction anyway, if relevant)

            #           - Using Pearson residuals bc deviance took nearly an hour
            #             with up to 96G RAM...
Sys.time()
    #[1] "2022-03-23 11:52:14 EDT"

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
     file=here("snRNAseq_mouse", "processed_data","SCE", "sce_updated_LS.rda"))



## Feature selection: Take top 2000 highly-deviant genes ('HDG's) for PCA
 #      (based on total binomial deviance, computed above)
hdgs.ls <- rownames(sce.ls)[order(rowData(sce.ls)$binomial_deviance, decreasing=T)][1:2000]
hdgs.symbols <- rowData(sce.ls)$gene_name[match(hdgs.ls, rowData(sce.ls)$gene_id)]
# Out of curiosity
c("Snap25", "Mbp", "Gad1") %in% hdgs.symbols   # all TRUE

# Run PCA:
Sys.time()
    #[1] "2022-03-23 12:44:41 EDT"
set.seed(109)
sce.ls <-  runPCA(sce.ls,
                  #exprs_values="binomial_deviance_residuals",
                  exprs_values="binomial_pearson_residuals",
                  subset_row=hdgs.ls, ncomponents=100,
                  name="GLMPCA_approx",
                  BSPARAM=BiocSingular::IrlbaParam())
Sys.time()
    #[1] "2022-03-23 14:22:15 EDT"      - this took ~1hr40min to compute...



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
sce.ls <- runUMAP(sce.ls, dimred="GLMPCA_approx",
                  n_dimred=50, name="UMAP")
# t-SNE
set.seed(109)
sce.ls <- runTSNE(sce.ls, dimred="GLMPCA_approx",
                  n_dimred=50, name="TSNE")

# Visualize top PCs and these 2D embeddings:
pdf(here("snRNAseq_mouse","plots","reducedDims_mouseLS-n4_noBatchCorrxn.pdf"))
plotReducedDim(sce.ls, dimred="GLMPCA_approx", colour_by="Sample",
               ncomponents=4, point_alpha=0.3, point_size=1.5)
# UMAPs, along with some other metrics
plotReducedDim(sce.ls, dimred="UMAP", colour_by="Sample",
               point_alpha=0.3, point_size=1.5)
plotReducedDim(sce.ls, dimred="UMAP", colour_by="sum",
               point_alpha=0.3, point_size=1.5)
plotReducedDim(sce.ls, dimred="UMAP", colour_by="doubletScore",
               point_alpha=0.3, point_size=1.5)
# TSNE
plotReducedDim(sce.ls, dimred="TSNE", colour_by="Sample",
               point_alpha=0.3, point_size=1.5)
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
     file=here("snRNAseq_mouse", "processed_data","SCE", "sce_updated_LS.rda"))



### Clustering ====================
# Perform graph-based clustering, as in Tran-Maynard, et al. Neuron 2021

reducedDim(sce.ls, "GLMPCA_50") <- reducedDim(sce.ls, "GLMPCA_approx")[ ,1:50]

snn.gr.glmpca <- buildSNNGraph(sce.ls, k=20, use.dimred="GLMPCA_50")
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
     file=here("snRNAseq_mouse", "processed_data","SCE", "graph_clusters_glmpca_LS-n3.rda"))

table(sce.ls$clusters.glmpca, sce.ls$Sample)
    # pretty even distribution across samples

## doubletScore distributions / cluster?
cellClust.idx <- splitit(sce.ls$clusters.glmpca)
sapply(cellClust.idx, function(x){round(quantile(sce.ls$doubletScore[x]), 2)})
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
sce.ls <- multiBatchNorm(sce.ls, batch=sce.ls$Sample)

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
     file=here("snRNAseq_mouse", "processed_data","SCE", "sce_updated_LS.rda"))


## Broad markers of interest:
markers.mathys.tran = list(
    'neuron' = c('SYT1', 'SNAP25', 'GRIN1'),
    'excit_neuron' = c('CAMK2A', 'NRGN','SLC17A7', 'SLC17A6', 'SLC17A8'),
    'inhib_neuron' = c('GAD1', 'GAD2', 'SLC32A1'),
    # Norepinephrine & serotonergic markers
    'neuron.NE' = c("TH", "DBH", "SLC6A2", "SLC18A2", "GCH1", "DDC"), #SLC6A3 - saw no DAT
    'neuron.5HT' = c("SLC6A4", "TPH1", "TPH2", "DDC"),
    # SERT, serotonin T (aka 5-HTT); 
    'monoamine.metab' = c("COMT", "MAOA", "MAOB"),
    # MSN markers
    'MSNs.pan' = c("PPP1R1B","BCL11B"),# "CTIP2")
    'MSNs.D1' = c("DRD1", "PDYN", "TAC1"),
    'MSNs.D2' = c("DRD2", "PENK"),
    ## Non-neuronal:
    'oligodendrocyte' = c('MBP', 'MOBP', 'PLP1'),
    'oligo_precursor' = c('PDGFRA', 'VCAN', 'CSPG4'),
    'microglia' = c('CD74', 'CSF1R', 'C3'),
    'astrocyte' = c('GFAP', 'TNC', 'AQP4', 'SLC1A2'),
    'endothelial' = c('CLDN5', 'FLT1', 'VTN'),
    # Post-hoc from Tran-Maynard, et al. Neuron 2021
    'differn_committed_OPC' = c("SOX4", "BCAN", "GPR17", "TNS3"),
    'Tcell' = c('SKAP1', 'ITK', 'CD247'),
    'Mural' = c('COL1A2', 'TBX18', 'RBPMS'),
    'Macro' = c('CD163', 'SIGLEC1', 'F13A1')
)

# Will have to 'make these mouse'
broadMarkers <- markers.mathys.tran
for(i in 1:length(broadMarkers)){
    broadMarkers[[i]] <- paste0(substr(broadMarkers[[i]], 1,1),
                                   tolower(substr(broadMarkers[[i]], 2, nchar(broadMarkers[[i]]))))
}

table(unname(unlist(broadMarkers)) %in% rowData(sce.ls)$gene_name) #all good

rownames(sce.ls) <- rowData(sce.ls)$gene_name

pdf(here("snRNAseq_mouse","plots",paste0("LS-n4_expression_broadMarkers_GLMPCA-graphClusters.pdf")),
    height=6, width=14)
for(i in 1:length(broadMarkers)){
    print(
        plotExpressionCustom(sce = sce.ls,
                             exprs_values = "logcounts",
                             features = broadMarkers[[i]], 
                             features_name = names(broadMarkers)[[i]], 
                             anno_name = "clusters.glmpca",
                             ncol=2, point_alpha=0.4, point_size=0.9,
                             scales="free_y") +  
            ggtitle(label=paste0("mouse LS (n4) clusters: ",
                                 names(broadMarkers)[[i]], " markers")) +
            theme(plot.title = element_text(size = 12),
                  axis.text.x = element_text(size=7))
    )
}
dev.off()



# for temp exploration:
# pdf(here("snRNAseq_mouse","plots",paste0("temp_markerExplore.pdf")),height=6, width=14)
#     print(
#         plotExpressionCustom(sce = sce.ls,
#                              exprs_values = "logcounts",
#                              features = c("Gfap", "Aldh1l1", "Aldoc", "S100b"),
#                              features_name = "", 
#                              anno_name = "clusters.glmpca",
#                              ncol=2, point_alpha=0.4, point_size=0.9,
#                              scales="free_y", swap_rownames="gene_name") +
#             theme(plot.title = element_text(size = 12),
#                   axis.text.x = element_text(size=7))
#     )
# dev.off()



## Annotate clusters ===
annotationTab.ls <- data.frame(cluster=c(1:35))
annotationTab.ls$cellType <- NA
annotationTab.ls$cellType[c(4,11,20,21,34)] <- paste0("Excit_", c("A","B","C","D", "E"))
annotationTab.ls$cellType[c(1,2,3,8, 12,15,23,25, 26,30)] <-
    paste0("Inhib_", c("A","B","C","D", "E","F","G","H", "I","J"))
annotationTab.ls$cellType[c(31)] <- paste0("Neuron.mixed_", c("A"))
annotationTab.ls$cellType[c(19)] <- "Neuron.Ppp1r1b"

annotationTab.ls$cellType[c(6,22)] <- paste0("Astro_", c("A","B"))
annotationTab.ls$cellType[c(13)] <- "Astro.OPC_COP"
annotationTab.ls$cellType[c(17,24)] <- paste0("Aqp4.Rbpms_", c("A","B"))
annotationTab.ls$cellType[c(16,27)] <- paste0("Endo_", c("A","B"))
annotationTab.ls$cellType[c(14,29)] <- paste0("Mural_", c("A","B"))
annotationTab.ls$cellType[c(9)] <- "Micro"
annotationTab.ls$cellType[c(7,18)] <- paste0("Oligo_", c("A","B"))
annotationTab.ls$cellType[10] <- c("OPC")
annotationTab.ls$cellType[33] <- c("OPC_COP")

annotationTab.ls$cellType[28] <- c("drop.lowNTx")
annotationTab.ls$cellType[32] <- c("drop.doublet")
annotationTab.ls$cellType[c(5,35)] <- paste0("ambig.glial_", c("A","B"))




sce.ls$cellType <- annotationTab.ls$cellType[match(sce.ls$clusters.glmpcamnn,
                                                   annotationTab.ls$cluster)]
sce.ls$cellType <- factor(sce.ls$cellType)

options(width=100)
table(sce.ls$cellType)




## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
    #[1] "2022-03-31 14:28:11 EDT"
proc.time()
    #     user    system   elapsed 
    #  446.909    16.777 15914.965 
options(width = 120)
session_info()
    #─ Session info ──────────────────────────────────────────────────────────────
    # setting  value
    # version  R version 4.1.2 Patched (2021-11-04 r81138)
    # os       CentOS Linux 7 (Core)
    # system   x86_64, linux-gnu
    # ui       X11
    # language (EN)
    # collate  en_US.UTF-8
    # ctype    en_US.UTF-8
    # tz       US/Eastern
    # date     2022-03-31
    # pandoc   2.13 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-4.1.x/bin/pandoc
    # 
    # ─ Packages ──────────────────────────────────────────────────────────────────
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
    # cli                    3.2.0    2022-02-14 [2] CRAN (R 4.1.2)
    # cluster                2.1.3    2022-03-28 [3] CRAN (R 4.1.2)
    # colorspace             2.0-3    2022-02-21 [2] CRAN (R 4.1.2)
    # cowplot                1.1.1    2020-12-30 [2] CRAN (R 4.1.2)
    # crayon                 1.5.1    2022-03-26 [2] CRAN (R 4.1.2)
    # DBI                    1.1.2    2021-12-20 [2] CRAN (R 4.1.2)
    # DelayedArray           0.20.0   2021-10-26 [2] Bioconductor
    # DelayedMatrixStats     1.16.0   2021-10-26 [2] Bioconductor
    # digest                 0.6.29   2021-12-01 [2] CRAN (R 4.1.2)
    # dplyr                  1.0.8    2022-02-08 [2] CRAN (R 4.1.2)
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
    # ggplot2              * 3.3.5    2021-06-25 [2] CRAN (R 4.1.0)
    # ggrepel                0.9.1    2021-01-15 [2] CRAN (R 4.1.0)
    # glue                   1.6.2    2022-02-24 [2] CRAN (R 4.1.2)
    # googledrive            2.0.0    2021-07-08 [2] CRAN (R 4.1.0)
    # gridExtra            * 2.3      2017-09-09 [2] CRAN (R 4.1.0)
    # gtable                 0.3.0    2019-03-25 [2] CRAN (R 4.1.0)
    # HDF5Array              1.22.1   2021-11-14 [2] Bioconductor
    # here                 * 1.0.1    2020-12-13 [2] CRAN (R 4.1.2)
    # igraph                 1.2.11   2022-01-04 [2] CRAN (R 4.1.2)
    # IRanges              * 2.28.0   2021-10-26 [2] Bioconductor
    # irlba                  2.3.5    2021-12-06 [2] CRAN (R 4.1.2)
    # jaffelab             * 0.99.31  2021-12-13 [1] Github (LieberInstitute/jaffelab@2cbd55a)
    # labeling               0.4.2    2020-10-20 [2] CRAN (R 4.1.0)
    # lattice                0.20-45  2021-09-22 [3] CRAN (R 4.1.2)
    # lifecycle              1.0.1    2021-09-24 [2] CRAN (R 4.1.2)
    # limma                  3.50.1   2022-02-17 [2] Bioconductor
    # locfit                 1.5-9.5  2022-03-03 [2] CRAN (R 4.1.2)
    # magrittr               2.0.3    2022-03-30 [2] CRAN (R 4.1.2)
    # Matrix                 1.4-1    2022-03-23 [3] CRAN (R 4.1.2)
    # MatrixGenerics       * 1.6.0    2021-10-26 [2] Bioconductor
    # matrixStats          * 0.61.0   2021-09-17 [2] CRAN (R 4.1.2)
    # metapod                1.2.0    2021-10-26 [2] Bioconductor
    # munsell                0.5.0    2018-06-12 [2] CRAN (R 4.1.0)
    # pillar                 1.7.0    2022-02-01 [2] CRAN (R 4.1.2)
    # pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 4.1.0)
    # purrr                  0.3.4    2020-04-17 [2] CRAN (R 4.1.0)
    # R.methodsS3            1.8.1    2020-08-26 [2] CRAN (R 4.1.0)
    # R.oo                   1.24.0   2020-08-26 [2] CRAN (R 4.1.0)
    # R.utils                2.11.0   2021-09-26 [2] CRAN (R 4.1.2)
    # R6                     2.5.1    2021-08-19 [2] CRAN (R 4.1.2)
    # rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 4.1.2)
    # RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 4.1.0)
    # Rcpp                   1.0.8.3  2022-03-17 [2] CRAN (R 4.1.2)
    # RCurl                  1.98-1.6 2022-02-08 [2] CRAN (R 4.1.2)
    # ResidualMatrix         1.4.0    2021-10-26 [1] Bioconductor
    # rhdf5                  2.38.1   2022-03-10 [2] Bioconductor
    # rhdf5filters           1.6.0    2021-10-26 [2] Bioconductor
    # Rhdf5lib               1.16.0   2021-10-26 [2] Bioconductor
    # rlang                  1.0.2    2022-03-04 [2] CRAN (R 4.1.2)
    # rprojroot              2.0.2    2020-11-15 [2] CRAN (R 4.1.0)
    # rsvd                   1.0.5    2021-04-16 [2] CRAN (R 4.1.2)
    # S4Vectors            * 0.32.4   2022-03-24 [2] Bioconductor
    # ScaledMatrix           1.2.0    2021-10-26 [2] Bioconductor
    # scales                 1.1.1    2020-05-11 [2] CRAN (R 4.1.0)
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
    # tibble                 3.1.6    2021-11-07 [2] CRAN (R 4.1.2)
    # tidyselect             1.1.2    2022-02-21 [2] CRAN (R 4.1.2)
    # utf8                   1.2.2    2021-07-24 [2] CRAN (R 4.1.0)
    # vctrs                  0.4.0    2022-03-30 [2] CRAN (R 4.1.2)
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
    # ─────────────────────────────────────────────────────────────────────────────

