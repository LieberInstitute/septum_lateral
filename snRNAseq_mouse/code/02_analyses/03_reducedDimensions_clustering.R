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


#pdf(here("snRNAseq_mouse","plots",paste0("LS-n4_expression_broadMarkers_GLMPCA-graphClusters.pdf")),
pdf(here("snRNAseq_mouse","plots",paste0("LS-n4_expression_broadMarkers_GLMPCA-graphClusters_annotated.pdf")),
    height=6, width=14)
for(i in 1:length(broadMarkers)){
    print(
        plotExpressionCustom(sce = sce.ls,
                             exprs_values = "logcounts",
                             features = broadMarkers[[i]], 
                             features_name = names(broadMarkers)[[i]], 
                             #anno_name = "clusters.glmpca",
                             anno_name = "cellType",
                             ncol=2, point_alpha=0.4, point_size=0.9,
                             scales="free_y", swap_rownames="gene_name") +  
            ggtitle(label=paste0("mouse LS (n4) clusters: ",
                                 names(broadMarkers)[[i]], " markers")) +
            theme(plot.title = element_text(size = 12),
                  axis.text.x = element_text(size=7))
    )
}
dev.off()



#for temp exploration:
# pdf(here("snRNAseq_mouse","plots",paste0("temp_markerExplore.pdf")),height=6, width=14)
# plotExpressionCustom(sce = sce.ls,
#                      exprs_values = "logcounts",
#                      # Ependymal cell markers
#                      features = c("Rarres2", "Ccdc153","Tmem212","S100b","Acta2","Foxj1"), 
#                      features_name = "custom-selected", 
#                      anno_name = "cellType",
#                      ncol=2, point_alpha=0.4, point_size=0.9,
#                      scales="free_y", swap_rownames="gene_name") +  
#     ggtitle(label=paste0("Lionel's custom markers of interest")) +
#     theme(plot.title = element_text(size = 12),
#           axis.text.x = element_text(size=7))
# dev.off()



## Annotate clusters ===
annotationTab.ls <- data.frame(cluster=c(1:35))
annotationTab.ls$cellType <- NA
annotationTab.ls$cellType[c(4,11,20,21,34,25, 26)] <- paste0("Excit_", c("A","B","C","D", "E","F","G"))
annotationTab.ls$cellType[c(1,2,3,8, 12,15,23,30)] <-
    paste0("Inhib_", c("A","B","C","D", "E","F","G","H"))
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


sce.ls$cellType <- annotationTab.ls$cellType[match(sce.ls$clusters.glmpca,
                                                   annotationTab.ls$cluster)]
sce.ls$cellType <- factor(sce.ls$cellType)

table(sce.ls$cellType)
# ambig.glial_A  ambig.glial_B   Aqp4.Rbpms_A   Aqp4.Rbpms_B        Astro_A 
#            63             76            106            146           4153 
#       Astro_B  Astro.OPC_COP   drop.doublet    drop.lowNTx         Endo_A 
#           171           1529            150            180            118 
#        Endo_B        Excit_A        Excit_B        Excit_C        Excit_D 
#            40            451            670            232             67 
#       Excit_E        Excit_F        Excit_G        Inhib_A        Inhib_B 
#            34            202            728            797            760 
#       Inhib_C        Inhib_D        Inhib_E        Inhib_F        Inhib_G 
#           191           4327            306           4364            363 
#       Inhib_H          Micro        Mural_A        Mural_B Neuron.mixed_A 
#            86            166            112             34             73 
#Neuron.Ppp1r1b        Oligo_A        Oligo_B            OPC        OPC_COP 
#            65           1478            105            464             53 


## Save
save(sce.ls,
     file=here("snRNAseq_mouse", "processed_data","SCE", "sce_updated_LS.rda"))



## re-print reducedDims with these new annotations ===
pdf(here("snRNAseq_mouse","plots","reducedDims_mouseLS-n4_graph-basedClusters_annotated.pdf"))
# UMAPs, along with some other metrics
plotReducedDim(sce.ls, dimred="UMAP", colour_by="Sample",
               point_alpha=0.3, point_size=1.5) +
    ggtitle("UMAP of LS (n=4)")
plotReducedDim(sce.ls, dimred="UMAP", colour_by="Sex",
               point_alpha=0.3, point_size=1.5) +
    ggtitle("UMAP of LS (n=4)")
plotReducedDim(sce.ls, dimred="UMAP", colour_by="sum",
               point_alpha=0.3, point_size=1.5) +
    ggtitle("UMAP of LS (n=4)")
plotReducedDim(sce.ls, dimred="UMAP", colour_by="doubletScore",
               point_alpha=0.3, point_size=1.5) +
    ggtitle("UMAP of LS (n=4)")
plotReducedDim(sce.ls, dimred="UMAP", colour_by="cellType",
               text_by="cellType", text_size=3,
               point_alpha=0.3, point_size=1.5) +
    guides(color=guide_legend(ncol=1)) +
    ggtitle("TSNE of LS (n=4), colored by annotated graph-based clusters")

# TSNE
plotReducedDim(sce.ls, dimred="TSNE", colour_by="cellType",
               text_by="cellType", text_size=3,
               point_alpha=0.3, point_size=1.5) +
    guides(color=guide_legend(ncol=1)) +
    ggtitle("UMAP of LS (n=4), colored by annotated graph-based clusters")
dev.off()


## Look at some distributions with these annotations
cellClust.idx <- splitit(sce.ls$cellType)
sapply(cellClust.idx, function(x){round(quantile(sce.ls$sum[x]), 2)})

sapply(cellClust.idx, function(x){round(quantile(sce.ls$doubletScore[x]), 2)})






## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
    #
proc.time()
    #     user    system   elapsed 
    #
options(width = 120)
session_info()
    #

