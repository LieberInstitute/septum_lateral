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

source("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/plotExpressionCustom.R")


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
        #                          batch=as.factor(sce.lc$Sample),
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


snn.gr.glmpcamnn <- buildSNNGraph(sce.ls, k=20, use.dimred="glmpca_mnn_50")
clusters.glmpcamnn <- igraph::cluster_walktrap(snn.gr.glmpcamnn)
table(clusters.glmpcamnn$membership)
        #



# Compute logcounts to visualize expression in the traditional way
sce.test <- multiBatchNorm(sce.test, batch=sce.test$Sample)

