### mouse LS snRNA-seq analysis
### Miscellaneous/manual exploration of data
### qsub -l bluejay,mf=60G,h_vmem=64G,h_fsize=40G
### Initiated LAR,MNT 03May2022
  #   - moved from '05_enrichmentTests_TrkB_others.R' on 01Jul2022

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

# ===



load(here("snRNAseq_mouse", "processed_data","SCE", "sce_updated_LS.rda"), verbose=T)
# sce.ls, annotationTab.ls, cell_colors.ls

# For ease of querying
rownames(sce.ls) <- uniquifyFeatureNames(rowData(sce.ls)$gene_id, rowData(sce.ls)$gene_name)


## Combinatorial expression ======
geneExprs <- assay(sce.ls, "counts")

lapply(cellType.idx, function(x){
  signif(prop.table(table(geneExprs["Slc17a7",x] > 0 &
                            geneExprs["Slc17a6",x] > 0)
  ),2)
})

lapply(cellType.idx, function(x){
  signif(prop.table(table(geneExprs["Slc17a7",x] > 0 &
                            geneExprs["Gad1",x] > 0)
  ),2)
})


lapply(cellType.idx, function(x){
  signif(prop.table(table(geneExprs["Slc17a6",x] > 0 &
                            geneExprs["Gad1",x] > 0)
  ),2)
})


lapply(cellType.idx, function(x){
  signif(prop.table(table(geneExprs["Slc17a7",x] > 0 &
                            geneExprs["Gad2",x] > 0)
  ),2)
})

lapply(cellType.idx, function(x){
  signif(prop.table(table(geneExprs["Slc17a6",x] > 0 &
                            geneExprs["Gad2",x] > 0)
  ),2)
})

lapply(cellType.idx, function(x){
  signif(prop.table(table(geneExprs["Gad1",x] > 0 &
                            geneExprs["Gad2",x] > 0)
  ),2)
})







sce.ls$MSN.D1 <- NA
sce.ls$MSN.D1[grep("MSN.D1", sce.ls$cellType)] <- TRUE

sce.ls$MSN.D1.0 <- NA
sce.ls$MSN.D1.0 <- sce.ls$MSN.D1 & geneExprs["DRD1", ]==0

# Now plot & color by this - is DRD2 expressed?
plotExpression(sce.ls, exprs_values = "logcounts", features=c("DRD1", "DRD2"),
               x="cellType", colour_by="MSN.D1.0", point_alpha=0.5, point_size=1.2, ncol=5,
               add_legend=T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25))


sce.ls$MSN.D2 <- NA
sce.ls$MSN.D2[grep("MSN.D2", sce.ls$cellType)] <- TRUE

sce.ls$MSN.D2.0 <- NA
sce.ls$MSN.D2.0 <- sce.ls$MSN.D2 & geneExprs["DRD2", ]==0

# Now plot & color by this - is DRD1 expressed?
plotExpression(sce.ls, exprs_values = "logcounts", features=c("DRD1", "DRD2"),
               x="cellType", colour_by="MSN.D2.0", point_alpha=0.5, point_size=1.2, ncol=5,
               add_legend=T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25))




# non-0 expression of BOTH? ===
table(geneExprs["DRD1", union(which(sce.ls$MSN.D1), which(sce.ls$MSN.D2))] > 0 &
        geneExprs["DRD2", union(which(sce.ls$MSN.D1), which(sce.ls$MSN.D2))] > 0)






## Reproducibility information ====
print('Reproducibility information:')
Sys.time()
#[1] "2022-05-10 15:48:33 EDT"
proc.time()
#     user    system   elapsed 
#    148.449    5.915 5428.779 
options(width = 120)
session_info()
#


