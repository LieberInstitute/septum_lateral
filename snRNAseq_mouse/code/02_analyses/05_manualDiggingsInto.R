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


