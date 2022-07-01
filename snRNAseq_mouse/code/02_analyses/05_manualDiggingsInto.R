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


## Want to know %(+) per cell type for Slc17a7, Slca7a6, Gad1 & Gad2 %'s =====
# Adapted from https://github.com/LieberInstitute/10xPilot_snRNAseq-human/blob/master/shiny_apps/00_clean_functions.R
sce.hold <- sce.ls

cellType.idx <- splitit(sce.ls$cellType.final)
rowdat.sce <- rowData(sce.ls)
for(i in names(cellType.idx)){
  message(Sys.time(), " computing propNucleiExprs for ", i)
  rowdat.sce[, paste0("propExprsIn.", i)] <- apply(
    assay(sce.ls, "counts")[, cellType.idx[[i]]],
    1,
    function(x){
      round(mean(x != 0), 3)
    }
  )
}
rowData(sce.ls) <- rowdat.sce
# For ease of querying
rownames(sce.ls) <- uniquifyFeatureNames(rowData(sce.ls)$gene_id, rowData(sce.ls)$gene_name)

rowData(sce.ls)[c("Slc17a7", "Slc17a6", "Gad1", "Gad2"), ]

# Save this for reference
rowdat.sce <- rowData(sce.ls)

write.table(rowdat.sce, file=here("snRNAseq_mouse","processed_data","tables",
                                  #"rowData_with_propExprsd_byCellType.tsv"),
                                  "rowData_with_propExprsd_byCellType.final.tsv"),
            sep="\t", quote = FALSE)




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


