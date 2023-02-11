### Lionel learning to plot UMAPs
###    qrsh -l bluejay,mf=60G,h_vmem=64G,h_fsize=40G
### Initiated LAR,ST 06Jan2023

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
library(bluster)

here()

#If you need to mamually set working directory
setwd("/dcs04/lieber/marmaypag/pilotLS_LIBD1070/")
getwd() #to check it's in the right place


# plotExpressionCustom for nicer aesthetics
source("/dcs04/lieber/lcolladotor/pilotLC_LIBD001/locus-c/code/analyses_sn/plotExpressionCustom.R")


## ===

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


cell_colors.ls <- c(tableau20, tableau10medium, cbPalette[-6])

load("/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/processed_data/SCE/sce_updated_LS.rda", verbose = T)
# sce.ls, annotationTab.ls, cell_colors.ls

#Check to see it's the right processed data
rowData(sce.ls)
colData(sce.ls)

#Drop the 2 populations we don't want shown
sce.ls <- sce.ls[, -grep("drop", sce.ls$cellType.final)]
sce.ls <- sce.ls[, -grep("mixed", sce.ls$cellType.final)]
sce.ls$cellType.final <- droplevels(sce.ls$cellType.final)

#to check if it dropped
levels(sce.ls$cellType.final)

#change Ens ID to gene names
rownames(sce.ls) <- uniquifyFeatureNames(rowData(sce.ls)$gene_id, rowData(sce.ls)$gene_name)


pdf("/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/plots/UMAP_Trpc4.pdf")
plotUMAP(sce.ls, colour_by="Trpc4")
dev.off()


pdf("/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/plots/UMAP_Snap25.pdf")
plotUMAP(sce.ls, colour_by="Snap25")
dev.off()





#Can I do multiple? Let's see...

#Update, it works!

pdf("/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/plots/Expr Cell Type Markers on UMAP.pdf", width = 9)

plotUMAP(sce.ls, colour_by = "cellType.final"
) +
  scale_color_manual(
    values = cell_colors.ls,
    labels = paste0(
      levels(sce.ls$cellType.final), " (",
      table(sce.ls$cellType.final), ")"
    )
  )

plotUMAP(sce.ls, colour_by="Snap25")

plotUMAP(sce.ls, colour_by="Gad1")

plotUMAP(sce.ls, colour_by="Gad2")

plotUMAP(sce.ls, colour_by="Slc17a7")

plotUMAP(sce.ls, colour_by="Slc17a6")

plotUMAP(sce.ls, colour_by="Mbp")

plotUMAP(sce.ls, colour_by="Slc1a2")

dev.off()

pdf("/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/plots/UMAP_Practice.pdf")
plotUMAP(sce.ls, colour_by = "cellType.final"
) +
  scale_color_manual(
    values = cell_colors.ls,
    labels = paste0(
      levels(sce.ls$cellType.final), " (",
      table(sce.ls$cellType.final), ")"
    )
  )
dev.off()

