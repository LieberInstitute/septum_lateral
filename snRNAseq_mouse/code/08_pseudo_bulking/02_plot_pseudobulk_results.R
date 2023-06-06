library("here")
library("SingleCellExperiment")
library("scater")
library("sessioninfo")



################### Load sce object and pseudobulk results ####################

load(
    here(
        "snRNAseq_mouse",
        "processed_data",
        "SCE",
        "sce_updated_LS.rda"
    ),
    verbose = TRUE
)
# Loading objects:
#   sce.ls
#   annotationTab.ls
#   cell_colors.ls

load(
    here(
        "snRNAseq_mouse",
        "processed_data",
        "SCE",
        "sce_pseudobulking_LS_and_broad.rda"
    ),
    verbose = TRUE
)
# Loading objects:
#   sce_pseudo_all
#   sce_pseudo_LS
#   sce_pseudo_LS.Sept
#   sce_pseudo_neuronal
#   sce_pseudo_neuronal.broad

###############################################################################



############################## Setting up colors ##############################

cell_colors_all <- c(1:length(levels(colData(sce_pseudo_all)$cellType.broad)))
names(cell_colors_all) <- levels(colData(sce_pseudo_all)$cellType.broad)
neuronal_names <- c("Chol", "LS", "Sept", "Str", "MS", "TNoS", "TT.IG.SH", "Thal", "IoC")
cell_colors_all[names(cell_colors_all) %in% neuronal_names] <- "#e68f00"
cell_colors_all[!names(cell_colors_all) %in% neuronal_names] <- "#006164"

cell_colors_LS <- cell_colors.ls[grep("LS",names(cell_colors.ls))]
cell_colors_LS[1:length(cell_colors_LS)] <- c("#D62728","#FF9896","#9467BD","#C5B0D5","#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F")

#ED665D”,  “#AD8BC9”
cell_colors_LS.Sept <- cell_colors.ls[grep("LS|Sept",names(cell_colors.ls))]
cell_colors_LS.Sept[1:length(cell_colors_LS.Sept)] <- c("#D62728","#FF9896","#9467BD","#C5B0D5","#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#ED665D", "#AD8BC9","#7F7F7F")

cell_colors_neur <- cell_colors.ls[grep("LS|Sept|MS|TT|TNoS",names(cell_colors.ls))]
cell_colors_neur[1:length(cell_colors_neur)] <- c("#D62728","#FF9896","#9467BD","#C5B0D5","#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#BCBD22", "#DBDB8D", "#ED665D", "#AD8BC9","#7F7F7F", "#6DCCDA", "#999999", "#E69F00", "#56B4E9")

###############################################################################



############################ Plot LS and broad PCs ############################

colData(sce_pseudo_all)$cell.group <- "non neuronal"
colData(sce_pseudo_all)$cell.group[colData(sce_pseudo_all)$cellType.broad %in% neuronal_names] <- "neuronal"

## Broad clusters
p_all <- plotPCA(
    sce_pseudo_all,
    colour_by = "cellType.broad",
    shape_by = "cell.group",
    ncomponents = 4,
    point_size = 2,
    label_format = c("%s %02i", " (%i%%)"),
    percentVar = (metadata(sce_pseudo_all)$PCA_var_explained) * 100
    ) + scale_color_manual("cellType.broad", values = cell_colors_all)

pdf(file = here("snRNAseq_mouse/plots/08_pseudo_bulking/sce_pseudo_broad_PCs.pdf"), width = 8, height = 8)
print(p_all)
dev.off()

## Broad clusters colored by cell.group
p_all <- plotPCA(
    sce_pseudo_all,
    colour_by = "cell.group",
    ncomponents = 4,
    point_size = 2,
    label_format = c("%s %02i", " (%i%%)"),
    percentVar = (metadata(sce_pseudo_all)$PCA_var_explained) * 100
    ) + scale_color_manual("cellType.broad", values = cell_colors_all)

pdf(file = here("snRNAseq_mouse/plots/08_pseudo_bulking/sce_pseudo_broad_PCs_v2.pdf"), width = 8, height = 8)
print(p_all)
dev.off()

## LS clusters
p_LS <- plotPCA(
    sce_pseudo_LS,
    colour_by = "cellType.final",
    ncomponents = 4,
    point_size = 2,
    label_format = c("%s %02i", " (%i%%)"),
    percentVar = (metadata(sce_pseudo_LS)$PCA_var_explained) * 100
    ) + scale_color_manual("cellType.final",values = cell_colors_LS)

pdf(file = here("snRNAseq_mouse/plots/08_pseudo_bulking/sce_pseudo_LS_PCs.pdf"), width = 8, height = 8)
print(p_LS)
dev.off()

## LS and Sept clusters
p_LS_Sept <- plotPCA(
    sce_pseudo_LS.Sept,
    colour_by = "cellType.final",
    ncomponents = 5,
    point_size = 2,
    label_format = c("%s %02i", " (%i%%)"),
    percentVar = (metadata(sce_pseudo_LS.Sept)$PCA_var_explained) * 100
    ) + scale_color_manual("cellType.final",values = cell_colors_LS.Sept)

pdf(file = here("snRNAseq_mouse/plots/08_pseudo_bulking/sce_pseudo_LS_Sept_PCs.pdf"), width = 10, height = 10)
print(p_LS_Sept)
dev.off()

p_neur <- plotPCA(
    sce_pseudo_neuronal,
    colour_by = "cellType.final",
    ncomponents = 5,
    point_size = 2,
    label_format = c("%s %02i", " (%i%%)"),
    percentVar = (metadata(sce_pseudo_neuronal)$PCA_var_explained) * 100
    ) + scale_color_manual("cellType.final", values = cell_colors_neur)

pdf(file = here("snRNAseq_mouse/plots/08_pseudo_bulking/sce_pseudo_neuronal_PCs.pdf"), width = 10, height = 10)
print(p_neur)
dev.off()

p_neur <- plotPCA(
    sce_pseudo_neuronal.broad,
    colour_by = "cellType.broad",
    ncomponents = 12,
    point_size = 1,
    label_format = c("%s %02i", " (%i%%)"),
    percentVar = (metadata(sce_pseudo_neuronal.broad)$PCA_var_explained) * 100
    )
#+ scale_color_manual("cellType.broad", values = cell_colors_all)

pdf(file = here("snRNAseq_mouse/plots/08_pseudo_bulking/sce_pseudo_neuronal_broad_PCs.pdf"), width = 14, height = 14)
print(p_neur)
dev.off()


###############################################################################

