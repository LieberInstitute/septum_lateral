library(SingleCellExperiment)
library(scater)
library(scran)
library(dplyr)
library(tidyr)
library(forcats)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)

load("/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/processed_data/SCE/sce_updated_LS.rda", verbose = T)
# sce.ls, annotationTab.ls, cell_colors.ls

# Check to see it's the right processed data
rowData(sce.ls)
colData(sce.ls)

# Drop the 2 populations we don't want shown
sce.ls <- sce.ls[, -grep("drop", sce.ls$cellType.final)]
sce.ls <- sce.ls[, -grep("mixed", sce.ls$cellType.final)]
sce.ls$cellType.final <- droplevels(sce.ls$cellType.final)

# to check if it dropped
levels(sce.ls$cellType.final)

# change Ens ID to gene names
rownames(sce.ls) <- uniquifyFeatureNames(rowData(sce.ls)$gene_id, rowData(sce.ls)$gene_name)

# Code adapted from: https://github.com/lmweber/locus-c/blob/main/code/analyses_snRNAseq/05b_cluster_identification.R

# --------------
# Color palettes
# --------------

# color palettes from scater package
tableau20 <- c(
    "#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
    "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
    "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
    "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5"
)

tableau10medium <- c(
    "#729ECE", "#FF9E4A", "#67BF5C", "#ED665D", "#AD8BC9",
    "#A8786E", "#ED97CA", "#A2A2A2", "#CDCC5D", "#6DCCDA"
)

colors <- unique(c(tableau20, tableau10medium))


# --------------
# Market Genes
# --------------

markers_all <- c(
    "Snap25", # Neuronal
    "Slc17a7", "Slc17a6", # Glutamatergic
    "Gad1", "Gad2", # GABAergcic
    "Mbp", "Mobp", # Oligo
    "Pdgfra", "Vcan", # Oligo precursors
    "Bcan", "Sox4", # OPC_COP
    "Slc1a2", "Aqp4", # Astrocyte
    "Flt1", "Cldn5", "Col1a2", # Endothelial/Mural
    "Cx3cr1", "Csf1r", # Microglia
    "Rbpms" # Ependymal
)

markers <- c(
    "Snap25", "Slc17a7", "Slc17a6", "Gad1", "Gad2", "Mbp", "Mobp", "Pdgfra",
    "Vcan", "Bcan", "Sox4", "Slc1a2", "Aqp4", "Flt1", "Cldn5",
    "Cx3cr1", "Csf1r", "Col1a2", "Rbpms"
)

# marker labels
marker_labels <- c(
    rep("neuron", 1),
    rep("excitatory", 2),
    rep("inhibitory", 2),
    rep("oligodendrocytes", 2),
    rep("OPC", 2),
    rep("OPC_COP", 2),
    rep("astrocytes", 2),
    rep("endothelial_mural", 3),
    rep("microglia", 2),
    rep("ependymal", 1)
)

marker_labels <-
    factor(marker_labels, levels = unique(marker_labels))

# colors: selected from tableau20 and tableau10medium
colors_markers <- list(marker = c(
    neuron = "black",
    excitatory = "#1F77B4",
    inhibitory = "#AEC7E8",
    oligodendrocytes = "#D62728",
    OPC = "#9467BD",
    OPC_COP = "#FF7F0E",
    astrocytes = "#9EDAE5",
    endothelial_mural = "#98DF8A",
    microglia = "#8C564B",
    ependymal = "#17BECF"
))

# cluster labels
cluster_pops <- list(
    excitatory = c(2, 29, 30, 31, 32, 33),
    inhibitory = c(6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 24, 25, 26, 27, 28),
    oligodendrocytes = 21,
    OPC = 22,
    OPC_COP = 23,
    astrocytes = 1,
    endothelial_mural = c(4, 19),
    microglia = 16,
    ependymal = 5,
    neuroblast = 20,
    choroid_plexus = 3
)

# cluster labels order
cluster_pops_order <- unname(unlist(cluster_pops))
# swap values and names of list
cluster_pops_rev <- rep(names(cluster_pops), times = sapply(cluster_pops, length))
names(cluster_pops_rev) <- unname(unlist(cluster_pops))
cluster_pops_rev <- cluster_pops_rev[as.character(sort(cluster_pops_order))]

cluster_pops_rev <- factor(cluster_pops_rev, levels = names(cluster_pops))

# second set of cluster labels
neuron_pops <- ifelse(
    cluster_pops_rev %in% c("excitatory", "inhibitory"),
    "neuronal", "non-neuronal"
) %>%
    factor(., levels = c("neuronal", "non-neuronal"))

# colors: selected from tableau20 and tableau10medium
colors_clusters <- list(population = c(
    excitatory = "#1F77B4",
    inhibitory = "#AEC7E8",
    oligodendrocytes = "gray60",
    OPC = "#D62728",
    OPC_COP = "#9467BD",
    astrocytes = "#FF7F0E",
    endothelial_mural = "#98DF8A",
    microglia = "#8C564B",
    ependymal = "#9EDAE5",
    neuroblast = "#17BECF",
    choroid_plexus = "#FF9E4A"
))

colors_neurons <- list(class = c(
    neuronal = "black",
    `non-neuronal` = "gray90"
))

n <- table(colData(sce.ls)$cellType.final)

# heatmap data

# using 'splitit' function from rafalib package
# code from Matthew N Tran
splitit <- function(x) split(seq(along = x), x)

cell_idx <- splitit(sce.ls$cellType.final)
dat <- as.matrix(logcounts(sce.ls))
rownames(dat) <- rowData(sce.ls)$gene_name


hm_mat <- t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers, i]))))

# row annotation
row_ha <- rowAnnotation(
    n = anno_barplot(as.numeric(n), gp = gpar(fill = "navy"), border = FALSE),
    class = neuron_pops,
    population = cluster_pops_rev,
    show_annotation_name = FALSE,
    col = c(colors_clusters, colors_neurons)
)

# column annotation
col_ha <- columnAnnotation(
    marker = marker_labels,
    show_annotation_name = FALSE,
    show_legend = FALSE,
    col = colors_markers
)


hm <- Heatmap(
    hm_mat,
    name = "mean\nlogcounts",
    column_title = "LS Clusters Mean Marker Expression",
    column_title_gp = gpar(fontface = "bold"),
    col = brewer.pal(n = 7, "OrRd"),
    right_annotation = row_ha,
    bottom_annotation = col_ha,
    row_order = cluster_pops_order,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = cluster_pops_rev,
    row_title = NULL,
    column_split = marker_labels,
    heatmap_legend_param = list(at = c(0, 2, 4, 6)),
    column_names_gp = gpar(fontface = "italic"),
    rect_gp = gpar(col = "gray50", lwd = 0.5)
)

# hm

pdf("/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/plots/Complex Heatmap All Clusters.pdf", width = 9)
hm
dev.off()


# Complex heatmap for just neuronal populations with region markers

# TT.IG.SH
# TNoS
# Thal
# Striatal
# MS
# LS
# IoC
# Chol

# Broad Marker Expression Heatmap for Fig 1
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

markers_all <- c(
    "Sv2b", "Pcsk5", "Satb2", # TT/IG/SH
    "Samd3", "Sema3a", # TNoS
    "Synpo2", # Thal
    "Bcl11b", "Rarb", "Tac1", "Penk", # Striatal
    "Trpc5", "Elavl2", # MS
    "Trpc4", "Homer2", "Ptpn3", "Trhde", "Cpne7", "Nrp1", # LS
    "Frmd3", "Ablim3", "Unc5d", # Sept
    "Pkib", "Drd3", # IoC
    "Chat" # Cholinergic
)

markers <- c(
    "Sv2b", "Pcsk5", "Satb2", "Samd3", "Sema3a", "Synpo2", "Bcl11b", "Rarb",
    "Tac1", "Penk", "Trpc5", "Elavl2", "Trpc4", "Homer2", "Ptpn3", "Trhde", "Cpne7",
    "Nrp1", "Frmd3", "Ablim3", "Unc5d", "Pkib", "Drd3", "Chat"
)



# marker labels
marker_labels <- c(
    rep("TT.IG.SH", 3),
    rep("TNoS", 2),
    rep("Thal", 1),
    rep("Striatal", 4),
    rep("MS", 2),
    rep("LS", 6),
    rep("Sept", 3),
    rep("IoC", 2),
    rep("Chol", 1)
)


marker_labels <-
    factor(marker_labels, levels = unique(marker_labels))

# colors: selected from tableau20 and tableau10medium
colors_markers <- list(marker = c(
    TT.IG.SH = "#8C564B",
    TNoS = "#1F77B4",
    Thal = "#AEC7E8",
    Striatal = "#D62728",
    MS = "#9467BD",
    LS = "#FF7F0E",
    Sept = "#17BECF",
    IoC = "#9EDAE5",
    Chol = "#98DF8A"
))

# cluster labels
cluster_pops <- list(
    TT.IG.SH = c(21, 22, 23),
    TNoS = 20,
    Thal = 19,
    Striatal = c(16, 17, 18),
    MS = c(12, 13),
    LS = c(3, 4, 5, 6, 7, 8, 9, 10, 11),
    Sept = c(14, 15),
    IoC = 2,
    Chol = 1
)


# cluster labels order
cluster_pops_order <- unname(unlist(cluster_pops))
# swap values and names of list
cluster_pops_rev <- rep(names(cluster_pops), times = sapply(cluster_pops, length))
names(cluster_pops_rev) <- unname(unlist(cluster_pops))
cluster_pops_rev <- cluster_pops_rev[as.character(sort(cluster_pops_order))]

cluster_pops_rev <- factor(cluster_pops_rev, levels = names(cluster_pops))

# second set of cluster labels
neuron_pops <- ifelse(
    cluster_pops_rev %in% c("LS", "MS", "Sept"),
    "Septal", "Non-Septal"
) %>%
    factor(., levels = c("Septal", "Non-Septal"))

# colors: selected from tableau20 and tableau10medium
colors_clusters <- list(population = c(
    TT.IG.SH = "#8C564B",
    TNoS = "#1F77B4",
    Thal = "#AEC7E8",
    Striatal = "#D62728",
    MS = "#9467BD",
    LS = "#FF7F0E",
    Sept = "#17BECF",
    IoC = "#9EDAE5",
    Chol = "#98DF8A"
))

colors_neurons <- list(class = c(
    Septal = "black",
    `Non-Septal` = "gray90"
))

n <- table(colData(neuron.sce.ls)$cellType.final)

# heatmap data

# using 'splitit' function from rafalib package
# code from Matthew N Tran
splitit <- function(x) split(seq(along = x), x)

cell_idx <- splitit(neuron.sce.ls$cellType.final)
dat <- as.matrix(logcounts(neuron.sce.ls))
rownames(dat) <- rowData(neuron.sce.ls)$gene_name


hm_mat <- t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers, i]))))

# row annotation
row_ha <- rowAnnotation(
    n = anno_barplot(as.numeric(n), gp = gpar(fill = "navy"), border = FALSE),
    class = neuron_pops,
    population = cluster_pops_rev,
    show_annotation_name = FALSE,
    col = c(colors_clusters, colors_neurons)
)

# column annotation
col_ha <- columnAnnotation(
    marker = marker_labels,
    show_annotation_name = FALSE,
    show_legend = FALSE,
    col = colors_markers
)


hm <- Heatmap(
    hm_mat,
    name = "mean\nlogcounts",
    column_title = "Neuronal Region Marker Mean Expression",
    column_title_gp = gpar(fontface = "bold"),
    col = brewer.pal(n = 7, "OrRd"),
    right_annotation = row_ha,
    bottom_annotation = col_ha,
    row_order = cluster_pops_order,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = cluster_pops_rev,
    row_title = NULL,
    column_split = marker_labels,
    # heatmap_legend_param = list(at = c( 0, 2, 4, 6)),
    column_names_gp = gpar(fontface = "italic"),
    rect_gp = gpar(col = "gray50", lwd = 0.5)
)

# hm

pdf("/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/plots/Complex Heatmap Neuronal Populations.pdf", width = 9)
hm
dev.off()


# Complex heatmap for just neuronal populations with region markers

# TT.IG.SH
# TNoS
# Thal
# Striatal
# MS
# LS
# IoC
# Chol

# Broad Marker Expression Heatmap for Fig 1

# Generate heatmap for only LS clusters
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




markers_all <- c(
    "Kcnmb2", "Gpr176", "Col19a1", # LS_C
    "Htr7", "Baiap3", # LS_D
    "Nos1", "Npy", "Sst", # LS_L
    "Stac2", "Grid2ip", "Trpc6", "Trpc3", # LS_M
    "Npas1", "Cnr1", "Kcnq4", "Reln", # LS_N
    "Pax6", "Gpc3", "Gabrg1", # LS_O
    "Rnf207", "St8sia6", # LS_P
    "Crhr2", "Tafa1", # LS_Q
    "Htr1b", "Slc18a2", "Col15a1", "Ano2", "Ntf3", # LS_R
    "Vipr2", "Cpa6", "Lgr5", # Sept_G
    "Baiap3", "Pdcd7", "Trpc5" # Sept_I
)

markers <- c(
    "Kcnmb2", "Gpr176", "Col19a1", "Htr7", "Baiap3", "Nos1", "Npy", "Sst", "Stac2", "Grid2ip", "Trpc6", "Trpc3",
    "Npas1", "Cnr1", "Kcnq4", "Reln", "Pax6", "Gpc3", "Gabrg1", "Rnf207", "St8sia6", "Crhr2", "Tafa1",
    "Htr1b", "Slc18a2", "Col15a1", "Ano2", "Ntf3", "Vipr2", "Cpa6", "Lgr5", "Baiap3", "Pdcd7", "Trpc5"
)


# marker labels
marker_labels <- c(
    rep("LS_C", 3),
    rep("LS_D", 2),
    rep("LS_L", 3),
    rep("LS_M", 4),
    rep("LS_N", 4),
    rep("LS_O", 3),
    rep("LS_P", 2),
    rep("LS_Q", 2),
    rep("LS_R", 5),
    rep("Sept_G", 3),
    rep("Sept_I", 3)
)


marker_labels <-
    factor(marker_labels, levels = unique(marker_labels))

# colors: selected from tableau20 and tableau10medium
colors_markers <- list(marker = c(
    LS_C = "#8C564B",
    LS_D = "#1F77B4",
    LS_L = "#AEC7E8",
    LS_M = "#D62728",
    LS_N = "#9467BD",
    LS_O = "#FF7F0E",
    LS_P = "#17BECF",
    LS_Q = "#9EDAE5",
    LS_R = "#98DF8A",
    Sept_G = "#8C564B",
    Sept_I = "#17BECF"
))

# cluster labels
cluster_pops <- list(
    LS_C = 1,
    LS_D = 2,
    LS_L = 3,
    LS_M = 4,
    LS_N = 5,
    LS_O = 6,
    LS_P = 7,
    LS_Q = 8,
    LS_R = 9,
    Sept_G = 10,
    Sept_I = 11
)


# cluster labels order
cluster_pops_order <- unname(unlist(cluster_pops))
# swap values and names of list
cluster_pops_rev <- rep(names(cluster_pops), times = sapply(cluster_pops, length))
names(cluster_pops_rev) <- unname(unlist(cluster_pops))
cluster_pops_rev <- cluster_pops_rev[as.character(sort(cluster_pops_order))]

cluster_pops_rev <- factor(cluster_pops_rev, levels = names(cluster_pops))

# second set of cluster labels
neuron_pops <- ifelse(
    cluster_pops_rev %in% c("Sept_G", "Sept_I"),
    "Broad Septal", "Lateral Septal"
) %>%
    factor(., levels = c("Broad Septal", "Lateral Septal"))

# colors: selected from tableau20 and tableau10medium
colors_clusters <- list(population = c(
    LS_C = "#8C564B",
    LS_D = "#1F77B4",
    LS_L = "#AEC7E8",
    LS_M = "#D62728",
    LS_N = "#9467BD",
    LS_O = "#FF7F0E",
    LS_P = "#17BECF",
    LS_Q = "#9EDAE5",
    LS_R = "#98DF8A",
    Sept_G = "#8C564B",
    Sept_I = "#17BECF"
))

colors_neurons <- list(class = c(
    `Broad Septal` = "black",
    `Lateral Septal` = "gray90"
))

n <- table(colData(ls.sce.ls)$cellType.final)

# heatmap data

# using 'splitit' function from rafalib package
# code from Matthew N Tran
splitit <- function(x) split(seq(along = x), x)

cell_idx <- splitit(ls.sce.ls$cellType.final)
dat <- as.matrix(logcounts(ls.sce.ls))
rownames(dat) <- rowData(ls.sce.ls)$gene_name


hm_mat <- t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers, i]))))

# row annotation
row_ha <- rowAnnotation(
    n = anno_barplot(as.numeric(n), gp = gpar(fill = "navy"), border = FALSE),
    class = neuron_pops,
    population = cluster_pops_rev,
    show_annotation_name = FALSE,
    col = c(colors_clusters, colors_neurons)
)

# column annotation
col_ha <- columnAnnotation(
    marker = marker_labels,
    show_annotation_name = FALSE,
    show_legend = FALSE,
    col = colors_markers
)


hm <- Heatmap(
    hm_mat,
    name = "mean\nlogcounts",
    column_title = "LS clusters mean marker expression",
    column_title_gp = gpar(fontface = "bold"),
    col = brewer.pal(n = 7, "OrRd"),
    right_annotation = row_ha,
    bottom_annotation = col_ha,
    row_order = cluster_pops_order,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_split = cluster_pops_rev,
    row_title = NULL,
    column_split = marker_labels,
    heatmap_legend_param = list(at = c(0, 2, 4, 6)),
    column_names_gp = gpar(fontface = "italic"),
    rect_gp = gpar(col = "gray50", lwd = 0.5)
)

# hm

pdf("/dcs04/lieber/marmaypag/pilotLS_LIBD1070/snRNAseq_mouse/plots/Complex Heatmap LS Cluster Markers.pdf", width = 9)
hm
dev.off()
