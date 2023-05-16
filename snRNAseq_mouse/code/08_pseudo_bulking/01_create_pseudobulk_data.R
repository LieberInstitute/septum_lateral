library("here")
library("spatialLIBD")
library("SingleCellExperiment")
library("scater")
library("jaffelab")
library("CATALYST")
library("sessioninfo")



############################### Load sce object ###############################

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

sce.ls

###############################################################################



###################### Function to perform pseudobulking ######################

do_pseudobulk <- function(sce, cell_cluster) {
    sce_pseudo <-
        registration_pseudobulk(sce,
            var_registration = cell_cluster,
            var_sample_id = "Sample",
            min_ncells = 10
        )

    ## Compute PCs
    pca <- prcomp(t(assays(sce_pseudo)$logcounts))
    metadata(sce_pseudo) <- list("PCA_var_explained" = (summary(pca))$importance[2, 1:20])
    metadata(sce_pseudo)
    pca_pseudo <- pca$x[, seq_len(20)]
    colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
    reducedDims(sce_pseudo) <- list(PCA = pca_pseudo)

    ## Compute some reduced dims
    set.seed(20230509)
    sce_pseudo <- scater::runMDS(sce_pseudo, ncomponents = 20)
    sce_pseudo <- scater::runPCA(sce_pseudo, name = "runPCA")

    return(sce_pseudo)
}

###############################################################################



############################## Filter sce object ##############################

## sce.ls.filter has all cell types except the ones that needed to be dropped
## sce.ls.LS has only LS clusters

all_names <- levels(colData(sce.ls)$cellType.broad)[-grep("drop|mixed", levels(colData(sce.ls)$cellType.final))]
all_names <- all_names[-grep("mixed", all_names)]
LS_names <- levels(colData(sce.ls)$cellType.final)[grep("LS", levels(colData(sce.ls)$cellType.final))]
neuronal_names <- c("Chol", "LS", "Sept", "Str", "MS", "TNoS", "TT.IG.SH", "Thal", "IoC")

sce.ls.filter <- filterSCE(sce.ls, cellType.broad %in% all_names)
sce.ls.LS <- filterSCE(sce.ls, cellType.final %in% LS_names)
sce.ls.neuronal <- filterSCE(sce.ls, cellType.broad %in% neuronal_names)

###############################################################################



############################# pseudobulking for LS ############################

## pseudobulking across LS clusters
sce_pseudo_LS <- do_pseudobulk(sce.ls.LS, "cellType.final")
# 2023-05-16 11:47:54 make pseudobulk object
# 2023-05-16 11:47:56 dropping 6 pseudo-bulked samples that are below 'min_ncells'.
# 2023-05-16 11:47:56 drop lowly expressed genes
# 2023-05-16 11:47:57 normalize expression
# Warning in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth = TRUE,  :
#   You're computing too large a percentage of total singular values, use a standard svd instead.
# Warning message:
# In check_numbers(k = k, nu = nu, nv = nv, limit = min(dim(x)) -  :
#   more singular values/vectors requested than available


## pseudobulking for broad cell type clusters
sce_pseudo_all <- do_pseudobulk(sce.ls.filter, "cellType.broad")
# 2023-05-16 11:49:46 make pseudobulk object
# 2023-05-16 11:50:03 dropping 5 pseudo-bulked samples that are below 'min_ncells'.
# 2023-05-16 11:50:03 drop lowly expressed genes
# 2023-05-16 11:50:04 normalize expression
# Warning in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth = TRUE,  :
#   You're computing too large a percentage of total singular values, use a standard svd instead.


## pseudobulking for neuronal broad clusters
sce_pseudo_neuronal <- do_pseudobulk(sce.ls.neuronal, "cellType.broad")
# 2023-05-16 11:53:41 make pseudobulk object
# 2023-05-16 11:53:50 dropping 1 pseudo-bulked samples that are below 'min_ncells'.
# 2023-05-16 11:53:50 drop lowly expressed genes
# 2023-05-16 11:53:50 normalize expression
# Warning in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth = TRUE,  :
#   You're computing too large a percentage of total singular values, use a standard svd instead.
# Warning message:
# In check_numbers(k = k, nu = nu, nv = nv, limit = min(dim(x)) -  :
#   more singular values/vectors requested than available

###############################################################################



#################### Save both pseudobulking results to rda ###################

save(sce_pseudo_LS, sce_pseudo_all, sce_pseudo_neuronal, file = here(
    "snRNAseq_mouse",
    "processed_data",
    "SCE",
    "sce_pseudobulking_LS_and_broad.rda"
))

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
