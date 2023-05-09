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



############################## Filter sce object ##############################

## sce.ls.filter has all cell types except the ones that needed to be dropped
## sce.ls.LS has only LS clusters

all_names <- levels(colData(sce.ls)$cellType.final)[-grep("drop|mixed", levels(colData(sce.ls)$cellType.final))]
LS_names <- levels(colData(sce.ls)$cellType.final)[grep("LS", levels(colData(sce.ls)$cellType.final))]

sce.ls.filter <- filterSCE(sce.ls, cellType.final %in% all_names)
sce.ls.LS <- filterSCE(sce.ls, cellType.final %in% LS_names)

###############################################################################



############################# pseudobulking for LS ############################

## pseudobulking across LS clusters
sce_pseudo_LS <-
    registration_pseudobulk(sce.ls.LS,
        var_registration = "cellType.final",
        var_sample_id = "Sample",
        min_ncells = 10
    )

dim(sce_pseudo_LS)
colData(sce_pseudo_LS)

## Compute PCs
pca <- prcomp(t(assays(sce_pseudo_LS)$logcounts))
#metadata(sce_pseudo_LS) <- list("PCA_var_explained" = jaffelab::getPcaVars(pca)[seq_len(20)])
#metadata(sce_pseudo_LS)
pca_pseudo <- pca$x[, seq_len(20)]
colnames(pca_pseudo) <- paste0("PC", sprintf("%02d", seq_len(ncol(pca_pseudo))))
reducedDims(sce_pseudo_LS) <- list(PCA = pca_pseudo)

## Compute some reduced dims
set.seed(20230509)
sce_pseudo_LS <- scater::runMDS(sce_pseudo_LS, ncomponents = 20)
sce_pseudo_LS <- scater::runPCA(sce_pseudo_LS, name = "runPCA")
# Warning in (function (A, nv = 5, nu = nv, maxit = 1000, work = nv + 7, reorth = TRUE,  :
#   You're computing too large a percentage of total singular values, use a standard svd instead.
# Warning message:
# In check_numbers(k = k, nu = nu, nv = nv, limit = min(dim(x)) -  :
#   more singular values/vectors requested than available

###############################################################################
