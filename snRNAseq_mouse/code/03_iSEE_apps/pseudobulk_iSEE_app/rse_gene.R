library("SingleCellExperiment")
library("here")
library(SummarizedExperiment)
library(edgeR)


here()
setwd(here("snRNAseq_mouse","code", "03_iSEE_apps", "pseudobulk_iSEE_app"))
load("rse_gene_TrkB_KO_LS_n8.Rdata", verbose = TRUE)
SummarizedExperiment::colData(rse_gene)

meta <- read.csv(here("snRNAseq_mouse","raw_data", "sample_info", "nsh-sampleList.txt"))


rownames(meta) = paste0("Sample.", meta$SampleID)
repmeta <- as(cbind(meta[colnames(rse_gene),], colData(rse_gene)), "DFrame")
colData(rse_gene) <- repmeta

## Genes, Computing logcounts
assays(rse_gene, withDimnames=FALSE)$logcounts <- edgeR::cpm(calcNormFactors(rse_gene, method = "TMM"), log = TRUE, prior.count = 0.5)

save(rse_gene, getRPKM, file = here("snRNAseq_mouse","code", "03_iSEE_apps", "pseudobulk_iSEE_app", "rse_gene_TrkB_KO_LS_n8_wm.Rdata"))
