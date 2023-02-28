library("here")
library("sessioninfo")
library("biomaRt")

## Load markers by Tran et al
load(
    here(
        "snRNAseq_mouse",
        "processed_data",
        "SCE",
        "markers-stats_LS-n4_findMarkers_33cellTypes.broad.rda"
    ),
    verbose = TRUE
)
# Loading objects:
#     markers.ls.t.pw.broad
#     markers.ls.t.1vAll.broad
#     medianNon0.ls.broad

lobstr::obj_size(markers.ls.t.pw.broad)
# 423.75 MB
lobstr::obj_size(markers.ls.t.1vAll.broad)
# 63.08 MB

class(markers.ls.t.1vAll.broad)
names(markers.ls.t.1vAll.broad)
class(markers.ls.t.1vAll.broad$LS)
names(markers.ls.t.1vAll.broad$LS)
class(markers.ls.t.1vAll.broad$LS$LS_enriched)
markers.ls.t.1vAll.broad$LS$LS_enriched

colnames(markers.ls.t.1vAll.broad$LS$LS_enriched)[1:3]<-c("t_stat", "p_value", "fdr")
markers.ls.t.1vAll.broad$LS$LS_enriched$ensembl<-rownames(markers.ls.t.1vAll.broad$LS$LS_enriched)

mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
genes <- markers.ls.t.1vAll.broad$LS$LS_enriched$ensembl
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","mgi_symbol"),values=genes,mart= mart)
LS_enriched<-merge(markers.ls.t.1vAll.broad$LS$LS_enriched,G_list,by.x="ensembl",by.y="ensembl_gene_id")
LS_enriched<-LS_enriched[,c(2,3,4,1,6)]
colnames(LS_enriched)[5]<-"gene"


modeling_results <- list(
    "enrichment" = LS_enriched
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
