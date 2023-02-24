library("here")
library("sessioninfo")

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

modeling_results <- list(
    "enrichment" = todo
)

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
