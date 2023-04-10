library("here")
library("sessioninfo")
library("dplyr")
library("data.table")
library("ggplot2")
library("EnhancedVolcano")
library("tibble")
library("tidyr")


################################ Load DE data #################################

sigGenes <- fread(
    file = "/dcl01/lieber/ajaffe/Keri/TrkBKO/tables/sigGenes_DE_FDR05.csv",
    header = TRUE,
    sep = ",",
    data.table = FALSE,
    stringsAsFactors = FALSE
)

outGenes <- fread(
    file = "/dcl01/lieber/ajaffe/Keri/TrkBKO/tables/allGenes_DE_noCut.csv",
    header = TRUE,
    sep = ",",
    data.table = FALSE,
    stringsAsFactors = FALSE
)

###############################################################################



################################ Volcano plot #################################

outGenes_plot <- outGenes %>%
    select(logFC, P.Value, adj.P.Val, ensemblID, Symbol)
rownames(outGenes_plot) <- outGenes_plot$ensemblID

volcano_plot <- EnhancedVolcano(outGenes_plot,
    lab = rownames(outGenes_plot),
    x = "logFC",
    y = "P.Value",
    pCutoff = 1e-02,
    FCcutoff = 0,
    legendPosition = "right",
    legendLabels = c(
        "Not sig.", "Not sig.", "FDR < 0.05",
        "FDR < 0.05"
    ),
    col = c("grey30", "grey30", "royalblue", "red2")
)

pdf("~/volcano.pdf", height = 10, width = 12)
print(volcano_plot)
dev.off()

###############################################################################



######################## Create csv with top-n DE genes #######################

topgenes <- function(genes, p_value_thresh, gene_set, top, csv_name) {
    if (gene_set == "+") {
        gene_df <- genes %>%
            filter(adj.P.Val < p_value_thresh, logFC > 0) %>%
            arrange(desc(abs(logFC))) %>%
            head(top) %>%
            dplyr::select(ensemblID, Symbol, EntrezID, logFC)
    } else {
        gene_df <- genes %>%
            filter(adj.P.Val < p_value_thresh, logFC < 0) %>%
            arrange(desc(abs(logFC))) %>%
            head(top) %>%
            dplyr::select(ensemblID, Symbol, EntrezID, logFC)
    }

    write.csv(
        x = gene_df,
        file = csv_name,
        quote = FALSE,
        row.names = FALSE
    )
}

topgenes(genes = sigGenes, p_value_thresh = 0.01, gene_set = "+", top = 100, csv_name = here("snRNAseq_mouse", "processed_data", "tables", "DEgenes_FDR01_pos_top100.csv"))

topgenes(genes = sigGenes, p_value_thresh = 0.01, gene_set = "-", top = 100, csv_name = here("snRNAseq_mouse", "processed_data", "tables", "DEgenes_FDR01_neg_top100.csv"))

topgenes(genes = sigGenes, p_value_thresh = 0.05, gene_set = "+", top = 100, csv_name = here("snRNAseq_mouse", "processed_data", "tables", "DEgenes_FDR05_pos_top100.csv"))

topgenes(genes = sigGenes, p_value_thresh = 0.05, gene_set = "-", top = 100, csv_name = here("snRNAseq_mouse", "processed_data", "tables", "DEgenes_FDR05_neg_top100.csv"))

###############################################################################



######################### Reproducibility information #########################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

###############################################################################
