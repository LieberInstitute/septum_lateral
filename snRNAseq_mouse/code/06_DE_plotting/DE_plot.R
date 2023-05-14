library("here")
library("sessioninfo")
library("dplyr")
library("data.table")
library("ggplot2")
library("EnhancedVolcano")
library("tibble")
library("tidyr")
library("scater")

topgenes <- function(genes, p_value_thresh, top) {
    upgene_df <- genes %>%
        filter(adj.P.Val < p_value_thresh, logFC > 0) %>%
        arrange(desc(abs(logFC))) %>%
        head(top) %>%
        dplyr::select(ensemblID, Symbol, EntrezID, logFC)
    downgene_df <- genes %>%
        filter(adj.P.Val < p_value_thresh, logFC < 0) %>%
        arrange(desc(abs(logFC))) %>%
        head(top) %>%
        dplyr::select(ensemblID, Symbol, EntrezID, logFC)
    gene_df <- rbind(upgene_df, downgene_df)

    return(gene_df)
}

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



################################# Format data #################################

outGenes_plot <- outGenes %>%
    select(logFC, P.Value, adj.P.Val, ensemblID, Symbol)

rownames(outGenes_plot) <- uniquifyFeatureNames(
    outGenes_plot$ensemblID,
    outGenes_plot$Symbol
)

###############################################################################



################################ Volcano plot #################################

downgene_df <- c("ENSMUSG00000022285", "ENSMUSG00000028524", "ENSMUSG00000029405", "ENSMUSG00000034796", "ENSMUSG00000034958", "ENSMUSG00000037492", "ENSMUSG00000046178", "ENSMUSG00000028176", "ENSMUSG00000033676")
upgene_df <- sigGenes %>%
    filter(adj.P.Val < 0.01, logFC > 0) %>%
    arrange(desc(abs(logFC))) %>%
    head(5) %>%
    select(ensemblID)
upgene_df <- as.vector(upgene_df$ensemblID)
selected <- c(downgene_df, upgene_df)

genes2plot <- sigGenes %>% filter(ensemblID %in% selected)
#genes2plot <- topgenes(genes = sigGenes, p_value_thresh = 0.01, top = 5)
n_overlaps <- ifelse(length(selected) > 0, Inf, 15)

volcano_plot <- EnhancedVolcano(outGenes_plot,
    x = "logFC",
    y = "P.Value",
    pCutoff = 1e-02,
    FCcutoff = 0,
    selectLab = genes2plot$Symbol,
    max.overlaps = n_overlaps,
    drawConnectors = TRUE,
    lab = rownames(outGenes_plot),
    labSize = 6.0,
    labCol = "black",
    labFace = "bold",
    boxedLabels = TRUE,
    colAlpha = 4 / 5,
    col = c("#a8b6cc", "#a8b6cc", "#1f449c", "#1f449c"),
     legendLabels = c(
            "Not sig.", "Not sig.", "FDR < 0.05",
            "FDR < 0.05"
        ),
     caption = paste0("total = ", nrow(outGenes_plot), " genes"),
    title = "",
    subtitle = "",
    legendPosition = "right"
) +
    ylim(c(0, 8)) +
    coord_flip()

pdf(here("snRNAseq_mouse",
        "plots",
        "06_DE_plotting/",
        "volcano_plot_flip_TrkB-KD.pdf"), height = 10, width = 14)
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
