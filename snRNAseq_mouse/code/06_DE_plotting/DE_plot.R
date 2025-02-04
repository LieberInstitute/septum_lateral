library("here")
library("sessioninfo")
library("dplyr")
library("data.table")
library("ggplot2")
library("EnhancedVolcano")
library("tibble")
library("tidyr")
library("scater")



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



############# Create groups of genes according to their function ##############

## This classification was made by Lionel Rodrqiguez based on what we want to
## talk about in the paper. The orange genes are LS specific, the green ones are
## related to neurodevelopment, and the blue group has plasticitiy/synaptic genes

glia_genes <- c("ENSMUSG00000007613", "ENSMUSG00000018008", "ENSMUSG00000018654", "ENSMUSG00000020143", "ENSMUSG00000026288")
orange_genes <- c("ENSMUSG00000022285", "ENSMUSG00000028524", "ENSMUSG00000029405", "ENSMUSG00000034796", "ENSMUSG00000034958", "ENSMUSG00000037492", "ENSMUSG00000046178")
green_genes <- c("ENSMUSG00000008658", "ENSMUSG00000020297", "ENSMUSG00000033676", "ENSMUSG00000049583", "ENSMUSG00000056755")
blue_genes <- c("ENSMUSG00000021373", "ENSMUSG00000021448", "ENSMUSG00000028176", "ENSMUSG00000037386", "ENSMUSG00000039419", "ENSMUSG00000041852", "ENSMUSG00000062296")

###############################################################################



#################### Genes to highlight and color selection ###################

## Selection of genes to highlight
downgene_df <- c(orange_genes, green_genes, blue_genes)
upgene_df <- glia_genes
#upgene_df <- as.vector(upgene_df$ensemblID)
selected <- c(downgene_df, upgene_df)
genes2plot <- sigGenes %>% filter(ensemblID %in% selected)

## Colors for the significant and not significant genes
keyvals <- ifelse(
    outGenes_plot$P.Value > 1e-02, "#f0e3d6", "#E2C6A7"
)

## Assigning colors for each groups of highlited genes
keyvals[outGenes_plot$ensemblID %in% orange_genes & outGenes_plot$ensemblID %in% selected] <- "#FB8500"
keyvals[outGenes_plot$ensemblID %in% green_genes & outGenes_plot$ensemblID %in% selected] <- "#A5C0DF"
keyvals[outGenes_plot$ensemblID %in% blue_genes & outGenes_plot$ensemblID %in% selected] <- "#789C25"
keyvals[outGenes_plot$ensemblID %in% upgene_df] <- "#006164"


## Legend names
names(keyvals)[keyvals == "#E2C6A7"] <- "FDR < 0.05"
names(keyvals)[keyvals == "#f0e3d6"] <- "Not significant"
names(keyvals)[keyvals == "#FB8500"] <- "Lateral Septum specific genes"
names(keyvals)[keyvals == "#789C25"] <- "Neurodevelopmental genes"
names(keyvals)[keyvals == "#A5C0DF"] <- "Plasticity genes"
names(keyvals)[keyvals == "#006164"] <- "Microglia specific genes"

###############################################################################



################################ Volcano plot #################################

volcano_plot <- EnhancedVolcano(outGenes_plot,
    x = "logFC",
    y = "P.Value",
    pCutoff = 1e-02,
    FCcutoff = 0,
    lab = rownames(outGenes_plot),
    selectLab = genes2plot$Symbol,
    labSize = 6.0,
    drawConnectors = TRUE,
    arrowheads = FALSE,
    max.overlaps = Inf,
    labCol = "black",
    pointSize = c(ifelse(outGenes_plot$ensemblID %in% selected, 6, 2)),
    colAlpha = c(ifelse(outGenes_plot$ensemblID %in% selected, 1, 2 / 5)),
    colCustom = keyvals,
    caption = paste0("total = ", nrow(outGenes_plot), " genes"),
    title = "",
    subtitle = "",
    legendPosition = "bottom"
) +
    ylim(c(0, 8)) +
    coord_flip()

pdf(here(
    "snRNAseq_mouse",
    "plots",
    "06_DE_plotting/",
    "volcano_plot_flip_TrkB-KD.pdf"
), height = 10, width = 14)
print(volcano_plot)
dev.off()


volcano_plot <- EnhancedVolcano(outGenes_plot,
    x = "logFC",
    y = "P.Value",
    pCutoff = 1e-02,
    FCcutoff = 0,
    lab = rownames(outGenes_plot),
    selectLab = c(""),
    labSize = 6.0,
    drawConnectors = TRUE,
    arrowheads = FALSE,
    max.overlaps = Inf,
    labCol = "black",
    pointSize = c(ifelse(outGenes_plot$ensemblID %in% selected, 6, 2)),
    colAlpha = c(ifelse(outGenes_plot$ensemblID %in% selected, 1, 2 / 5)),
    colCustom = keyvals,
    caption = paste0("total = ", nrow(outGenes_plot), " genes"),
    title = "",
    subtitle = "",
    legendPosition = "bottom"
) +
    ylim(c(0, 8)) +
    coord_flip()

pdf(here(
    "snRNAseq_mouse",
    "plots",
    "06_DE_plotting/",
    "volcano_plot_flip-nonames_TrkB-KD.pdf"
), height = 10, width = 14)
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
