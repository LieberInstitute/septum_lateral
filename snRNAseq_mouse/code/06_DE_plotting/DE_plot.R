library("here")
library("sessioninfo")
library("dplyr")
library("data.table")
library("ggplot2")
library("EnhancedVolcano")


######################
#### Load DE data ####
######################

sigGene <- fread(file = "/dcl01/lieber/ajaffe/Keri/TrkBKO/tables/sigGenes_DE_FDR05.csv", header = TRUE, sep = ",", data.table = FALSE, stringsAsFactors = FALSE)


########################################
#### Create csv with top-n DE genes ####
########################################

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
