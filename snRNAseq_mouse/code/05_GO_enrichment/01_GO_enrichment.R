## Call libraries
library("here")
library("sessioninfo")
library("dplyr")
library("data.table")
library("clusterProfiler")
library("org.Mm.eg.db")
library("ggplot2")

######################
#### Load DE data ####
######################

sigGene <- fread(file = "/dcl01/lieber/ajaffe/Keri/TrkBKO/tables/sigGenes_DE_FDR05.csv", header = TRUE, sep = ",", data.table = FALSE, stringsAsFactors = FALSE)

geneUniverse <- fread(file = "/dcl01/lieber/ajaffe/Keri/TrkBKO/tables/allGenes_DE_noCut.csv", header = TRUE, sep = ",", data.table = FALSE, stringsAsFactors = FALSE)

sigGene <- split(sigGene$EntrezID, sign(sigGene$logFC))
sigGene <- lapply(sigGene, function(x) x[!is.na(x)])

geneUniverse <- as.character(geneUniverse$EntrezID)
geneUniverse <- geneUniverse[!is.na(geneUniverse)]


#############################################
#### Run GO and KEGG enrichment analysis ####
#############################################

go <- compareCluster(sigGene,
    fun = "enrichGO",
    universe = geneUniverse,
    OrgDb = org.Mm.eg.db,
    ont = "ALL",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.05,
    readable = TRUE
)

kegg <- compareCluster(sigGene,
    fun = "enrichKEGG",
    universe = geneUniverse,
    organism = "mmu",
    pvalueCutoff = 0.1,
    qvalueCutoff = 0.05
)


#################################################
#### Save rda with gene ontology enrichemnts ####
#################################################

save(go, kegg, file = here(
    "snRNAseq_mouse",
    "processed_data",
    "SCE",
    "gene_ontology_enrichment_objects.rda"
))

# bp_out<-rbind(as.data.frame(filter(go, ONTOLOGY == "BP", Cluster == -1)) %>% dplyr::select(!"geneID") %>% head(50), as.data.frame(filter(go, ONTOLOGY == "BP", Cluster == 1)) %>% dplyr::select(!"geneID") %>% head(50))
# cc_out<-rbind(as.data.frame(filter(go, ONTOLOGY == "CC", Cluster == -1)) %>% dplyr::select(!"geneID") %>% head(50), as.data.frame(filter(go, ONTOLOGY == "CC", Cluster == 1)) %>% dplyr::select(!"geneID") %>% head(50))
# mf_out<-rbind(as.data.frame(filter(go, ONTOLOGY == "MF", Cluster == -1)) %>% dplyr::select(!"geneID") %>% head(50), as.data.frame(filter(go, ONTOLOGY == "MF", Cluster == 1)) %>% dplyr::select(!"geneID") %>% head(50))
# kegg_out<-rbind(as.data.frame(filter(kegg, Cluster == -1)) %>% dplyr::select(!"geneID") %>% head(n = 50), as.data.frame(filter(kegg, Cluster == 1)) %>% dplyr::select(!"geneID") %>% head(n = 50))
# write.csv(bp_out, here("snRNAseq_mouse/processed_data/tables/GOterms-BP_top50.csv"))
# write.csv(cc_out, here("snRNAseq_mouse/processed_data/tables/GOterms-CC_top50.csv"))
# write.csv(mf_out, here("snRNAseq_mouse/processed_data/tables/GOterms-MF_top50.csv"))
# write.csv(kegg_out, here("snRNAseq_mouse/processed_data/tables/KEGGterms_top50.csv"))


############################
#### Plot BP, CC and MF ####
############################

plot_go <- function(ont, title_p, path, filename) {
    dotplot_1 <- ggplot(filter(go, ONTOLOGY == ont), aes(Cluster, Description)) +
        theme_bw() +
        geom_point(aes(color = p.adjust, size = Count)) +
        scale_color_gradientn(
            colours = c("#f7ca64", "#46bac2", "#7e62a3"),
            trans = "log10",
            guide = guide_colorbar(reverse = TRUE, order = 1)
        ) +
        scale_size_continuous(range = c(2, 10)) +
        xlab("Cluster") +
        ylab("") +
        ggtitle(title_p)

    ggsave(filename = filename, path = path, dotplot_1, height = 6, width = 5)
}

plot_go(ont = "BP", title_p = "Biological Process", filename = "GOenrichment_BP.pdf", path = here("snRNAseq_mouse", "plots/", "05_GO_enrichment/"))

plot_go(ont = "CC", title_p = "Cellular Component", filename = "GOenrichment_CC.pdf", path = here("snRNAseq_mouse", "plots/", "05_GO_enrichment/"))

plot_go(ont = "MF", title_p = "Molecular Function", filename = "GOenrichment_MF.pdf", path = here("snRNAseq_mouse", "plots/", "05_GO_enrichment/"))


###################
#### Plot KEGG ####
###################

kegg@compareClusterResult$Description <- gsub(kegg@compareClusterResult$Description, pattern = " \\- Mus musculus \\(house mouse\\)", replacement = "")

dotplot_1 <- ggplot(kegg, aes(Cluster, Description)) +
    theme_bw() +
    geom_point(aes(color = p.adjust, size = Count)) +
    scale_color_gradientn(
        colours = c("#f7ca64", "#46bac2", "#7e62a3"),
        trans = "log10",
        guide = guide_colorbar(reverse = TRUE, order = 1)
    ) +
    scale_size_continuous(range = c(2, 10)) +
    xlab("Cluster") +
    ylab("") +
    ggtitle("KEGG")

ggsave(filename = "KEGGenrichment.pdf", path = here("snRNAseq_mouse", "plots/", "05_GO_enrichment/"), dotplot_1, height = 6, width = 5)


#####################################
#### Reproducibility information ####
#####################################

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
