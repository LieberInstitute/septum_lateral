## Call libraries
library("here")
library("sessioninfo")
library("dplyr")
library("data.table")
library("clusterProfiler")
library("org.Mm.eg.db")

######################
#### Load DE data ####
######################

sigGene <- fread(file = "/dcl01/lieber/ajaffe/Keri/TrkBKO/tables/sigGenes_DE_FDR05.csv", header = TRUE, sep = ",", data.table = FALSE, stringsAsFactors = FALSE)

geneUniverse <- fread(file = "/dcl01/lieber/ajaffe/Keri/TrkBKO/tables/allGenes_DE_noCut.csv", header = TRUE, sep = ",", data.table = FALSE, stringsAsFactors = FALSE)

go <- compareCluster(sigGene, fun = "enrichGO",
                universe = geneUniverse, OrgDb = org.Mm.eg.db,
                ont = "ALL", pAdjustMethod = "BH",
                pvalueCutoff  = 0.1, qvalueCutoff  = 0.05,
				readable= TRUE)
