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

