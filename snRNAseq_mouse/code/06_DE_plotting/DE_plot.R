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
