#############################
#### Reshape 1vsAll data ####
#############################

reshape_1vsAll <- function(OnevsAll) {
    ## Select genes with non0median == TRUE
    non0med_genes <- lapply(OnevsAll, function(x) {
        rownames(x[[2]][x[[2]]$non0median == TRUE, ])
    })

    non0med_genes <- unique(unlist(non0med_genes))
    non0med_genes <- non0med_genes[order(non0med_genes)]

    ## Change pvalues. fdrs and t-stats for genes non0median == FALSE
    OnevsAll_modified <- lapply(OnevsAll, function(celltype) {
        enriched <- celltype[[2]][non0med_genes, ]
        enriched$std.logFC[!enriched$non0median] <- 0
        enriched$log.p.value[!enriched$non0median] <- log(1)
        enriched$log.FDR[!enriched$non0median] <- log(1)
        return(enriched)
    })

    ## Unlog p-values and FDRs
    OnevsAll_modified <- lapply(OnevsAll_modified, function(enriched) {
        res <- enriched
        res$FDR <- exp(enriched$log.FDR)
        res$p.value <- exp(enriched$log.p.value)
        res <- res[, c("std.logFC", "p.value", "FDR")]
        return(res)
    })

    ## Change column names
    OnevsAll_modified <- mapply(function(res, names_ct) {
        colnames(res) <- paste0(c("t_stat_", "p_value_", "fdr_"), names_ct)
        res$ensembl <- rownames(res)
        return(res)
    }, OnevsAll_modified, names(OnevsAll_modified))

    ## Convert to data.frame
    modeling_result_enrichment <-
        as.data.frame(Reduce(
            function(...) {
                merge(..., by = "ensembl")
            },
            OnevsAll_modified
        ))

    ## Add names from mgi data base
    mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
    genes <- modeling_result_enrichment$ensembl
    gene_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "mgi_symbol"), values = genes, mart = mart)
    modeling_result_enrichment <- merge(modeling_result_enrichment, gene_list, by.x = "ensembl", by.y = "ensembl_gene_id")
    l_names <- length(colnames(modeling_result_enrichment))
    modeling_result_enrichment <- modeling_result_enrichment[, c(2:(l_names - 1), 1, l_names)]
    modeling_result_enrichment <- dplyr::rename(modeling_result_enrichment, gene = mgi_symbol)

    return(modeling_result_enrichment)
}
