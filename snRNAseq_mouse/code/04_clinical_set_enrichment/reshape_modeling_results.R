############################# Reshape 1vsAll data #############################

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
    OnevsAll_modified <-
        as.data.frame(Reduce(
            function(...) {
                merge(..., by = "ensembl")
            },
            OnevsAll_modified
        ))

    ## Add names from mgi data base
    modeling_result_enrichment <- add_gene_names(OnevsAll_modified)

    return(modeling_result_enrichment)
}



############################## Reshape 1vs1 data ##############################

reshape_1vs1 <- function(OnevsOne) {
    ## Change pvalues. fdrs and t-stats for genes non0median == FALSE
    OnevsOne_modified <- lapply(OnevsOne_modified, function(enriched) {
        enriched <- as.data.frame(enriched)
        enriched[enriched$non0median == FALSE, grep("FC", names(enriched))] <- 0
        enriched[enriched$non0median == FALSE, grep("value", names(enriched))] <- log(1)
        enriched[enriched$non0median == FALSE, grep("FDR", names(enriched))] <- log(1)
        enriched <- enriched %>% dplyr::select(-non0median)
        return(enriched)
    })

    ## Un log pvalues and FDR
    OnevsOne_modified <- lapply(OnevsOne_modified, function(enriched) {
        enriched <- enriched %>% mutate_at(vars(contains("FDR")), exp)
        enriched <- enriched %>% mutate_at(vars(contains("value")), exp)
        return(enriched)
    })

    ## Change column names
    OnevsOne_modified <- lapply(OnevsOne_modified, function(enriched) {
        names(enriched) <- gsub(names(enriched), pattern = "stats\\.", replacement = "__")
        return(enriched)
    })

    ## Convert to data frame
    OnevsOne_modified <- as.data.frame(OnevsOne_modified)

    ## Change column names
    names(OnevsOne_modified) <- gsub(names(OnevsOne_modified), pattern = "\\.\\_\\_", replacement = "-")
    names(OnevsOne_modified) <- sapply(
        lapply(
            strsplit(names(OnevsOne_modified), "\\.log"),
            rev
        ),
        paste,
        collapse = "_"
    )
    names(OnevsOne_modified) <- gsub(names(OnevsOne_modified), pattern = "FC", replacement = "t_stat")
    names(OnevsOne_modified) <- gsub(names(OnevsOne_modified), pattern = "\\.p\\.value", replacement = "p_value")
    names(OnevsOne_modified) <- gsub(names(OnevsOne_modified), pattern = "\\.FDR", replacement = "fdr")
    OnevsOne_modified$ensembl <- rownames(OnevsOne_modified)
    rownames(OnevsOne_modified) <- NULL

    ## Add names from mgi data base
    modeling_result_1vs1 <- add_gene_names(OnevsOne_modified)

    return(modeling_result_1vs1)
}



####################### Add gene symbols from ensemblID #######################

add_gene_names <- function(reshaped_df) {
    mart <- useDataset("mmusculus_gene_ensembl", useMart("ensembl"))
    genes <- reshaped_df$ensembl
    gene_list <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", "mgi_symbol"), values = genes, mart = mart)
    reshaped_df <- merge(reshaped_df, gene_list, by.x = "ensembl", by.y = "ensembl_gene_id", all.x = TRUE)
    l_names <- length(colnames(reshaped_df))
    reshaped_df <- reshaped_df[, c(2:(l_names - 1), 1, l_names)]
    reshaped_df <- dplyr::rename(reshaped_df, gene = mgi_symbol)

    return(reshaped_df)
}
