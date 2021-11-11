s_gsea <- function(annot_res,
                   gmt,
                   gsea_pvalue = 0.2) {

    # 1. Obtain geneLists
    ## Create a list to store as many geneList as conditions
    geneLists <- list()
    
    ## Create two temporal objects
    temp_df <- NULL
    temp_gs <- NULL
    
    ## Create geneLists
    ### Iterate over the annotated results list
    for (i in seq_along(annot_res)) {
        # Get the data frame with gene expression and sort by fold change
        temp_df <- annot_res[[i]] %>%
            #Filter genes whose log2FoldChange == NA
            filter(!is.na(.data$log2FoldChange)) %>%
            arrange(desc(.data$log2FoldChange))
        
        # Create a vector with the log2FoldChange
        temp_gs <- temp_df$log2FoldChange
        # Get the gene symbols of each log2FoldChange
        names(temp_gs) <- as.character(temp_df$hgnc_symbol)
        
        # Store the geneList
        geneLists[[i]] <- temp_gs
        names(geneLists)[i] <- paste0('geneList_', names(annot_res)[i])
    }
    
    # 3. Perform GSEA
    ## Create an empty list to store GSEA results
    GSEA.res <- list()
    temp_gsea <- NULL
    
    ## Execute GSEA
    ### Iterate over the annotated results list
    for (i in seq_along(annot_res)) {
        temp_gsea <- GSEA(geneList = geneLists[[i]],
                          TERM2GENE = gmt,
                          pvalueCutoff = gsea_pvalue)
        # Delay the next process 0.5 seconds
        Sys.sleep(0.5)
        # Save the GSEA results
        GSEA.res[[i]] <- temp_gsea
        names(GSEA.res)[i] <- paste0('GSEA_', names(annot_res)[i])
    }
    
    # 3. GSEAmining
    ### Iterate over the annotated results list
    for (i in seq_along(annot_res)) {
        # Filter gene sets to analyse the top ones
        gs.filt <- GSEA.res[[i]]@result %>%
            dplyr::arrange(desc(.data$NES)) %>%
            dplyr::mutate(group = ifelse(test = .data$NES > 0,
                                         yes = 'Positive',
                                         no = 'Negative')) %>%
            dplyr::group_by(.data$group) %>%
            dplyr::filter(.data$p.adjust < gsea_pvalue) %>%
            dplyr::top_n(., n = 20, wt = abs(.data$NES)) %>%
            dplyr::ungroup(.) %>%
            dplyr::select(.data$ID, .data$NES,
                          .data$p.adjust, .data$core_enrichment)
        
        gs.cl <- GSEAmining::gm_clust(df = gs.filt)
        
        gsea <- list(gs.filt = gs.filt,
                     gs.cl = gs.cl)
        
        
    }
    
    return(gsea)
    
}