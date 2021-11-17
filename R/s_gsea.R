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
        # Get the data frame with gene expression and sort by fold change
        temp_df <- annot_res %>%
            #Filter genes whose log2FoldChange == NA
            filter(!is.na(.data$log2FoldChange)) %>%
            arrange(desc(.data$log2FoldChange))
        
        # Create a vector with the log2FoldChange
        temp_gs <- temp_df$log2FoldChange
        # Get the gene symbols of each log2FoldChange
        names(temp_gs) <- as.character(temp_df$hgnc_symbol)

    
    # 3. Perform GSEA
    ## Create an empty list to store GSEA results
    temp_gsea <- NULL
    
    ## Execute GSEA
    ### Iterate over the annotated results list
        temp_gsea <- GSEA(geneList = temp_gs,
                          TERM2GENE = gmt,
                          pvalueCutoff = gsea_pvalue)
    
    # 3. GSEAmining
        # Filter gene sets to analyse the top ones
        gs.filt <- temp_gsea@result %>%
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
        
        gsea <- list(gsea_result = temp_gsea@result,
                     gs.filt = gs.filt,
                     gs.cl = gs.cl)
        
    
    return(gsea)
    
}
