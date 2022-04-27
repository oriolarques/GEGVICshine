s_gsea <- function(annot_res,
                   gmt,
                   gsea_pvalue = 0.2) {

    # 0. If samples come from mouse (by the presence of one extra column)
    if('human_ortholog' %in% colnames(annot_res) == TRUE){
        # By each data frame in annot_res
        annot_res <- annot_res %>%
            # Substitute mouse gene symbol with the human homolog symbol
            dplyr::mutate(hgnc_symbol = human_ortholog) %>%
            # Filter those missing gene symbols
            dplyr::filter(hgnc_symbol != '') %>%
            # Remove duplicated genes
            dplyr::distinct(hgnc_symbol, .keep_all = TRUE)
    }
        
    # 1. Obtain geneLists
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

    
    # 2. Perform GSEA
    ## Create an empty list to store GSEA results
    temp_gsea <- NULL
    
    ## Execute GSEA
    ### Iterate over the annotated results list
    temp_gsea <- clusterProfiler::GSEA(geneList = temp_gs,
                      TERM2GENE = gmt,
                      pvalueCutoff = gsea_pvalue)
            
    # Check if GSEA result is empty and print a message
    if (nrow(temp_gsea@result) == 0) {
        gsea <- list(no_genesets = temp_gsea@result)
    } else {
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
            dplyr::select(.data$ID, .data$NES, .data$p.adjust, 
                          .data$leading_edge, .data$core_enrichment)
        
        gs.cl <- GSEAmining::gm_clust(df = gs.filt)
        
        gsea <- list(gsea_result = temp_gsea@result,
                     gs.filt = gs.filt,
                     gs.cl = gs.cl)
    }
    
    return(gsea)
    
}
