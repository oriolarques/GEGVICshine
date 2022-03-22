s_deconv <- function(gene_expression,
                     indications = NULL,
                     cibersort_R = NULL,
                     cibersort_LM22 = NULL,
                     tumor = TRUE,
                     rmgenes = NULL,
                     scale_mrna = TRUE,
                     expected_cell_types = NULL) {
    
    library(immunedeconv)  # load immunedeconv to avoid Error: object 'xCell.data' not found
    
    
    # Define immunedeconv methods except CIBERSORT
    idc_methods <- c('quantiseq', 'timer',
                     'mcp_counter', 'xcell', 'epic')
    
    # Create necessary objects
    temp_df <- NULL # Data frame to store the predictions for each method
    ic_pred <- NULL # Final data frame to store all predictions
    
    # Iterate over the methods vector
    for (i in seq_along(idc_methods)) {
        # Obtain predictions for each method
        temp_df <- immunedeconv::deconvolute(gene_expression = gene_expression,
                                             method = idc_methods[i],
                                             indications = indications,
                                             tumor = tumor)
        # Add the method name to the cell_type column
        temp_df <- temp_df %>%
            dplyr::mutate(cell_type = paste0(cell_type, '_', idc_methods[i]))
        
        # Store predictions in ic_pred object
        ## In the first iteration ic_pred = temp_df
        if (i == 1) {
            ic_pred <- temp_df
            ## Afterwards add rows to ic_pred object
        } else {
            ic_pred <- bind_rows(ic_pred, temp_df)
        }
        
    }
    
    # If the user has access to CIBERSOFT
    if (is.null(cibersort_R) == FALSE){
        # Set the path to CIBERSORT and its matrix
        set_cibersort_binary(cibersort_R)
        set_cibersort_mat(cibersort_LM22)
        
        # Obtain predictions for CIBERSORT
        temp_df <- immunedeconv::deconvolute(gene_expression = gene_expression,
                                             method = 'cibersort')
        # Add the method name to the cell_type column
        temp_df <- temp_df %>%
            dplyr::mutate(cell_type = paste0(cell_type, '_CIBERSORT'))
        
        # Add predictions in ic_pred object
        ic_pred <- bind_rows(ic_pred, temp_df)
        
        # Obtain predictions for CIBERSORT-ABS
        temp_df <- immunedeconv::deconvolute(gene_expression = gene_expression,
                                             method = 'cibersort_abs')
        # Add the method name to the cell_type column
        temp_df <- temp_df %>%
            dplyr::mutate(cell_type = paste0(cell_type, '_CIBERSORT-ABS'))
        
        # Add predictions in ic_pred object
        ic_pred <- bind_rows(ic_pred, temp_df)
        
    }
    
    # Separate method name from cell_type
    ic_pred <- ic_pred %>%
        tidyr::separate(cell_type, into = c('cell_type', 'method'), sep = '_') %>%
        dplyr::mutate(method = toupper(method)) %>%
        # Set MCP to MCP_COUNTER
        dplyr::mutate(method = ifelse(method == 'MCP', 'MCP_COUNTER', method))
    
    return(ic_pred)
    
}
