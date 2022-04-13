s_oncoplot <- function(muts,
                       metadata,
                       response,
                       top_genes = 10,
                       col.names = TRUE,
                       colors = c('black' ,'orange')
                       ){
    
    # Transform response to symbol for later use (instead of enquote)
    response <- rlang::sym(response)
    
    # Process input as MAF file -----------------------------------------------
    maf <- maftools::read.maf(maf = muts,
                              clinicalData = metadata %>%
                                  dplyr::rename('Tumor_Sample_Barcode' = 'Samples'))
    
    # Oncoplot ----------------------------------------------------------------

    # Determine the order of the samples to separate between groups in the plot
    samples_order <- metadata %>%
        dplyr::arrange(!!response) %>%
        dplyr::pull(Samples)
    
    # Trick to not use quasiquotation with maftools oncoplot to define annotationColor
    ## Get the colname of the response variable so it is quoted
    quoted.resp <- metadata %>%
        dplyr::select(!!response) %>%
        colnames(.)
    ## Get the levels in the response variable
    resp.levels <- metadata %>%
        dplyr::select(!!response) %>%
        dplyr::pull(!!response)
    ## Associate user input colors to an object
    cols <- colors
    # Name the colors vector with the response variable names
    names(cols) <- unique(resp.levels)
    ## Create a list with the named color vector
    cols.list <- list(cols)
    ## Add the quoted named of the response variable as name of the list
    names(cols.list) <- quoted.resp
    
    # Plot oncoplot
    maftools::oncoplot(maf = maf,
                       top = top_genes,
                       clinicalFeatures = quoted.resp,
                       genes = NULL,
                       sampleOrder = samples_order,
                       showTumorSampleBarcodes = col.names,
                       sortByAnnotation = TRUE,
                       annotationColor = cols.list)
    par(mfrow = c(1,1))
    
    
    
}