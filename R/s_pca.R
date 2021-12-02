s_pca <- function(counts,
                   genes_id,
                   metadata,
                   response,
                   design,
                   colors = c('black', 'orange')) {
    
    # Preprocess counts data
    counts <- preprocess_ge_counts(counts = counts,
                                   genes_id = genes_id)
    
    # Preprocess metadata
    metadata <- preprocess_ge_meta(metadata = metadata)
    
    # Create DESeq2Dataset object
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                          colData = metadata,
                                          design = formula(paste('~', design,
                                                                 collapse = " ")))
    # Generate  normalized counts
    dds <- DESeq2::estimateSizeFactors(dds)
    
    # Transform normalized counts for data visualization
    vsd <- DESeq2::vst(dds, blind=FALSE)
    
    ## Enquote response variable
    response <- rlang::sym(response)
    # Trick to not use quasiquotation in plotPCA intgroup argument
    ## Get the colname of the response variable so it is quoted
    quoted.resp <- metadata %>%
        dplyr::select(!!response) %>%
        colnames(.)
    
    # plot PCA
    DESeq2::plotPCA(vsd, intgroup = quoted.resp) +
        scale_color_manual(values = colors) +
        labs(title = 'Principal Component Analysis') +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
    
}
