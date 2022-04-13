s_single <- function(counts,
                     metadata,
                     genes_id,
                     response,
                     design,
                     biomart,
                     gsva_gmt,
                     method,
                     colors = c('black', 'orange'),
                     row.names = TRUE,
                     col.names = TRUE) {
    
    # Enquote response variable
    response <- rlang::enquo(response)
    
    # Preprocess counts data
    counts <- preprocess_ge_counts(counts = counts,
                                   genes_id = genes_id)
    
    # Preprocess metadata
    metadata <- preprocess_ge_meta(metadata = metadata,
                                   counts = counts) %>%
        dplyr::select(!!response)
    
    # Create DESeq2Dataset object
    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                          colData = metadata,
                                          design = formula(paste('~', design,
                                                                 collapse = " ")))
    # Generate  normalized counts
    dds <- DESeq2::estimateSizeFactors(dds)
    
    # Transform normalized counts for data visualization
    vsd <- DESeq2::vst(dds, blind=FALSE)
    
    # Get expression data
    exprs.mat <- SummarizedExperiment::assay(vsd)
    
    # Annotate results
    ## First evaluate if genes are annotated already as hgnc_symbol
    if (genes_id == 'hgnc_symbol'){
        # If so:
        exprs.mat.annot <- as.data.frame(exprs.mat) %>%
            # Rownames to column with the name indicated in the genes_id parameter
            tibble::rownames_to_column(genes_id) %>%
            # Filter those missing gene symbols
            dplyr::filter(hgnc_symbol != '') %>%
            # Remove duplicated genes
            dplyr::distinct(hgnc_symbol, .keep_all = TRUE) %>%
            # Return hgnc_symbol column as rownames
            tibble::column_to_rownames('hgnc_symbol')
        
    } else {
        # If genes are identified as entrezgene_id or ensembl_gene_id:
        exprs.mat.annot <- as.data.frame(exprs.mat) %>%
            # Rownames to column with the name indicated in the genes_id parameter
            tibble::rownames_to_column(genes_id) %>%
            # Join the data frame with the GRCh38_p13 biomaRt table stored in data
            dplyr::inner_join(x = .,
                              y = biomart %>%
                                  mutate(entrezgene_id = as.character(entrezgene_id)),
                              by = genes_id) %>%
            # From all annotation columns keep only hgnc symbol column
            dplyr::select(hgnc_symbol,
                          everything(),
                          -c(entrezgene_id, ensembl_gene_id,
                             transcript_length, refseq_mrna)) %>%
            # Filter those missing gene symbols
            dplyr::filter(hgnc_symbol != '') %>%
            # Remove duplicated genes
            dplyr::distinct(hgnc_symbol, .keep_all = TRUE) %>%
            # Return hgnc_symbol column as rownames
            tibble::column_to_rownames('hgnc_symbol')
    }
    
    # If samples come from mouse (by the presence of one extra column)
    if('human_ortholog' %in% colnames(exprs.mat.annot) == TRUE){
        exprs.mat.annot <- exprs.mat.annot %>%
            # Rownames to column with the name indicated in the genes_id parameter
            tibble::rownames_to_column('hgnc_symbol') %>%
            # Substitute mouse gene symbol with the human homolog symbol
            dplyr::mutate(hgnc_symbol = human_ortholog) %>%
            # Filter those missing gene symbols
            dplyr::filter(hgnc_symbol != '') %>%
            # Remove duplicated genes
            dplyr::distinct(hgnc_symbol, .keep_all = TRUE) %>%
            # Remove ortholog column
            dplyr::select(-human_ortholog) %>%
            # Return hgnc_symbol column as rownames
            tibble::column_to_rownames('hgnc_symbol')
    }
    
    # Read genesets
    if(gsva_gmt == 'Hallmark'){
        gmt <- hallmark.gmt
    } else {
        gmt <- GSEABase::getGmt(gsva_gmt)
    }
    
    # Calculate GSVA/ssGSEA
    gsva_temp <- GSVA::gsva(expr = as.matrix(exprs.mat.annot),
                            gset.idx.list = gmt,
                            method = method)
    
    # Plot Heatmap
    ## Define response level group colors in a list
    temp_color <- colors
    
    resp.levels <- metadata %>%
        dplyr::select(!!response) %>%
        dplyr::pull(!!response)
    
    names(temp_color) <- levels(resp.levels)
    
    # Get the quoted name of the response variable
    quoted.resp <- metadata %>%
        dplyr::select(!!response) %>%
        colnames(.)
    
    pheat.anno.color <- list(temp_color)
    # Name the list
    names(pheat.anno.color) <- quoted.resp
    
    # Plot heatmap
    heat.map <- pheatmap(gsva_temp,
                         scale = 'row',
                         #cutree_rows = 2,
                         #cutree_cols = 2,
                         color = colorRampPalette(c('blue', 'grey', 'red'))(10),
                         show_rownames = row.names,
                         show_colnames = col.names,
                         annotation_col = metadata,
                         annotation_colors = pheat.anno.color,
                         annotation_names_col = FALSE,
                         clustering_method = 'ward.D',
                         main = paste0('Samples clustering by ', toupper(method)))
        
    # Return results table
    gsva_out <- list(gsva_table = gsva_temp,
                     heatmap = ggplotify::as.ggplot(heat.map))
    return(gsva_out)
    
}