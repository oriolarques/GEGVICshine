s_mut_load <- function(muts,
                        metadata,
                        response,
                        compare = NULL,
                        p_label = 'p.format',
                        colors = c('black', 'orange')){
    
    # Get nonsynonymous mutations
    mut.load <- muts %>%
        dplyr::filter(grepl('Missense_Mutation|Nonsense_Mutation', Variant_Classification)) %>%
        # Group by patient
        dplyr::group_by(Tumor_Sample_Barcode) %>%
        # Sum the number of mutations per patient
        dplyr::summarise(mut_load = n()) %>%
        # Join mutations with metadata to get the response groups
        dplyr::right_join(x = .,
                          y = metadata,
                          by = c('Tumor_Sample_Barcode' = 'Samples'))
    
    
    # Transform response to symbol for later use (instead of enquote)
    response <-  rlang::sym(response)
    
    # Plot mutational load by groups
    p <- ggplot(mut.load, aes(x = !!response,
                              y = mut_load,
                              col = !!response)) +
        # Geometric objects
        geom_violin() +
        geom_boxplot(width = 0.2, outlier.shape = NA) +
        geom_point(alpha = 0.5, position = position_jitter(0.1)) +
        
        # Define colors
        scale_color_manual(values = colors) +
        
        # Title
        ggtitle('Mutational Load:\nBetween groups comparison per population') +
        
        # Themes
        theme_bw() +
        theme(
            plot.title = element_text(size = 15, hjust = 0.5, face = 'bold'),
            axis.text = element_text(size=15, face = "bold"),
            axis.text.x.bottom = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = 'bottom',
            legend.title = element_text(face='bold', size =12),
            legend.text = element_text(size =12)
        )
    
    # Add p-values for comparisons
    if(is.null(compare) == FALSE){
        library('ggplot2') # load to avoid Error in `ggpubr::stat_compare_means()`  could not find function "after_stat"
        
        p <- p +
            ggpubr::stat_compare_means(method = compare,
                                       label = p_label,
                                       label.y.npc = 0.95,
                                       label.x.npc = 0.3,
                                       show.legend = FALSE)
        
    }
    
    # Return the plot
    res.mut.load <- list(mut.load.table = mut.load,
                         mut.load.plot = p)
    return(res.mut.load)
}