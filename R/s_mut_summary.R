s_mut_summary <- function(muts,
                          metadata,
                          top_genes = 10){
    
    # Process input as MAF file -----------------------------------------------
    maf <- maftools::read.maf(maf = muts,
                              clinicalData = metadata %>%
                                  dplyr::rename('Tumor_Sample_Barcode' = 'Samples'))
    
    # Plot MAF Summary --------------------------------------------------------
    maftools::plotmafSummary(maf = maf,
                             addStat = 'median',
                             titvRaw = FALSE,
                             top = top_genes)
    par(mfrow = c(1,1))

}