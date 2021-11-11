library(shiny)
library(GEGVIC)
library(clusterProfiler)
library(dplyr)
library(GSEAmining)
library(ggplot2)
library(ggrepel)

source(file = 'R/s_gsea.R', local = TRUE)
source(file = 'R/s_volcano.R', local = TRUE)

options(shiny.maxRequestSize = 20*1024^2)


#############################################################################        
#############################################################################        
# UI -----------------------------------------------------          
#############################################################################  
#############################################################################        

ui <- fluidPage(
    
    tabsetPanel(
        
        tabPanel(title = 'Parameters',
  #############################################################################        
        # User parameters -----------------------------------------------------          
  #############################################################################        
                 ## Inputs ----------------------------------------------------
                 # Input: RNA-seq raw counts
                 fileInput(inputId = "counts", 
                           label = 'RNA-seq raw counts',
                           accept = c('.csv')),
                 # Input: Metadata
                 fileInput(inputId = "metadata", 
                           label = 'Metadata',
                           accept = c('.csv')),
                    # response
                    ## It will be dynamic UI component
                    selectInput(inputId = 'response',
                                label = 'Select response variable',
                                choices = ""),
                 # Input: Genetic variations
                 fileInput(inputId = "muts", 
                           label = 'Genetic Variations',
                           accept = c('.csv')),
                 # Input: gmt
                 fileInput(inputId = "gmt", 
                           label = 'Gene sets (as .gmt file)',
                           accept = c('.gmt')),
        
                # cibersort
                ## !!!!!!
                ## !!!!!!
                ## !!!!!!
                ## !!!!!!
                 ## Parameters ------------------------------------------------
                 # genes_id
                 selectInput(inputId = 'genes_id', 
                             label = 'Genes ID', 
                             choices = list('Genes ID' = c("hgnc_symbol",
                                                           "entrezgene_id", 
                                                           "ensembl_gene_id"))),
                 # design
                 textInput(inputId = 'design',
                           label = 'Design formula', 
                           placeholder = 'Cell + Treatment + Cell:Treatment'),                
                 # colors
                 textInput(inputId = 'colors',
                           label = 'Colors: Indicate the color for each sample group
                           separated by commas',
                           placeholder = 'black, orange'), 
                 # ref_level
                 textInput(inputId = 'ref_level',
                           label = 'Reference level: Name of the grouping variable',
                           placeholder = 'Non_Responders'), 
                 # shrink
                 selectInput(inputId = 'shrink', 
                             label = 'Shrinkage method', 
                             choices = list('Shrinkage' = c("apeglm", "ashr", 
                                                            "normal", "none")),
                             selected = 'apeglm'),
                 # Input: biomart
                 selectInput(inputId = "biomart", 
                             label = 'BiomaRt database',
                             choices = list('Genome_version' = c('ensembl_biomart_GRCh38_p13',
                                                                 'ensembl_biomart_GRCh37'))),
                 # fold_change
                 numericInput(inputId = 'fold_change', 
                              label = 'Fold Change', 
                              value = 2, 
                              min = 0, 
                              step = 0.5),
                 # p.adj
                 numericInput(inputId = 'p.adj', 
                              label = 'Adjusted p-value for gene expression data', 
                              value = 0.05, 
                              min = 0, 
                              max = 1,
                              step = 0.1),
                 # gsea_pvalue
                 numericInput(inputId = 'gsea_pvalue', 
                              label = 'Adjusted p-value for GSEA', 
                              value = 0.2, 
                              min = 0, 
                              max = 1,
                              step = 0.1),
                 # indications
                 selectInput(inputId = 'indications', 
                             label = 'Cancer type method', 
                             choices = list('Cancer type' = c("kich", "blca", "brca",
                                                              "cesc", "gbm", "hnsc",
                                                              "kirp", "lgg", "lihc",
                                                              "luad", "lusc", "prad",
                                                              "sarc", "pcpg", "paad",
                                                              "tgct", "ucec", "ov",
                                                              "skcm", "dlbc", "kirc", 
                                                              "acc", "meso", "thca",
                                                              "uvm", "ucs", "thym",
                                                              "esca", "stad", "read",
                                                              "coad", "chol"))),
                 # compare
                 selectInput(inputId = 'compare',
                             label = 'Select means comparison method',
                             choices = list('Comparison' = c('t.test',
                                                             'wilcox.test',
                                                             'anova',
                                                             'kruskal.test'))),
                 # p_label
                 selectInput(inputId = 'compare',
                             label = 'Select means comparison method',
                             choices = list('Comparison' = c('p.format',
                                                             'p.signif'))),
                 # top_genes
                 numericInput(inputId = 'top_genes', 
                              label = 'Number of genes for oncoplot', 
                              value = 10, 
                              min = 1, 
                              max = 50,
                              step = 1),
                 # gbuild
                 selectInput(inputId = 'gbuild',
                             label = 'Select genomic build',
                             choices = list('Comparison' = c('BSgenome.Hsapiens.UCSC.hg19',
                                                             'BSgenome.Hsapiens.UCSC.hg38'))),
                 # mut_sigs
                 selectInput(inputId = 'gbuild',
                             label = 'Select mutational signatures matrix',
                             choices = list('Comparison' = c('COSMIC_v2_SBS_GRCh37',
                                                             'COSMIC_v2_SBS_GRCh38',
                                                             'COSMIC_v3.2_SBS_GRCh37',
                                                             'COSMIC_v3.2_SBS_GRCh38',
                                                             'COSMIC_v3.2_DBS_GRCh37',
                                                             'COSMIC_v3.2_DBS_GRCh38'))),
            
            actionButton(inputId = "Execute", label = 'Execute')
                 ),
        
  #############################################################################        
        # GE_module -----------------------------------------------------         
  #############################################################################        
        tabPanel(title = 'GE_module',
            # PCA
            wellPanel(
              
              tags$h1('PCA'),
              
              #Plot PCA
              plotOutput(outputId = 'pca')
            ),
            
            # Differential expressed genes table
            tags$h1('Differentially Expressed genes'),
            
            tabsetPanel(id = 'de_tables')#,
            
            # Volcano plot
            #wellPanel(
              #Plot Volcano
            #  plotOutput(outputId = 'volcano')
            #),
            
            # GSEA
            #wellPanel(
              #Plot gsea cluster
            #  plotOutput(outputId = 'gsea_clust'),
              #Plot gsea wordclouds
           #   plotOutput(outputId = 'gsea_word')
            #)
            
        ),

  #############################################################################        
        # GV_module -----------------------------------------------------         
  #############################################################################        
        tabPanel(title = 'GV_module',
                 ),

  #############################################################################        
        # IC_module -----------------------------------------------------         
  #############################################################################        
        tabPanel(title = 'IC_module',
                
                 )

        
    )
    

)

#############################################################################        
#############################################################################        
# SERVER -----------------------------------------------------          
############################################################################# 
#############################################################################        
server <- function(input, output) {
    
    #############################################################################        
    # Automatic detection of colnames in the metadata---------------------------         
    #############################################################################       
    
    # Create reactives values to use later on
    reactives <- reactiveValues(
        temp_meta = NULL
    )
    
    # React to metadata upload
    observeEvent(input$metadata, {
        # Make it empty if no metadata is loaded
        if (is.null(input$metadata))
            return()
        # Ensure the file is uploaded and available
        req(input$metadata)
        # Read the name in a new variable 
        file <- input$metadata
        # Read metdata in the reactive object already created
        reactives$temp_meta <- read.csv(file = file$datapath, 
                                       header = TRUE, sep = ',')
        # In case the sep value is ';'
        if (ncol(reactives$temp_meta) == 1) {
            reactives$temp_meta <- read.csv(file = file$datapath, 
                     header = TRUE, sep = ';')
        }
        # Update the selection input for response
        updateSelectInput(inputId = 'response',
                          label = 'Select response variable', 
                          choices  = colnames(reactives$temp_meta))
    })
    
   
    #############################################################################        
    # Click button ---------------------------         
    #############################################################################       
    
    observeEvent(input$Execute, {

        #######################################################################        
        # Read and process inputs ---------------------------         
        #######################################################################  
        
        # Read the counts input -----------------------------------------------
        counts <- reactive({
            # Ensure the file is uploaded and available
            req(input$counts)
            # Read the name in a new variable 
            file <- input$counts
            # Get the datapath
            ext <- tools::file_ext(file$datapath)
            # Validate that the file is a .csv
            validate(need(ext == "csv", "Please upload a csv file"))
            # Create a temporal object to store the data
            temp_df <- read.csv(file = file$datapath, 
                                header = TRUE, sep = ',')
            # if separator is ','
            if (ncol(temp_df) != 1) {
                read.csv(file = file$datapath, 
                         header = TRUE, sep = ',')
                # if not, it is ';'
            } else {
                read.csv(file = file$datapath, 
                         header = TRUE, sep = ';')
            }
            
        })
        
        # Read the metadata input ---------------------------------------------
        metadata <- reactive({
            # Ensure the file is uploaded and available
            req(input$metadata)
            # Read the name in a new variable 
            file <- input$metadata
            # Get the datapath
            ext <- tools::file_ext(file$datapath)
            # Validate that the file is a .csv
            validate(need(ext == "csv", "Please upload a csv file"))
            # Create a temporal object to store the data
            temp_df <- read.csv(file = file$datapath, 
                                header = TRUE, sep = ',')
            # if separator is ','
            if (ncol(temp_df) != 1) {
                read.csv(file = file$datapath, 
                         header = TRUE, sep = ',')
                # if not, it is ';'
            } else {
                read.csv(file = file$datapath, 
                         header = TRUE, sep = ';')
            }
            
        })
        
        # Read genetic variations input ---------------------------------------
        muts <- reactive({
          # Ensure the file is uploaded and available
          req(input$muts)
          # Read the name in a new variable 
          file <- input$muts
          # Get the datapath
          ext <- tools::file_ext(file$datapath)
          # Validate that the file is a .csv
          validate(need(ext == "csv", "Please upload a csv file"))
          # Create a temporal object to store the data
          temp_df <- read.csv(file = file$datapath, 
                              header = TRUE, sep = ',')
          # if separator is ','
          if (ncol(temp_df) != 1) {
            read.csv(file = file$datapath, 
                     header = TRUE, sep = ',')
            # if not, it is ';'
          } else {
            read.csv(file = file$datapath, 
                     header = TRUE, sep = ';')
          }
          
        })
        
        # Read gmt input ------------------------------------------------------
        gmt <- reactive({
          # Ensure the file is uploaded and available
          req(input$gmt)
          # Read the name in a new variable 
          file <- input$gmt
          # Get the datapath
          ext <- tools::file_ext(file$datapath)
          # Validate that the file is a .gmt
          validate(need(ext == "gmt", "Please upload a gmt file"))
          
          clusterProfiler::read.gmt(file$datapath)
          
          
        })
        
        
        
        
        # Colors: Convert colors input from text to vector --------------------
        colors <- reactive({
            unlist(strsplit(input$colors,split = ','))
        })
        
        # reference level  ----------------------------------------------------
        ref_level <- reactive({
            c(input$response, input$ref_level)
        })
        
        
        # reference level  ----------------------------------------------------
        ref_level <- reactive({
          c(input$response, input$ref_level)
        })
        
        
        # biomart
        biomart <- reactive({
          if(input$biomart == 'ensembl_biomart_GRCh38_p13'){
            GEGVIC::ensembl_biomart_GRCh38_p13
          } else {
            GEGVIC::ensembl_biomart_GRCh37
          }
          
          
        })
        
        
        # indications  --------------------------------------------------------
        
        
        
        
        #######################################################################        
        # Generate outputs ---------------------------         
        #######################################################################  
        
        # GE_module -----------------------------------------------------------  
        
        # Output: PCA
        output$pca <- renderPlot({
            GEGVIC::ge_pca(counts = counts(),
                           genes_id = input$genes_id,
                           metadata = metadata(),
                           design = input$design,
                           colors = colors())
        })
        
        
        # Differential gene expression calculation
        results.dds <- GEGVIC::ge_diff_exp(counts = counts(),
                                           genes_id = input$genes_id,
                                           metadata = metadata(),
                                           design = input$design,
                                           ref_level = ref_level(),
                                           shrink = input$shrink)
        # Annotate gene results
        annot.res <- GEGVIC::ge_annot(results_dds = results.dds,
                                      genes_id = input$genes_id,
                                      biomart = biomart())
        
        
        # Output: Differential gene expression table
        #for(i in seq_along(annot.res)){
        #  n <- names(annot.res)[i]
        #  appendTab(inputId = 'de_tables',
        #            tabPanel(title = n,
        #                     tags$h4(paste(n)),
        #                     dataTableOutput(outputId = n)), 
        #            select = TRUE)
        #  output[[n]] <- renderDataTable({
        #    annot.res[[i]]
        #  })
        #}
        
        lapply(seq_along(annot.res), function(x){
          
          n <- names(annot.res)[x]
          
          appendTab(inputId = 'de_tables',
                    
                    tabPanel(title = n,
                             tags$h4(paste(n)),
                             # Differential gene expression table
                             wellPanel(
                               dataTableOutput(outputId = paste0(n, '_table'))
                               ),
                             # Volcano plot
                             wellPanel(
                               #Plot Volcano
                               plotOutput(outputId = paste0(n, '_volcano'))
                               ),
                             # GSEA
                             wellPanel(
                               #Plot gsea cluster
                               plotOutput(outputId = paste0(n, '_gsea_clust')),
                               #Plot gsea wordclouds
                               plotOutput(outputId = paste0(n, '_gsea_word'))
                               )),
                    select = TRUE
                  )
          
              # Outpot: Differential gene expression table
              output[[paste0(n, '_table')]] <- renderDataTable({
                annot.res[[x]]
              })
              
              # Output: Volcano Plot
              output[[paste0(n, '_volcano')]] <- renderPlot({
                s_volcano(annot_res = annot.res[[x]],
                          fold_change = input$fold_change,
                          p.adj = input$p.adj)
                
              })
              
              # Output: GSEA
              gsea <- s_gsea(annot_res = annot.res[[x]],
                             gmt = gmt(),
                             gsea_pvalue = input$gsea_pvalue)
              
              output[[paste0(n, '_gsea_clust')]] <- renderPlot({
                
                GSEAmining::gm_dendplot(df = gsea$gs.filt,
                                        hc = gsea$gs.cl)
              })
              
              output[[paste0(n, '_gsea_word')]] <- renderPlot({
                
                GSEAmining::gm_enrichterms(df = gsea$gs.filt,
                                           hc = gsea$gs.cl)
              })
        })
          
          
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
