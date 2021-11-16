# Required packages
library(shiny)
library(GEGVIC)
library(clusterProfiler)
library(dplyr)
library(tidyr)
library(tibble)
library(GSEAmining)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggplotify)
library(rlang)
library(pheatmap)
library(deconstructSigs)
library(shinyFiles)
library(shinythemes)

# Functions inside the App
source(file = 'R/s_gsea.R', local = TRUE)
source(file = 'R/s_volcano.R', local = TRUE)
source(file = 'R/s_mut_summary.R', local = TRUE)
source(file = 'R/s_oncoplot.R', local = TRUE)
source(file = 'R/s_mut_load.R', local = TRUE)
source(file = 'R/s_mut_signatures.R', local = TRUE)
source(file = 'R/s_plot_comp_samples.R', local = TRUE)
source(file = 'R/s_plot_comp_celltypes.R', local = TRUE)
source(file = 'R/s_score.R', local = TRUE)

# Increment input file size available for gmt files.
options(shiny.maxRequestSize = 20*1024^2)


#############################################################################        
#############################################################################        
# UI -----------------------------------------------------          
#############################################################################  
#############################################################################        

ui <- navbarPage(
    title = 'GEGVICshine',
    
    theme = shinytheme('cosmo'),
    
    tabPanel(title = 'Parameters',
  #############################################################################        
        # User parameters -----------------------------------------------------          
  #############################################################################    
      tabsetPanel(
        
        id = 'parameters_panel',
        
        ## Modules selection ----------------------------------------------
        tabPanel(
          title = '1. Modules selection',
          
          wellPanel(
            tags$h3('Select all the modules you want to run'),
            tags$hr(),
            checkboxInput(inputId = 'ge_module', 
                          label = 'Analyse Gene Expression (GE_module)', 
                          value = FALSE),
            checkboxInput(inputId = 'gv_module', 
                          label = 'Analyse Genetic Variatons (GV_module)', 
                          value = FALSE),
            checkboxInput(inputId = 'ic_module', 
                          label = 'Analyse Immune Composition (IC_module)', 
                          value = FALSE)
            ),
          actionButton(inputId = 'to_input', 
                       label = 'Next section',
                       class = 'btn-primary')
        ),
        
        ## Inputs ----------------------------------------------------
        # Input: RNA-seq raw counts
        tabPanel(
          title = '2. Inputs',
          wellPanel(
            tags$h3('Upload necessary inputs'),
            tags$hr(),
            fileInput(inputId = "counts", 
                      label = 'RNA-seq raw counts (GE, IC)',
                      accept = c('.csv')),
            # Input: Metadata
            fileInput(inputId = "metadata", 
                      label = 'Metadata (GE, GV, IC)',
                      accept = c('.csv')),
              # response
              ## It will be dynamic UI component
              selectInput(inputId = 'response',
                          label = 'Select response variable (GV, IC)',
                          choices = ""),
            # Input: Genetic variations
            fileInput(inputId = "muts", 
                      label = 'Genetic Variations (GV)',
                      accept = c('.csv')),
            # Input: gmt
            fileInput(inputId = "gmt", 
                      label = 'Gene sets (as .gmt file) (GE)',
                      accept = c('.gmt')),
        
                # cibersort
                tags$p(tags$strong('Select folder in your computer containing:')),
                tags$p(tags$strong('CIBERSORT.R and LM22.txt files (IC)')),
                ## Create a button
                shinyDirButton(id = 'cibersort', 
                               label = 'Browse folders', 
                               title = 'Select folder'),
                ## Generate output to show the user that the path has been obtained
                textOutput(outputId = 'cibersort.path'),
                tags$p()
            ),
          actionButton(inputId = 'to_parameters', 
                       label = 'Next section',
                       class = 'btn-primary')
        ),
  
  
            ## Parameters -----------------------------------------------------
          tabPanel(
            title = '3. Parameters',
            wellPanel(
              tags$h3('Define necessary parameters'),
              tags$hr(),
                 # genes_id
                 selectInput(inputId = 'genes_id', 
                             label = 'Genes ID (GE, IC)', 
                             choices = list('Genes ID' = c("hgnc_symbol",
                                                           "entrezgene_id", 
                                                           "ensembl_gene_id"))),
                 # design
                 textInput(inputId = 'design',
                           label = 'Design formula (GE)', 
                           placeholder = 'Cell + Treatment + Cell:Treatment'),                
                 # colors
                 textInput(inputId = 'colors',
                           label = 'Colors: Indicate the color for each sample group
                           separated by commas (GE, GV, IC)',
                           placeholder = 'black, orange'), 
                 # ref_level
                 textInput(inputId = 'ref_level',
                           label = 'Reference level: Name of the grouping variable (GE)',
                           placeholder = 'Non_Responders'), 
                 # shrink
                 selectInput(inputId = 'shrink', 
                             label = 'Shrinkage method (GE)', 
                             choices = list('Shrinkage' = c("apeglm", "ashr", 
                                                            "normal", "none")),
                             selected = 'apeglm'),
                 # Input: biomart
                 selectInput(inputId = "biomart", 
                             label = 'BiomaRt database (GE, IC)',
                             choices = list('Genome_version' = c('ensembl_biomart_GRCh37',
                                                                 'ensembl_biomart_GRCh38_p13',
                                                                 'ensembl_biomart_GRCm38_p6',
                                                                 'ensembl_biomart_GRCm39'))),
                 # fold_change
                 numericInput(inputId = 'fold_change', 
                              label = 'Fold Change (GE)', 
                              value = 2, 
                              min = 0, 
                              step = 0.5),
                 # p.adj
                 numericInput(inputId = 'p.adj', 
                              label = 'Adjusted p-value for gene expression data (GE)', 
                              value = 0.05, 
                              min = 0, 
                              max = 1,
                              step = 0.05),
                 # gsea_pvalue
                 numericInput(inputId = 'gsea_pvalue', 
                              label = 'Adjusted p-value for GSEA (GE)', 
                              value = 0.2, 
                              min = 0, 
                              max = 1,
                              step = 0.05),
                 # indications
                 selectInput(inputId = 'indications', 
                             label = 'Cancer types: TCGA Study Abbreviations. (IC)', 
                             choices = list('Cancer type' = sort(c("kich", "blca", "brca",
                                                              "cesc", "gbm", "hnsc",
                                                              "kirp", "lgg", "lihc",
                                                              "luad", "lusc", "prad",
                                                              "sarc", "pcpg", "paad",
                                                              "tgct", "ucec", "ov",
                                                              "skcm", "dlbc", "kirc", 
                                                              "acc", "meso", "thca",
                                                              "uvm", "ucs", "thym",
                                                              "esca", "stad", "read",
                                                              "coad", "chol")))),
                 # compare
                 selectInput(inputId = 'compare',
                             label = 'Select means comparison method (GV, IC)',
                             choices = list('Comparison' = c('t.test',
                                                             'wilcox.test',
                                                             'anova',
                                                             'kruskal.test'))),
                 # p_label
                 selectInput(inputId = 'p_label',
                             label = 'Select means comparison method (GV, IC)',
                             choices = list('Comparison' = c('p.format',
                                                             'p.signif'))),
                 # top_genes
                 numericInput(inputId = 'top_genes', 
                              label = 'Number of genes for oncoplot (GV)', 
                              value = 10, 
                              min = 1, 
                              max = 50,
                              step = 1),
                 # gbuild
                 selectInput(inputId = 'gbuild',
                             label = 'Select genomic build (GV)',
                             choices = list('Comparison' = c('BSgenome.Hsapiens.UCSC.hg19',
                                                             'BSgenome.Hsapiens.UCSC.hg38',
                                                             'BSgenome.Mmusculus.UCSC.mm10',
                                                             'BSgenome.Mmusculus.UCSC.mm39'))),
                 # mut_sigs
                 selectInput(inputId = 'mut_sigs',
                             label = 'Select mutational signatures matrix (GV)',
                             choices = list('Mutational signature matrix' = 
                                                           c('COSMIC_v2_SBS_GRCh37',
                                                             'COSMIC_v2_SBS_GRCh38',
                                                             'COSMIC_v3.2_SBS_GRCh37',
                                                             'COSMIC_v3.2_SBS_GRCh38',
                                                             'COSMIC_v3.2_DBS_GRCh37',
                                                             'COSMIC_v3.2_DBS_GRCh38',
                                                             'COSMIC_v2_SBS_mm9',
                                                             'COSMIC_v2_SBS_mm10',
                                                             'COSMIC_v3.2_SBS_mm9',
                                                             'COSMIC_v3.2_SBS_mm10',
                                                             'COSMIC_v3.2_DBS_mm9',
                                                             'COSMIC_v3.2_DBS_mm10')))
                 ),
            actionButton(inputId = 'to_execute', 
                         label = 'Next section',
                         class = 'btn-primary')
            
            ),
          
          
              ## Execution button -----------------------------------------------
          tabPanel(
            title = '4. Run GEGVIC',
            wellPanel(
              tags$h4('Run GEGVIC if you have filled all inputs and parameters'),
              tags$strong('Check the results visiting each module section using 
                          top level navigation bar.'),
              tags$hr(),
              actionButton(inputId = "Execute", 
                           label = 'Execute')
            )
          ) 
                  
      )
    ), # End Parameters
        
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
            
            # Dynamic tab to show per comparison a tab with:
            ## 1. Differential expressed genes table
            ## 2. Volcano plot
            ## 3. GSEA (cluster + wordclouds)
            tags$h1('Differentially Expressed genes'),
            
            tabsetPanel(id = 'de_tables')
            
        ), # End GE_module

  #############################################################################        
        # GV_module -----------------------------------------------------         
  #############################################################################        
        tabPanel(title = 'GV_module',
                 
                 wellPanel(
                   tags$h1('Mutations summary'),
                   plotOutput('mut_sum')
                   
                 ),
                 
                 wellPanel(
                   tags$h1('Oncoplot'),
                   plotOutput('oncoplot')
                 ),
                 
                 wellPanel(
                   tags$h1('Mutational Load'),
                   plotOutput('mut_load')
                 ),
                 
                 wellPanel(
                   tags$h1('Mutational Signatures'),
                   plotOutput('mut_sigs_bar'),
                   plotOutput('mut_sigs_heat')
                 )
                 
                 
                 
                 ), # End GV_module

  #############################################################################        
        # IC_module -----------------------------------------------------         
  #############################################################################        
        tabPanel(title = 'IC_module',
                 
                 wellPanel(
                   tags$h1('Immune composition by sample'),
                   plotOutput('ic_samples')
                 ),
                 
                 wellPanel(
                   tags$h1('Immune composition by cell type'),
                   plotOutput('ic_type')
                 ),
                 
                 wellPanel(
                   tags$h1('Immune score'),
                   plotOutput('ic_phenogram'),
                   plotOutput('ic_score')
                 )
                
              ) # End IC_module

)

#############################################################################        
#############################################################################        
# SERVER -----------------------------------------------------          
############################################################################# 
#############################################################################        
server <- function(input, output, session) {
  
  
  #############################################################################        
  # Buttons to switch sections in parameters tab ---------------------------         
  #############################################################################       
  
  observeEvent(input$to_input, {
    updateTabsetPanel(session, 
                      inputId = "parameters_panel",
                      selected = '2. Inputs')
  })
  
  observeEvent(input$to_parameters, {
    updateTabsetPanel(session, 
                      inputId = "parameters_panel",
                      selected = '3. Parameters')
  })
  
  observeEvent(input$to_execute, {
    updateTabsetPanel(session, 
                      inputId = "parameters_panel",
                      selected = '4. Run GEGVIC')
  })
    
  
  
  
  
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
    
    
    # Read the CIBERSORT.R file -------------------------------------------
    #cibersort <- reactive({
    # Read the name in a new variable 
    # file <- input$cibersort
    # Get the datapath
    #ext <- tools::file_ext(file$datapath)
    #return(ext)
    #})
    
    shinyFiles::shinyDirChoose(input, 
                               id = 'cibersort', 
                               roots = getVolumes()(), 
                               session = session)
    cibersort.path <- reactive({
      paste0(as.character(parseDirPath(roots = getVolumes()()
                                       ,selection = input$cibersort)),
             '/')
    })
    
    output$cibersort.path <- renderText({
      cibersort.path()
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
        
        # cibersort -----------------------------------------------------------
        #ciber.path <- reactive({
        #  as.character(parseDirPath(roots = getVolumes()(),
        #                            selection = input$cibersort))
        #})
        
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
        indications <- reactive({
          
          rep(input$indications, nrow(metadata()))
          
        })
        
        
        
        #######################################################################        
        # Generate outputs ---------------------------         
        #######################################################################  
        
        # GE_module -----------------------------------------------------------  
        
        if(input$ge_module == TRUE) {
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
          
          
          # Output: Dynamic tabs depending on the number of comparisons
          ## Use lapply instead of a loop so shiny doesn't just show the values from last comparison
          lapply(seq_along(annot.res), function(x){
            # Store the names of the comparisons
            n <- names(annot.res)[x]
            ## Create a tab for each comparison with the appropriate name
            appendTab(inputId = 'de_tables',
                      # In each tab there will be three panels
                      tabPanel(title = n,
                               tags$h4(paste(n)),
                               # 1. Differential gene expression table
                               wellPanel(
                                 dataTableOutput(outputId = paste0(n, '_table'))
                                 ),
                               # 2. Volcano plot
                               wellPanel(
                                 #Plot Volcano
                                 plotOutput(outputId = paste0(n, '_volcano'))
                                 ),
                               # 3. GSEA
                               wellPanel(
                                 #Plot gsea cluster
                                 plotOutput(outputId = paste0(n, '_gsea_clust')),
                                 #Plot gsea wordclouds
                                 plotOutput(outputId = paste0(n, '_gsea_word'))
                                 )),
                      select = TRUE
                    )
            
                # Generate Output: Differential gene expression table
                output[[paste0(n, '_table')]] <- renderDataTable({
                  annot.res[[x]]
                })
                
                # Generate Output: Volcano Plot
                output[[paste0(n, '_volcano')]] <- renderPlot({
                  s_volcano(annot_res = annot.res[[x]],
                            fold_change = input$fold_change,
                            p.adj = input$p.adj)
                  
                })
                
                # Generate Output: GSEA
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
        }
        
        # GV_module -----------------------------------------------------------
          
        if(input$gv_module == TRUE){
        
          output$mut_sum <- renderPlot({
            
            s_mut_summary(muts = muts(),
                          metadata = metadata(),
                          top_genes = input$top_genes)
            
          })
          
          output$oncoplot <- renderPlot({
            
            s_oncoplot(muts = muts(),
                       metadata = metadata(),
                       response = input$response,
                       top_genes = input$top_genes,
                       colors = colors())
            
          })
          
          output$mut_load <- renderPlot({
            s_mut_load(muts = muts(),
                       metadata = metadata(), 
                       response = input$response,
                       compare = input$compare,
                       p_label = input$p_label,
                       colors = colors())
            
            
            
          })
          
          
          mut.sigs <- s_mut_signatures(muts = muts(),
                                       metadata = metadata(), 
                                       response = input$response,
                                       gbuild = input$gbuild,
                                       mut_sigs = input$mut_sigs,
                                       colors = colors())
          
          output$mut_sigs_bar <- renderPlot({
            
            mut.sigs$mut_sig_barplot
            
            
          })
          
          output$mut_sigs_heat <- renderPlot({
            
            mut.sigs$mut_sig_heatmap
            
            
          })
          
        }
        
        # IC_module -----------------------------------------------------------  
        
        if(input$ic_module == TRUE){
          
          tpm <- GEGVIC::ic_raw_to_tpm(counts = counts(),
                                       genes_id = input$genes_id,
                                       biomart = biomart())
          
          ic.pred <- GEGVIC::ic_deconv(gene_expression = tpm,
                                       indications = indications(),
                                       cibersort = cibersort.path(),
                                       tumor = TRUE,
                                       rmgenes = NULL,
                                       scale_mrna = TRUE,
                                       expected_cell_types = NULL)
          
          
          output$ic_samples <- renderPlot({
            
            s_plot_comp_samples(df = ic.pred,
                                 metadata = metadata(),
                                 response = input$response,
                                 compare = input$compare,
                                 p_label = input$p_label,
                                 colors = colors())
            
          })
          
          output$ic_type <- renderPlot({
            
            s_plot_comp_celltypes(df = ic.pred,
                                  metadata = metadata(),
                                  response = input$response)
            
          })
          
          
          ic.IPS_IPG <- s_score(tpm = tpm,
                                metadata = metadata(),
                                response = input$response,
                                compare = input$compare,
                                p_label = input$p_label,
                                colors = colors())
          
          
          output$ic_phenogram <- renderPlot({
            ic.IPS_IPG$immunophenoGram
            
          })
          
          output$ic_score<- renderPlot({
            ic.IPS_IPG$immunophenoScore
            
          })
          
          
          
        }
          
    }) # End button
}

# Run the application 
shinyApp(ui = ui, server = server)
