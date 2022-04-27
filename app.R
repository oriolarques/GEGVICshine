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
library(patchwork)
library(gridExtra)
library(rlang)
library(pheatmap)
library(deconstructSigs)
library(DT)
library(shinyFiles)
library(shinythemes)
library(tm)

# Functions inside the App
source(file = 'R/s_deconv.R', local = TRUE)
source(file = 'R/s_gsea.R', local = TRUE)
source(file = 'R/s_volcano.R', local = TRUE)
source(file = 'R/s_mut_summary.R', local = TRUE)
source(file = 'R/s_oncoplot.R', local = TRUE)
source(file = 'R/s_mut_load.R', local = TRUE)
source(file = 'R/s_mut_signatures.R', local = TRUE)
source(file = 'R/s_plot_comp_samples.R', local = TRUE)
source(file = 'R/s_plot_comp_celltypes.R', local = TRUE)
source(file = 'R/s_score.R', local = TRUE)
source(file = 'R/s_single.R', local = TRUE)


# Increment input file size available for gmt files.
options(shiny.maxRequestSize = 100*1024^2)


#############################################################################        
#############################################################################        
# UI -----------------------------------------------------          
#############################################################################  
#############################################################################        

ui <- navbarPage(
    title = 'GEGVICshine',
    
    theme = shinytheme('cosmo'),
    
    selected = 'Parameters',
    
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
                       class = 'btn-primary'),
          br(),
          tags$hr(),
          # Download sample data
          tags$h3('Sample Data'),
          downloadButton(outputId = 'down_sample_data', 
                         label = 'Download Sample Data')
         
        ),
        
        ## Inputs ----------------------------------------------------
        tabPanel(
          title = '2. Inputs',
          wellPanel( 
            tags$h3('Upload necessary inputs'),
            tags$hr(),
            # Input: RNA-seq raw counts
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
                          label = 'Select response variable (GE, GV, IC)',
                          choices = ""),
            # Input: Genetic variations
            fileInput(inputId = "muts", 
                      label = 'Genetic Variations (GV)',
                      accept = c('.csv')),
            # Input: gmt
            fileInput(inputId = "gmt", 
                      label = 'Gene sets (as .gmt file) (GE)',
                      accept = c('.gmt'))
            ),
          wellPanel(
            tags$h3('Optional: CIBERSORT'),
            tags$hr(),
            # cibersort
            fileInput(inputId = "cibersort_R", 
                      label = 'CIBERSORT: .R file (IC)',
                      accept = c('.R')),
            fileInput(inputId = "cibersort_LM22", 
                      label = 'CIBERSORT: LM22.txt file (IC)',
                      accept = c('.txt'))
            ),
          actionButton(inputId = 'to_parameters', 
                       label = 'Next section',
                       class = 'btn-primary'),
          
          ## Add javascript code so when opening tab it scrolls to the top
          tags$script("$(document).ready(function () {
                        $('#to_parameters').on('click', function (e) {
                        window.scrollTo(0, 0)
                        });                 
                      });")
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
                             choices = list('Genes ID' = c("",
                                                           "hgnc_symbol",
                                                           "entrezgene_id", 
                                                           "ensembl_gene_id")),
                             selected = ),
                 # design
                 textInput(inputId = 'design',
                           label = 'Design formula (GE)', 
                           placeholder = 'Example: Cell + Treatment + Cell:Treatment'),                
                 # ref_level
                 selectInput(inputId = 'ref_level',
                             label = 'Reference level: Name of the grouping variable (GE)',
                             choices = ""), 
                 # colors
                 textInput(inputId = 'colors',
                           label = 'Colors: Indicate the color for each sample group
                           separated by commas (GE, GV, IC)',
                           placeholder = 'Example: black, orange'),
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
                 # gsva_gmt
                 selectInput(inputId = 'gsva_gmt',
                             label = 'Select genesets for GSVA (GE)',
                             choices = list('Gene set collection' = c('Hallmark',
                                                                      'Same as for GSEA'))),
                 # gsva_method
                 selectInput(inputId = 'gsva_method',
                             label = 'GSVA method (GE)',
                             choices = list('Gene set collection' = c('gsva',
                                                                      'ssgsea',
                                                                      'zscore'))),
                 
                 tags$strong('Advanced plot options'),
                     # gsva_rownames
                     checkboxInput(inputId = 'gsva_rownames', 
                                   label = 'Show row names in GSVA plot (GE)', 
                                   value = TRUE, 
                                   width = NULL),
                     # gsva_colnames
                     checkboxInput(inputId = 'gsva_colnames', 
                                   label = 'Show column names in GSVA plot (GE)', 
                                   value = TRUE, 
                                   width = NULL),
                    # mut_colnames
                    checkboxInput(inputId = 'mut_colnames', 
                                  label = 'Show column names in Oncoplot, Mutational Signature analyses and Immune Cell Composition (GV, IC)', 
                                  value = TRUE, 
                                  width = NULL),
                    # ic_samples_points
                    checkboxInput(inputId = 'ic_samples_points', 
                                  label = 'Show points in immune prediction plots (IC)', 
                                  value = TRUE, 
                                  width = NULL),
                 # indications
                 selectInput(inputId = 'indications', 
                             label = 'Cancer types: TCGA Study Abbreviations (IC)', 
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
              plotOutput(outputId = 'pca', 
                         width = '80%')
            ),
            
            # Dynamic tab to show per comparison a tab with:
            ## 1. Differential expressed genes table
            ## 2. Volcano plot
            ## 3. GSEA (cluster + wordclouds)
            tags$h1('Differentially Expressed Genes'),
            
            tabsetPanel(id = 'de_tables'),
            
            # GSVA
            wellPanel(
              tags$h1('GSVA'),
              ## Heatmap
              tags$h4('Heatmap'),
              plotOutput(outputId = 'gsva',height = '1000px'),
              ## Table with the results
              tags$h4('Table of Gene Sets Enrichment'),
              dataTableOutput(outputId = 'gsva_table'),
              downloadButton(outputId = 'down_gsva_table', 
                             label = 'Download table')
              
            )
            
        ), # End GE_module

  #############################################################################        
        # GV_module -----------------------------------------------------         
  #############################################################################        
        tabPanel(title = 'GV_module',
                 
                 wellPanel(
                   tags$h1('Mutations Summary'),
                   plotOutput(outputId = 'mut_sum', 
                              width = '75%',
                              height = '600px')
                   
                 ),
                 
                 wellPanel(
                   tags$h1('Oncoplot'),
                   plotOutput(outputId = 'oncoplot',
                              width = '75%',
                              height = '600px')
                 ),
                 
                 wellPanel(
                   tags$h1('Mutational Load'),
                   plotOutput(outputId = 'mut_load', 
                              width = '75%',
                              height = '600px'),
                   ## Table with the results
                   tags$h4('Mutational Load Table'),
                   dataTableOutput(outputId = 'mut_load_table'),
                   downloadButton(outputId = 'down_mut_load_table', 
                                  label = 'Download table')
                 ),
                 
                 wellPanel(
                   tags$h1('Mutational Signatures'),
                   plotOutput(outputId = 'mut_sigs_bar', 
                              width = '80%',
                              height = '700px'),
                   plotOutput(outputId = 'mut_sigs_heat', 
                              width = '80%',
                              height = '700px'),
                   ## Table with the results
                   tags$h4('Results Table'),
                   dataTableOutput(outputId = 'mut_sig_table'),
                   downloadButton(outputId = 'down_mut_sig_table', 
                                  label = 'Download table')
                 )
              
               ), # End GV_module

  #############################################################################        
        # IC_module -----------------------------------------------------         
  #############################################################################        
        tabPanel(title = 'IC_module',
                 
                 # Table with the prediction results
                 wellPanel(
                   tags$h1('Summary of Predicted Immune Cell Populations'),
                   dataTableOutput(outputId = 'ic_pred_table'),
                   downloadButton(outputId = 'down_ic_pred_table', 
                                  label = 'Download table')
                 ),
                 # Plot Immune Composition: Cell Types Comparison by Groups
                 wellPanel(
                   tags$h1('Immune Composition: Cell Types Comparison by Groups'),
                   plotOutput(outputId = 'ic_samples',
                              height = '1000px')
                 ),
                 # Plot Immune Composition: Cell Types Comparison within Samples
                 wellPanel(
                   tags$h1('Immune Composition: Cell Types Comparison within Samples'),
                   plotOutput(outputId = 'ic_type',
                              height = '800px')
                 ),
                 # Plot IPG and IPS
                 wellPanel(
                   tags$h1('Immune Score'),
                   #plotOutput(outputId = 'ic_phenogram',
                   #          height = '800px'),
                   plotOutput(outputId = 'ic_score', 
                              width = '75%',
                              height = '600px'),
                   ## Table with the results
                   tags$h4('Results Table'),
                   dataTableOutput(outputId = 'ips_table'),
                   downloadButton(outputId = 'down_ips_table', 
                                  label = 'Download table'),
                   # Immunophenogram 
                   tags$h4('Immunophenogram'),
                   downloadButton(outputId = 'ic_phenogram', 
                                  label = 'Download report')
                 )
                
              ), # End IC_module

  #############################################################################        
  # About and Getting started --------------------------------------------------         
  #############################################################################        
  tabPanel(title = 'About',
           fluidPage(
               htmlOutput(outputId = 'about')
           )
  ),
  
  tabPanel(title = 'Getting started',
           fluidPage(
               htmlOutput(outputId = 'manual')
           )
           
  ), # End Manual section
  
  # Align these two tabs to the right
  tags$head(tags$style('.navbar-nav :nth-child(5) {float:right}
                          .navbar-nav :nth-child(6) {float:right}'))
)

#############################################################################        
#############################################################################        
# SERVER -----------------------------------------------------          
############################################################################# 
#############################################################################        
server <- function(input, output, session) {
  
  #############################################################################        
  # Include the GEGVICshine About and Getting started page in HTML ------------         
  #############################################################################       

  output$about <- renderUI({
      tags$iframe(seamless="seamless",
                  src = 'gegvic_shine_about.html', 
                  width = '100%',  
                  height = 1000,  style = "border:none;")
  })  
  
  output$manual <- renderUI({
    tags$iframe(seamless="seamless",
                src = 'gegvic_shine_manual.html', 
                width = '100%',  
                height = 1000,  style = "border:none;")
  })
  
  
  #############################################################################        
  # Buttons to switch sections in parameters tab ------------------------------        
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
  # Download button Sample Data ------------------------------------------------         
  #############################################################################       
  
  output$down_sample_data <- downloadHandler(
    filename = function(){
      'sample_data.zip'
    },
    content = function(file){
      file.copy('www/sample_data.zip', file)
    },
    contentType = "application/zip"
  )
  
    #############################################################################        
    # Automatic detection of colnames in the metadata---------------------------         
    #############################################################################       
    
    # Create reactives values to use later on
    reactives <- reactiveValues(
        temp_meta = NULL
    )
    
    # response: React to metadata upload ---------------------------------------
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
                          label = 'Select response variable (GE, GV, IC)', 
                          choices  = colnames(reactives$temp_meta))
        
    })
    
    # ref_level --------------------------------------
    observeEvent(input$to_parameters, {
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
      sel.column <- rlang::sym(input$response)
      sel.column <- reactives$temp_meta %>% 
        dplyr::select(sel.column)
      
      # Update the selection input for response
      updateSelectInput(inputId = 'ref_level',
                        label = 'Reference level: Name of the grouping variable (GE)', 
                        choices  = unique(sel.column))
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
        

        # biomart  ------------------------------------------------------------
        biomart <- reactive({
          if(input$biomart == 'ensembl_biomart_GRCh38_p13'){
            GEGVIC::ensembl_biomart_GRCh38_p13
          } else if(input$biomart == 'ensembl_biomart_GRCh37') {
            GEGVIC::ensembl_biomart_GRCh37
          } else if(input$biomart == 'ensembl_biomart_GRCm38_p6') {
            GEGVIC::ensembl_biomart_GRCm38_p6
          } else {
            GEGVIC::ensembl_biomart_GRCm39
          }
        })
        
        # gsva_gmt  -----------------------------------------------------------
        gsva_gmt <- reactive({
          if(input$gsva_gmt == 'Hallmark'){
            'Hallmark'
          } else {
            # Ensure the file is uploaded and available
            req(input$gmt)
            # Read the name in a new variable 
            file <- input$gmt
            # Get the datapath
            ext <- tools::file_ext(file$datapath)
            # Validate that the file is a .gmt
            validate(need(ext == "gmt", "Please upload a gmt file"))
            # export data_path
            file$datapath
          }
        })
        
        # indications  --------------------------------------------------------
        indications <- reactive({
          
          rep(input$indications, nrow(metadata()))
          
        })
        
        # Cibersort  ----------------------------------------------------------
        cibersort_R <- reactive({
          # Read the name in a new variable 
          file <- input$cibersort_R
          # export data_path
          file$datapath
        })
        cibersort_LM22 <- reactive({
          # Read the name in a new variable 
          file <- input$cibersort_LM22
          # export data_path
          file$datapath
        })
        
        # Create a Progress object ---------------------------------------------
        progress <- shiny::Progress$new()
        # Make sure it closes when we exit this reactive, even if there's an error
        on.exit(progress$close())
        
      
        #######################################################################        
        # Generate outputs ---------------------------         
        #######################################################################  
        
        withProgress(message = NULL,
                     detail = NULL,
                     min = 0,
                     max = 10, 
                     expr = {
        
        # GE_module -----------------------------------------------------------  
        
        if(input$ge_module == TRUE) {
          
          # Output: PCA
          output$pca <- renderPlot({
              s_pca(counts = counts(),
                    genes_id = input$genes_id,
                    metadata = metadata(),
                    response = input$response,
                    design = input$design,
                    colors = colors())
          })
          
          
          # Differential gene expression calculation
          setProgress(value = 1, message = 'Calculating differential gene expression')
          
          results.dds <- GEGVIC::ge_diff_exp(counts = counts(),
                                             genes_id = input$genes_id,
                                             metadata = metadata(),
                                             design = input$design,
                                             ref_level = ref_level(),
                                             shrink = input$shrink)
          
          setProgress(value = 2, message = 'Annotating gene symbols')
          
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
                               # 1. Differential gene expression table
                               wellPanel(
                                 tags$h4('Table of Differentially Expressed Genes'),
                                 dataTableOutput(outputId = paste0(n, '_table')),
                                 downloadButton(outputId = paste0('down_',n, '_table'), 
                                                label = 'Download table'),
                                 ),
                               # 2. Volcano plot
                               wellPanel(
                                 tags$h4('Volcano Plot'),
                                 #Plot Volcano
                                 plotOutput(outputId = paste0(n, '_volcano'), 
                                            width = '75%',
                                            height = '600px')
                                 ),
                               # 3. GSEA
                               wellPanel(
                                 tags$h4('GSEA'),
                                 # Add text in case no gene set is found
                                 htmlOutput(outputId = paste0(n, '_no_gsea')),
                                 # Add table
                                 tags$h4('Table of Enriched Gene Sets'),
                                 dataTableOutput(outputId = paste0(n, '_gsea_table')),
                                 downloadButton(outputId = paste0('down_',n, '_gsea_table'), 
                                                label = 'Download table'),
                                 # Plot bubble graph
                                 tags$h4('Bubble Plot'),
                                 plotOutput(outputId = paste0(n, '_bubble_plot'),
                                            width = '75%',
                                            height = '800px'),
                                 # Plot gsea cluster
                                 tags$h4('Gene Set Clusters'),
                                 plotOutput(outputId = paste0(n, '_gsea_clust'), 
                                            width = '75%'#, height = '800px'
                                            ),
                                 # Plot gsea wordclouds
                                 tags$h4('Gene Set Enriched Terms'),
                                 plotOutput(outputId = paste0(n, '_gsea_word'), 
                                            width = '75%',
                                            height = '800px'),
                                 # Plot gsea leading edge
                                 tags$h4('Gene Set Core Enrichment (Leading Edge)'),
                                 plotOutput(outputId = paste0(n, '_gsea_core'),
                                            width = '75%',
                                            height = '800px')
                                 )),
                      select = TRUE
                    )
            
                # Generate Output: Differential gene expression table
                output[[paste0(n, '_table')]] <- DT::renderDataTable(
                  annot.res[[x]] %>% 
                    dplyr::mutate_if(is.numeric, function(x) round(x, 4)),
                  options = list(scrollX = TRUE,
                                 pageLength = 10)
                )
                
                output[[paste0('down_', n, '_table')]] <- downloadHandler(
                  filename = function(){
                    paste0(n, '_diff_exprs_table.csv')
                    },
                  content = function(f){
                    write.csv(x = annot.res[[x]],
                              file = f)
                  }
                )
                
                # Generate Output: Volcano Plot
                output[[paste0(n, '_volcano')]] <- renderPlot({
                  s_volcano(annot_res = annot.res[[x]],
                            fold_change = input$fold_change,
                            p.adj = input$p.adj)
                  
                })
                
                # Generate Output: GSEA
                setProgress(value = 3, message = 'Calculating GSEA')
                ## Perform GSEA
                gsea <- s_gsea(annot_res = annot.res[[x]],
                               gmt = gmt(),
                               gsea_pvalue = input$gsea_pvalue)
                ## In case no p-value is too restrictive
                if (names(gsea)[1] == 'no_genesets') {
                  output[[paste0(n, '_no_gsea')]] <- renderText({
                    
                    '<b> No gene sets are enriched under specific pvalueCutoff!! </b>'
                    
                  })
                  
                 
                } else {
                  setProgress(value = 7, message = 'Plotting GSEA')
                  # Output: GSEA Table
                  output[[paste0(n, '_gsea_table')]] <- DT::renderDataTable(
                    
                    gsea$gsea_result %>% 
                      dplyr::mutate_if(is.numeric, function(x) round(x, 4)) %>% 
                      tibble::rownames_to_column('name') %>% 
                      dplyr::select(-name), 
                    options = list(scrollX = TRUE,
                                   pageLength = 10,
                                   columnDefs = list(list(targets = c(10,11),
                                                          render = JS(
                                                            "function(data, type, row, meta) {",
                                                            "return type === 'display' && data.length > 6 ?",
                                                            "'<span title=\"' + data + '\">' + data.substr(0, 6) + '...</span>' : data;",
                                                            "}")))))
                  # Add download button
                  output[[paste0('down_', n, '_gsea_table')]] <- downloadHandler(
                    filename = function(){
                      paste0(n, '_gsea_table.csv')
                    },
                    content = function(f){
                      write.csv(x = gsea$gsea_result %>% 
                                  dplyr::mutate_if(is.numeric, function(x) round(x, 4)) %>% 
                                  tibble::rownames_to_column('name') %>% 
                                  dplyr::select(-name),
                                file = f)
                    }
                  )
                  # Output: Bubble plot
                  output[[paste0(n, '_bubble_plot')]] <- renderPlot({
                    gsea$gs.filt %>%
                      separate(leading_edge , into= 'tags', sep=',') %>%
                      separate(tags, into = c('tags', 'core_perc'), sep='=') %>%
                      separate(core_perc, into = 'core_perc', sep='%') %>%
                      mutate(core_perc = as.numeric(core_perc)) %>%
                      dplyr::select(-tags) %>%
                      ggplot(.,
                             aes(x= NES,
                                 y=reorder(ID, NES),
                                 size= core_perc,
                                 colour = p.adjust))+
                      geom_point(alpha=0.5)+
                      geom_vline(xintercept = 0)+
                      scale_color_gradient(low = "#FF9900", high = "#FF3300")+
                      labs(size='% of genes in\n leading edge', colour = 'p.adjust')+
                      theme_bw()+
                      theme(panel.grid = element_blank(),
                            axis.text = element_text(size=12, face = "bold"),
                            axis.title.y = element_blank(),
                            axis.title.x = element_text(size=15),
                            axis.text.x = element_text(size=10),
                            legend.title = element_text(face='bold', size =8),
                            legend.text = element_text(size =7))  
                  })
                  # Output: GSEA cluster
                  output[[paste0(n, '_gsea_clust')]] <- renderPlot({
                    GSEAmining::gm_dendplot(df = gsea$gs.filt,
                                            hc = gsea$gs.cl)
                  })
                  # Output: Plot enriched terms in gene sets names
                  output[[paste0(n, '_gsea_word')]] <- renderPlot({
                         
                    GSEAmining::gm_enrichterms(df = gsea$gs.filt,
                                               hc = gsea$gs.cl)
                  })
                  # Output: Plot enriched cores (leading edge analysis)
                  output[[paste0(n, '_gsea_core')]] <- renderPlot({
                    
                    GSEAmining::gm_enrichcores(df = gsea$gs.filt,
                                               hc = gsea$gs.cl)
                  })
                }
                
          })
          
          # Output: GSVA
          gsva.out <- s_single(counts = counts(),
                               metadata = metadata(),
                               genes_id = input$genes_id,
                               response = input$response,
                               design = input$design,
                               biomart = biomart(),
                               gsva_gmt = gsva_gmt(),
                               method = input$gsva_method,
                               colors = colors(),
                               row.names = input$gsva_rownames,
                               col.names = input$gsva_colnames)
        
          ## Table of predicted gene set values
          output[['gsva_table']] <- DT::renderDataTable(
            as.data.frame(gsva.out$gsva_table) %>% 
              dplyr::mutate_if(is.numeric, function(x) round(x, 3)),
            options = list(scrollX = TRUE,
                           pageLength = 10)
          )
          
          output[['down_gsva_table']] <- downloadHandler(
            filename = function(){
              'table_gsva.csv'
            },
            content = function(f){
              write.csv(x =  gsva.out$gsva_table,
                        file = f)
            }
          )
          ## Heatmap
          output$gsva <- renderPlot({
            gsva.out$heatmap
          })
        }
        
        # GV_module -----------------------------------------------------------
          
        if(input$gv_module == TRUE){
        
          # Output: Mutational summary
          setProgress(value = 4, message = 'Summarizing mutations')
          
          output$mut_sum <- renderPlot({
            
            s_mut_summary(muts = muts(),
                          metadata = metadata(),
                          top_genes = input$top_genes)
            
          })
          
          # Output: Oncoplot
          setProgress(value = 5, message = 'Generating Oncoplot')
          
          output$oncoplot <- renderPlot({
          
            s_oncoplot(muts = muts(),
                       metadata = metadata(),
                       response = input$response,
                       top_genes = input$top_genes,
                       col.names = input$mut_colnames,
                       colors = colors())
            
          })
          
          # Mutational load
          setProgress(value = 6, message = 'Calculating mutational load')
          
          ## Output: Mutational Load plot
          mut.load <- s_mut_load(muts = muts(),
                                 metadata = metadata(), 
                                 response = input$response,
                                 compare = input$compare,
                                 p_label = input$p_label,
                                 colors = colors())
                      
          output$mut_load <- renderPlot({
            
            mut.load$mut.load.plot
            
          })
          
          ## Output: Mutational Load Table
          output[['mut_load_table']] <- DT::renderDataTable(
            as.data.frame(mut.load$mut.load.table) %>% 
              dplyr::mutate_if(is.numeric, function(x) round(x, 3)),
            options = list(scrollX = TRUE,
                           pageLength = 10)
          )
          
          output[['down_mut_load_table']] <- downloadHandler(
            filename = function(){
              'table_mut_load.csv'
            },
            content = function(f){
              write.csv(x =  mut.load$mut.load.table,
                        file = f)
            }
          )

          
          # Mutational signatures
          setProgress(value = 7, message = 'Predicting mutational signatures')
          
          ## Predictions
          mut.sigs <- s_mut_signatures(muts = muts(),
                                       metadata = metadata(), 
                                       response = input$response,
                                       gbuild = input$gbuild,
                                       mut_sigs = input$mut_sigs,
                                       colors = colors(),
                                       col.names = input$mut_colnames)
          ## Output: Bar plot
          output$mut_sigs_bar <- renderPlot({
            
            mut.sigs$mut_sig_barplot
            
            
          })
          ## Output: Heatmap
          output$mut_sigs_heat <- renderPlot({
            
            mut.sigs$mut_sig_heatmap
            
          })
          ## Output: Predictions Table
          output[['mut_sig_table']] <- DT::renderDataTable(
            as.data.frame(mut.sigs$mut_sig_table) %>% 
              dplyr::mutate_if(is.numeric, function(x) round(x, 3)),
            options = list(scrollX = TRUE,
                           pageLength = 10)
          )
          
          output[['down_mut_sig_table']] <- downloadHandler(
            filename = function(){
              'table_mut_sig.csv'
            },
            content = function(f){
              write.csv(x =  mut.sigs$mut_sig_table,
                        file = f)
            }
          )
        }
        
        # IC_module -----------------------------------------------------------  
        
        if(input$ic_module == TRUE){
          
          # TPM calculation
          setProgress(value = 8, message = 'Calculating TPM')
          
          tpm <- GEGVIC::ic_raw_to_tpm(counts = counts(),
                                       genes_id = input$genes_id,
                                       biomart = biomart())
          

          # Immune prediction
          setProgress(value = 9, message = 'Predicting immune cell populations')
          
          ic.pred <- s_deconv(gene_expression = tpm,
                                      indications = indications(),
                                      cibersort_R = cibersort_R(),
                                      cibersort_LM22 = cibersort_LM22(),
                                      tumor = TRUE,
                                      rmgenes = NULL,
                                      scale_mrna = TRUE,
                                      expected_cell_types = NULL)
          
          # Generate Output: Table of predicted immune cell populations
          output[['ic_pred_table']] <- DT::renderDataTable(
            ic.pred %>% 
              dplyr::mutate_if(is.numeric, function(x) round(x, 3)),
            options = list(scrollX = TRUE,
                           pageLength = 10)
          )
          
          output[['down_ic_pred_table']] <- downloadHandler(
            filename = function(){
              'table_immune_prediction.csv'
            },
            content = function(f){
              write.csv(x = ic.pred,
                        file = f)
            }
          )
          
          
          ## Output: Plot immune prediction comparison of cell types by groups
          output$ic_samples <- renderPlot({
            
            s_plot_comp_samples(df = ic.pred,
                                metadata = metadata(),
                                response = input$response,
                                compare = input$compare,
                                p_label = input$p_label,
                                colors = colors(),
                                points = input$ic_samples_points)
            
          })
          
          ## Output: Plot immune prediction comparison within samples
          output$ic_type <- renderPlot({
            
            s_plot_comp_celltypes(df = ic.pred,
                                  metadata = metadata(),
                                  response = input$response,
                                  col.names = input$mut_colnames)
            
          })
          
          
          setProgress(value = 10, message = 'Calculating immune score')
          
          # Calculate IPS and IPG
          ic.IPS_IPG <- s_score(tpm = tpm,
                                metadata = metadata(),
                                response = input$response,
                                compare = input$compare,
                                p_label = input$p_label,
                                colors = colors())
          
          ## Output: Plot IPG
          #output$ic_phenogram <- renderPlot({
          #  ic.IPS_IPG$immunophenoGram
          #})
          
          ## Output: Plot IPS
          output$ic_score<- renderPlot({
            ic.IPS_IPG$immunophenoScore
          })
          
         
          ## Output: Predictions Table
          output[['ips_table']] <- DT::renderDataTable(
            as.data.frame(ic.IPS_IPG$ips_table) %>% 
              dplyr::mutate_if(is.numeric, function(x) round(x, 3)),
            options = list(scrollX = TRUE,
                           pageLength = 10)
          )
          
          output[['down_ips_table']] <- downloadHandler(
            filename = function(){
              'table_ips.csv'
            },
            content = function(f){
              write.csv(x =  ic.IPS_IPG$ips_table,
                        file = f)
            }
          )
          
          ## Output: Download IPG report in pdf
          output[['ic_phenogram']] <- downloadHandler(
            filename = function(){
              'immunophenogram_report.pdf'
            },
            content = function(f){
              ggsave(filename = f, 
                     plot = ic.IPS_IPG$immunophenoGram, 
                     device = 'pdf')
            }
          )
        }
        
      }) # End progress bar
      
    }) # End button
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
