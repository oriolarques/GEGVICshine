library(shiny)
library(GEGVIC)
library(tidyverse)

# Define UI 
ui <- fluidPage(
    
    tabsetPanel(
        
        tabPanel(title = 'Parameters',
                 ## Inputs ----------------------------------------------------
                 # Input: RNA-seq raw counts
                 fileInput(inputId = "counts", 
                           label = 'RNA-seq raw counts',
                           accept = c('.csv')),
                 # Input: Metadata
                 fileInput(inputId = "metadata", 
                           label = 'Metadata',
                           accept = c('.csv')),
                 # Input: Genetic variations
                 fileInput(inputId = "muts", 
                           label = 'Genetic Variations',
                           accept = c('.csv')),
                 ## Parameters ------------------------------------------------
                 # genes_id
                 selectInput(inputId = 'genes_id', 
                             label = 'Genes ID', 
                             choices = list('Genes ID' = c("hgnc_symbol",
                                                           "entrezgene_id", 
                                                           "ensembl_gene_id"))),
                 # design
                 textInput(inputId = 'design',
                           label = 'Design'),                
                 # colors
                 textInput(inputId = 'colors',
                           label = 'Colors: Enter a vector (comma separated)'), 
                 # ref_level
                 textInput(inputId = 'ref_level',
                           label = 'Reference level: Name of the grouping variable, 
                           Name of the group to be considered the reference level 
                           (comma separated)'), 
                 # shrink
                 selectInput(inputId = 'shrink', 
                             label = 'Shrinkage method', 
                             choices = list('Shrinkage' = c("apeglm", "ashr", 
                                                            "normal", "none")),
                             selected = 'apeglm'),
                 # Input: biomart
                 fileInput(inputId = "biomart", 
                           label = 'BiomaRt database',
                           accept = c('.csv')),
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
                 # Input: gmt
                 fileInput(inputId = "gmt", 
                           label = 'Gene sets (as .gmt file)',
                           accept = c('.gmt')),
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
                 # cibersort
                 
                

            #ic
                 #cibersort = NULL,
                 #response,
                 #compare = NULL,
                 #p_label = 'p.format',
            ##gv
                 #muts, 
                 #response,
                 #top_genes = 10,
                 #specific_genes = NULL,
                 #compare = NULL,
                 #p_label = 'p.format',
                 #gbuild = 'BSgenome.Hsapiens.UCSC.hg19',
                 #mut_sigs = COSMIC_v2_SBS_GRCh37,
                 #tri.counts.method = 'default'
            
            actionButton(inputId = "Execute", label = 'Execute')
                 ),
        tabPanel(title = 'GE_module',
           
            
            #Plot PCA
            plotOutput(outputId = 'pca')
            
        )
        
        
        
        
    )
    
    
    
    
)

# Define server logic 
server <- function(input, output) {
    
    observeEvent(input$Execute, {
        # Read the counts input
        counts <- reactive({
            # Ensure the file is uploaded and available
            req(input$counts)
            # Read the name in a new variable 
            file <- input$counts
            # Get the datapath
            ext <- tools::file_ext(file$datapath)
            # Validate that the file is a .csv
            validate(need(ext == "csv", "Please upload a csv file"))
            # 
            temp_df <- read.csv(file = file$datapath, 
                                header = TRUE, sep = ',')
            
            if (ncol(temp_df) != 1) {
                read.csv(file = file$datapath, 
                         header = TRUE, sep = ',')
            } else {
                read.csv(file = file$datapath, 
                         header = TRUE, sep = ';')
            }
            
        })
        
        # Read the metadata input
        metadata <- reactive({
            # Ensure the file is uploaded and available
            req(input$metadata)
            # Read the name in a new variable 
            file <- input$metadata
            # Get the datapath
            ext <- tools::file_ext(file$datapath)
            # Validate that the file is a .csv
            validate(need(ext == "csv", "Please upload a csv file"))
            # 
            temp_df <- read.csv(file = file$datapath, 
                                header = TRUE, sep = ',')
            
            if (ncol(temp_df) != 1) {
                read.csv(file = file$datapath, 
                         header = TRUE, sep = ',')
            } else {
                read.csv(file = file$datapath, 
                         header = TRUE, sep = ';')
            }
            
        })
        
        # colors input 
        colors <- reactive({
            
            unlist(strsplit(input$colors,split = ','))
            
        })
        
        output$pca <- renderPlot({
            GEGVIC::ge_pca(counts = counts(),
                           genes_id = input$genes_id,
                           metadata = metadata(),
                           design = input$design,
                           colors = colors())
            
        })
        
        
        
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)
