library(plotly)
library(shinyWidgets)
library(reshape2)

source("chooser.R")

#setPowerPointStyle()

options(warn=-1) 

load("IFN_TNF_1h")
immune_DB = read.csv('immune_selection/DB_immune_selection.txt', sep = '\t')
immune_DB = immune_DB[,1:12]
immune_DB = immune_DB[,-9]

ui = fluidPage(
  navbarPage("DrugSynergy: Resolving Combinatorial Treatments",
             
             tabsetPanel(id = "tabs",
                         
                         tabPanel("Import Data",
                                  
                                  sidebarLayout(
                                    
                                    sidebarPanel(width = 4,
                                                 
                                                 h4("upload data"),
                                                 fileInput("file", label = ""),
                                                 h4("set columns"),
                                                 dropdownButton(
                                                   uiOutput("dynamic"),
                                                   #switchInput(inputId = "id", value = FALSE),
                                                   circle = TRUE, status = "danger", icon = icon("gear"), width = "300px",
                                                   tooltip = tooltipOptions(title = "Define Samples")
                                                 ),
                                                 br(),
                                                 actionButton("start", "Quality control", icon("play")),
                                                 br(),
                                                 br(),
                                                 h4("or choose from Database"), 
                                                 br(),
                                                 actionButton("go_to_DB", "Go to Immune X+Y DB"),
                                                 #selectInput("geo", label = "", choices = make.names(names(immune_DB$GEO.ID), unique = T)),
                                                 textOutput("result")
                                                 
                                    ),
                                    
                                    mainPanel(
                                      
                                      textOutput('summary_text'),
                                      
                                      br(),
                                      
                                      dataTableOutput("input_data_table"),
                                      
                                      br(),
                                      textOutput("test")
                                      
                                      
                                    )
                                  )
                         ),
                         
                         tabPanel("Quality Control",
                                  
                                  sidebarLayout(
                                    
                                    sidebarPanel(width = 4,
                                                 
                                                 actionButton("start1", "Go to analysis", icon("share"))
                                    ),
                                    
                                    mainPanel(
                                      
                                      plotlyOutput("pca_plot")
                                      
                                    )
                                    )
                                  ),         
                         
                         tabPanel("Run analysis",
                                  
                                  sidebarLayout(
                                    
                                    sidebarPanel(width = 4,
                                                 
                                                 
                                                 radioButtons("platform", "Data type:",
                                                              c("microarray" = 'ma',
                                                                "RNA-seq (counts)" = 'rnaseq')),
                                                 
                                                 sliderInput("pval", "Significance threshold:", 
                                                              min = 0.001, max = 0.1, value = 0.05),
                                                 
                                                 sliderInput("minfc", "Minimum fold-change (log2):", 
                                                             min = 0.1, max = 1, value = 0.5),
                                                    
                                                 actionButton("start2", "Run analysis", icon("play"))
                                                 
                                    ),
                                    mainPanel(
                                      
                                    )
                                  )
                         ),    
                         
                         
                         tabPanel("Browse interactions",
                                  
                                  #sidebarLayout(
                                    
                                   # sidebarPanel(width = 4,
                                                 
                                    #             plotlyOutput("sankey")
                                    #),
                                    
                                    #mainPanel(
                                  br(),
                                  uiOutput("imageGrid"),
                                  tags$script(HTML(
                                    "$(document).on('click', '.clickimg', function() {",
                                    "  Shiny.onInputChange('clickimg', $(this).data('value'));",
                                    "});"
                                  )),
                                  
                                  
                                  br(),
                                    column(12, 
                                           fluidRow(
                                             
                                             column(3, br(), 
                                                    selectInput("case", "Select case:", paste(1:17)),
                                                    imageOutput('case_img')),
                                             
                                             column(1, br()),
                                             
                                             column(4, br(), 
                                                    radioButtons("interaction", "Top interactions:", 
                                                                 c("synergistic" = "P", "antagonistic" = "N"), selected = "P", inline = T), 
                                                    actionButton("explore1", 'functional enrichment', icon("pie-chart")),
                                                    #downloadButton("downloadData", "Download data"),
                                                    br(),
                                                    dataTableOutput("interactions_table")),
                                             column(1, br()),
                                             column(4, #br(), 
                                                    #selectInput("gene", "Select gene:", textOutput("genes")),
                                                    uiOutput("choose_gene"),
                                                    #br(),
                                                    plotlyOutput('gene_plot', width = "100%", height = "280px"))
                                           )

                                  )
                                  ),
                         
                         tabPanel(
                           'Functions/Pathways',
                           dataTableOutput('enrich_tab')
                         ),
                         
                         tabPanel(
                           'Immune X + Y Database',
                           dataTableOutput('immune_DB')
                         ),
                         
                         tabPanel(
                           'Help',
                           tags$video(src = 'test_app.mov')
                           
                           #uiOutput("imageGrid")
                           #tags$script(HTML(
                          #   "$(document).on('click', '.clickimg', function() {",
                          #   "  Shiny.onInputChange('clickimg', $(this).data('value'));",
                          #   "};"
                          # ))
                         
                         )
                         
                         
                         
                         )
             )
)

