ui = fluidPage(navbarPage(title = "Synergistic and Antagonistic Interaction Learner",
  
  tabsetPanel(
    
    id = "tabs",
    
    tabPanel("Welcome",
             
             sidebarLayout(
               
               sidebarPanel(width = 4,
                            
                            h5("Import Data"),
                            
                            actionButton("go_to_upload", "From File", style='width:200px'),
                            
                            br(),
                            
                            br(),
                            
                            actionButton("go_to_DB", "From Immune X+Y Database", style='width:200px'),
                            
                            br(),
                            
                            br(),
                            
                            actionButton("go_to_GEO", "From Gene Omnibus", style='width:200px'),
                            
                            br(),
                            
                            br(),
                            
                            br(),
                            
                            h4("Need help?"),
                            
                            actionButton("go_to_demo", "Watch Demo", icon("play"))
                            
                            ),
               
               mainPanel(
                 
                 h2("Here Logo + Text/images")
                 
                 )
               
               )
             
             ),
    
               
                      #br(),
               #h4("Identify synergistic and antagonistic interactions in 
              #         high-throughput combinatorial treatments"),
               #div(style="text-align:justify","Identify
              #     synergistic and antagonistic interactions in
              #     high-throughput combinatorial treatments")
            
    tabPanel("Upload Data",
             
             sidebarLayout(
               
               sidebarPanel(width = 4,
                            
                            h4("Upload file"),
                            
                            fileInput("file", label = ""),
                            
                            h4("set columns"),
                            
                            dropdownButton(
                              
                              uiOutput("set_columns"),
                              
                              circle = TRUE,
                              
                              status = "danger",
                              
                              icon = icon("gear"),
                              
                              width = "300px",
                              
                              tooltip = tooltipOptions(title = "Define Samples")
                              
                              ),
                            
                            br(),
                            
                            actionButton("start", "Quality control", icon("play")),
                            
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
    
    tabPanel("Check quality and run",
             
             sidebarLayout(
               
               sidebarPanel(width = 4,
                            
                            radioButtons("platform",
                                         "Data type:",
                                         c("microarray" = 'ma',
                                           "RNA-seq (counts)" = 'rnaseq')
                                         ),
                            
                            sliderInput("pval",
                                        "Significance threshold:",
                                        min = 0.001,
                                        max = 0.1,
                                        value = 0.05),
                            
                            sliderInput("minfc",
                                        "Minimum fold-change (log2):",
                                        min = 0.1,
                                        max = 1,
                                        value = 0.5),
                            
                            actionButton("start2", "Run analysis", icon("play"))
                            
                            ),
               
               mainPanel(plotlyOutput("pca_plot"))
               
               )
             ),
    
    tabPanel("Browse results",
             
             tags$style(type="text/css",
                        ".shiny-output-error { visibility: hidden; }",
                        ".shiny-output-error:before { visibility: hidden; }"
             ),
             
             br(),
             
             uiOutput("imageGrid"),
             
             uiOutput("imageGrid1"),
             
             br(),
             
             column(12,
                    
                    fluidRow(
                      
                      column(3, br(),
                             
                             selectInput("case", "Select case:", paste(1:17)),
                             
                             imageOutput('case_img')
                             
                             ),
                      
                      column(1, br()),
                      
                      column(4, br(),
                             
                             radioButtons("interaction",
                                          
                                          "Top interactions:",
                                          
                                          c("synergistic" = "P", "antagonistic" = "N"),
                                          
                                          selected = "P",
                                          
                                          inline = T
                                          
                                          ),
                             
                             actionButton("explore1", 'functional enrichment', icon("pie-chart")),
            
                             #downloadButton("downloadData", "Download data"),
                             
                             br(),
                             
                             dataTableOutput("interactions_table")
                             
                             ),
                      
                      column(1, br()),
                      
                      column(4,
                             
                             uiOutput("choose_gene"),
                             
                             plotlyOutput('gene_plot', width = "100%", height = "280px"),
                             
                             textOutput("t1")
                             
                             )
                      
                      ),
                    
                    dataTableOutput('enrich_tab')
                    
                    )
             
             ),
    
    tabPanel('Immune X + Y Database',
             
             dataTableOutput('immune_DB')
             
             ),
    
    tabPanel('Demo',
             
             embed_youtube("mIwhl8g5uKU", allowfullscreen = TRUE, width = 500)
             
             )
    )
  
  )
  )