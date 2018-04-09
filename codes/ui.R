ui = fluidPage(navbarPage(title = "Synergistic and Antagonistic Interaction Learner",
  
  tabsetPanel(
    
    id = "tabs",
    
    tabPanel("Welcome",
             
             #actionButton("go_to_GEO", "From Gene Omnibus", style='height:300px'),
             #actionBttn("go_to_GEO", "From Gene Omnibus", style='jelly', no_outline = FALSE, size = "lg", icon = 'case_1.png'),
             #tags$button(
            #  id = "web_button",
            #   class = "btn action_button",
            #   img(src = "case_1.png",
            #       height = "200px"),
            #       tags$style(HTML('color: #4d3a7d;'))
            # ),
             br(),
             h4("Here description and logo"),
             br(),
             br(),
             br(),
             br(),
             br(),
             h4("Demo"),
             embed_youtube("mIwhl8g5uKU", allowfullscreen = TRUE, width = 300, height = 250)
             
             #actionButton("go_to_demo", "Watch Demo", icon("play"))
    ),
    
    tabPanel("1. Upload Data",
             
             sidebarLayout(
               
               sidebarPanel(width = 4,
                            
                            #h5("Import Data"),
                            
                            radioButtons("upload",
                                         "Import Data:",
                                         c("from file" = 'go_to_upload',
                                           "from Database" = 'go_to_DB',
                                           "from Gene Omnibus" = 'go_to_GEO')
                            ),
                            
                            # actionButton("go_to_upload", "From File", style='width:200px'),
                            # 
                            # br(),
                            # 
                            # br(),
                            # 
                            # actionButton("go_to_DB", "From Immune X+Y Database", style='width:200px'),
                            # 
                            # br(),
                            # 
                            # br(),
                            # 
                            # actionButton("go_to_GEO", "From Gene Omnibus", style='width:200px'),
                            
                            br(),
                            
                            br(),
                            
                            br()
                            
                            ),
               
               mainPanel(
                 
                 uiOutput("main")
                 
                 #fileInput("file", label = ""),
                 
                 #textOutput('summary_text'),
                 
                 #dataTableOutput("input_data_table"),
                 
                 #actionButton("start", "Quality control", icon("play"))
                 
                 #h2("Here Logo + Text/images")
                 
                 )
               
               )
             
             ),
    
               
                      #br(),
               #h4("Identify synergistic and antagonistic interactions in 
              #         high-throughput combinatorial treatments"),
               #div(style="text-align:justify","Identify
              #     synergistic and antagonistic interactions in
              #     high-throughput combinatorial treatments")
            
    # tabPanel("Upload Data",
    #          
    #          sidebarLayout(
    #            
    #            sidebarPanel(width = 4,
    #                         
    #                         h4("Upload file"),
    #                         
    #                         #fileInput("file", label = ""),
    #                         
    #                         h4("set columns"),
    #                         
    #                         dropdownButton(
    #                           
    #                           uiOutput("set_columns"),
    #                           
    #                           circle = TRUE,
    #                           
    #                           status = "danger",
    #                           
    #                           icon = icon("gear"),
    #                           
    #                           width = "300px",
    #                           
    #                           tooltip = tooltipOptions(title = "Define Samples")
    #                           
    #                           ),
    #                         
    #                         br(),
    #                         
    #                         #actionButton("start", "Quality control", icon("play")),
    #                         
    #                         textOutput("result")
    #                         
    #                         ),
    #            
    #            mainPanel(
    #              
    #              #textOutput('summary_text'),
    #              
    #              #br(),
    #              
    #              #dataTableOutput("input_data_table"),
    #              
    #              br(),
    #              
    #              textOutput("test")
    #              
    #              )
    #            
    #            )
    #          
    #          ),
    
    tabPanel("2. Check Quality and Run",
             
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
    
    tabPanel("3. Browse Results",
             
             tags$style(type="text/css",
                        ".shiny-output-error { visibility: hidden; }",
                        ".shiny-output-error:before { visibility: hidden; }"
             ),
             
             br(),
             
             uiOutput("imageGrid"),
             
             #uiOutput("imageGrid1"),
             
             br(),
             
             column(12,
                    
                    fluidRow(
                      
                      column(3, br(),
                             
                             textOutput("click_case"),
                             #textOutput("click_case2"),
                             
                             selectInput("case", "Select case:", paste(1:17)),
                             
                             imageOutput('case_img'),
                             
                             imageOutput('case_img_new')
                             
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
             
             )
    
    # tabPanel('Immune X + Y Database',
    #          
    #          textOutput('myText')
    #          
    #          #dataTableOutput('immune_DB')
    #          
    #          ),
    
    # tabPanel('Demo',
    #          
    #          embed_youtube("mIwhl8g5uKU", allowfullscreen = TRUE, width = 500)
    #          
    #          )
    )
  
  )
  )