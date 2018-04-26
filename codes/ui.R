ui = fluidPage(navbarPage(title = "Synergistic and Antagonistic Interaction Learner",
  
  #                         tags$head(
  #                           tags$style(HTML("
  #                                           .shiny-options-group { 
  # height: auto;
  #                                           width: 920px;
  #                                           -webkit-column-count: 6; /* Chrome, Safari, Opera */ 
  #                                           -moz-column-count: 6;    /* Firefox */ 
  #                                           column-count: 1; 
  #                                           -webkit-column-fill: balance;
  #                                           -moz-column-fill: balance;
  #                                           column-fill: balance;
  #                                           margin-top: 2px;
  #                                           } 
  #                                           
  #                                           .control-label {
  #                                           padding-bottom: 2px;
  #                                           }
  #                                           
  #                                           div.radio {
  #                                           margin-top: 2px;
  #                                           margin-bottom: 2px;
  #                                           padding-bottom: 2px;
  #                                           }
  #                                           .option-header {
  #                                           color: #79d;
  #                                           text-transform: uppercase;
  #                                           margin-bottom: 2px;
  #                                           }
  #                                           "))
  #                           ),
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
             h4("Combinatorial treatments as easy as 1-2-3"),
             br(),
             br(),
             h5("Watch Demo"),
             embed_youtube("mIwhl8g5uKU", allowfullscreen = TRUE, width = 300, height = 250)
             
             #actionButton("go_to_demo", "Watch Demo", icon("play"))
    ),
    
    tabPanel("1. Upload Data",
             
             sidebarLayout(
               
               sidebarPanel(width = 3,
                            
                            #h5("Import Data"),
                            
                            radioButtons("upload",
                                         "Import Data:",
                                         c("from file" = 'go_to_upload',
                                           "from Database" = 'go_to_DB')
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

                            
                            actionButton("start2", "Run analysis", icon("play"))
                            
                            ),
               
               mainPanel(
                 
                 h5(textOutput("x_samples_1")),
                 
                 h5(textOutput("y_samples_1")),
                 

                 plotlyOutput("pca_plot"))
               
               )
             ),
    
    tabPanel("3. Browse Results",
             
             tags$style(type="text/css",
                        ".shiny-output-error { visibility: hidden; }",
                        ".shiny-output-error:before { visibility: hidden; }"
             ),
             
             h5(textOutput("x_samples_2")),
             h5(textOutput("y_samples_2")),
             uiOutput("imageGrid"),
             
             #uiOutput("imageGrid1"),
             
             br(),
             
             column(12,
                    
                    fluidRow(
                      
                      column(3, br(),
                             
                             #textOutput("click_case"),
                             #textOutput("click_case2"),
                             
                             selectInput("case", "Select case:", paste(1:17)),
                             
                             imageOutput('case_img'),
                             textOutput("tt")
                            
                             #imageOutput('case_img_new')
                             
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
                             
                             br(),
                             
                             dataTableOutput("interactions_table"),
                             
                             downloadButton("downloadData", "Download data")
                             
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
    
    tabPanel('Combinatorial Interaction Database',
             br(),
             searchInput(
               inputId = "id", 
               label = "Enter a Gene Symbol:", 
               placeholder = "e.g. IL17", 
               btnSearch = icon("search"), 
               btnReset = icon("remove"), 
               width = "30%"
             ))
             
    
    
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