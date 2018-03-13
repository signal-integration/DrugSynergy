library(plotly)
library(shinyWidgets)

source("setPowerPointStyle.R")
source("chooser.R")

setPowerPointStyle()


load("GEOmatrix (1)")
load("GSE43452_filtered_matched")

my_results_2.5h = my_data_filtered_matched
genes = sort(as.character(my_results_2.5h$gene))

ui = fluidPage(
  
  navbarPage("DrugSynergy: Resolving Combinatorial Treatments",
             
             tabsetPanel(id = "tabs",
                         
                         tabPanel("Import Data",
                                  
                                  sidebarLayout(
                                    
                                    sidebarPanel(width = 4,
                                                 
                                                 h4("upload file"),
                                                 fileInput("file", label = ""),
                                                 h4("or enter GEO code"), 
                                                 selectInput("geo", label = "", choices = make.names(names(GEOmatrix[,1]), unique = T)),
                                                 h4("set columns"),
                                                 dropdownButton(
                                                   uiOutput("dynamic"),
                                                   #switchInput(inputId = "id", value = FALSE),
                                                   circle = TRUE, status = "danger", icon = icon("gear"), width = "300px",
                                                   tooltip = tooltipOptions(title = "Define Samples")
                                                 ),
                                                 
                                                 br(""),
                                                 
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
                                      
                                      #plotlyOutput("sankey")
                                      
                                    )
                                  )
                         ),    
                         
                         
                         tabPanel("Browse interactions",
                                  
                                  #sidebarLayout(
                                    
                                   # sidebarPanel(width = 4,
                                                 
                                    #             plotlyOutput("sankey")
                                    #),
                                    
                                    #mainPanel(
                                  
                                    column(7
                                           , fluidRow(
                                             column(6, plotlyOutput("sankey"))
                                             , column(5, br(), selectInput("case", "Choose case:", paste(1:17)), imageOutput('case_img'))
                                           )

                                    ),
                                  column(5
                                         , fluidRow(
                                           column(6, br(), selectInput("gene", "Top synergies:", genes, selected = 'FOXO3'))
                                           , column(6, style = "background-color:green", div(style = "height:400px;"))
                                         )
                                         
                                  )
                                  
                                  
                                    #, column(6, style = "background-color:blue;", div(style = "height:400px;"))

                              

                                      # fluidRow(
                                      #   column(4, br(), plotlyOutput("sankey")
                                      #          ), 
                                      #   
                                      #   column(4, br(), selectInput("case", "Choose case:", paste(1:17)), 
                                      #          imageOutput('case_img')
                                      #          )),
                                      # 
                                      # 
                                      #   column(4, style = "background-color:red;",
                                      #          fluidRow(column(12, style = "background-color:red;")
                                      #                   ) 
                                      #                 
                                      #   
                                      #          
                                      #          )
                                  #),
                                        
                                        
                                        
                                        #)
                                  
                                       
                                      #   column(4, 
                                      #          
                                      #          fluidRow(
                                      #            
                                      #            column(4, style = "background-color:red;")#,
                                      #            #column(4, style = "background-color:blue;")
                                      #            
                                      #            
                                      #          )
                                      # )
                         
                                  
                                  
                                    #)
                                      #selectInput("case", "Choose case:", paste(1:17)),
                                      #fluidRow(column(4, offset = 0, HTML("<div style='height: 210px;'>"), 
                                      #                imageOutput('case_img'), HTML("</div>")))
                                      
                                    #)
             
                                  
                         ),
                         
                         # 
                         # 
                         #                         
                         # 
                         # tabPanel("Browse Responses",
                         #          
                         #          sidebarLayout(
                         #            
                         #            sidebarPanel(width = 4,
                         #                         
                         #                         h4("Cells = U87"),
                         #                         h4("X = Temozolomide"),
                         #                         h4("Y = Y15"),
                         #                         
                         #                         selectInput("case", "Choose case:", paste(1:17)),
                         #                         fluidRow(
                         #                           column(4, offset = 0, HTML("<div style='height: 210px;'>"), imageOutput('case_img'), HTML("</div>")) 
                         #                         ),
                         #                         br(),
                         #                         br()
                         #            ),
                         #            mainPanel(
                         #              
                         #              fluidRow(
                         #                column(4, offset = 0, actionButton("explore1", 'explore gene set', icon("pie-chart")), HTML("<div style='height: 200px;'>"), imageOutput('prof_1_img'), HTML("</div>")), 
                         #                column(4, offset = 0, actionButton("explore2", 'explore gene set', icon("pie-chart")), HTML("<div style='height: 200px;'>"), imageOutput('prof_2_img'), HTML("</div>")),
                         #                column(4, offset = 0, actionButton("explore3", 'explore gene set', icon("pie-chart")), HTML("<div style='height: 200px;'>"), imageOutput('prof_3_img'), HTML("</div>"))
                         #              ),
                         #              
                         #              br(),
                         #              
                         #              fluidRow(
                         #                column(4, offset = 0, actionButton("explore4", 'explore gene set', icon("pie-chart")), HTML("<div style='height: 200px;'>"), imageOutput('prof_4_img'), HTML("</div>")), 
                         #                column(4, offset = 0, actionButton("explore5", 'explore gene set', icon("pie-chart")), HTML("<div style='height: 200px;'>"), imageOutput('prof_5_img'), HTML("</div>")),
                         #                column(4, offset = 0, actionButton("explore6", 'explore gene set', icon("pie-chart")), HTML("<div style='height: 200px;'>"), imageOutput('prof_6_img'), HTML("</div>"))
                         #              ),     
                         #              
                         #              br(),
                         #              
                         #              fluidRow(
                         #                column(4, offset = 0, actionButton("explore7", 'explore gene set', icon("pie-chart")), HTML("<div style='height: 200px;'>"), imageOutput('prof_7_img'), HTML("</div>")), 
                         #                column(4, offset = 0, actionButton("explore8", 'explore gene set', icon("pie-chart")), HTML("<div style='height: 200px;'>"), imageOutput('prof_8_img'), HTML("</div>")),
                         #                column(4, offset = 0, actionButton("explore9", 'explore gene set', icon("pie-chart")), HTML("<div style='height: 200px;'>"), imageOutput('prof_9_img'), HTML("</div>"))
                         #              ),
                         #              
                         #              dataTableOutput("tab")
                         #              
                         #            )
                         #          )
                         # ),
                         # 
                         # # tabPanel(
                         # #   
                         # #   title = "Gene Sets",
                         # #   
                         # #   sidebarLayout(
                         # #     
                         # #     sidebarPanel(
                         # #       width = 3,
                         # #       h4(""),
                         # #       
                         # #       actionButton("enrich", "Enrichment Analysis", icon("play"))
                         # #       #radioButtons("type", " ", choices = c("genes", "enrichment"), selected = "genes")
                         # #       
                         # #     ),
                         # #     
                         # #     mainPanel(
                         # #       
                         # #       br(),
                         # #       
                         # #       br(),
                         # #       
                         # #       fluidRow(
                         # #         
                         # #         downloadButton("downloadData", "Download data"),
                         # #         
                         # #         br(),
                         # #         br(),
                         # #         
                         # #         dataTableOutput("mytab")
                         # #         
                         # #       )
                         # #       
                         # #     )
                         # #   )
                         # #   
                         # # ),
                         # 
                         # tabPanel('Function/Pathways',
                         #          h4("Cells = U87"),
                         #          h4("X = Temozolomide"),
                         #          h4("Y = Y15"),
                         #          
                         #          dataTableOutput('enrich_tab')
                         #          ),
                         # 
                         # tabPanel(
                         #   
                         #   'Look by gene',
                         #   
                         #   sidebarPanel(width = 3,
                         #                h4(""),
                         #                selectInput("gene", "Select gene:",
                         #                            genes, selected = 'FOXO3')),
                         #   
                         #   mainPanel(
                         #     h4("Cells = U87"),
                         #     h4("X = Temozolomide"),
                         #     h4("Y = Y15"),
                         #     
                         #     plotOutput('gene_plot')
                         #     
                         #   )
                         #   ),
                         
                         tabPanel(
                           'X + Y Database',
                           dataTableOutput('francois')
                         )
                         
                         )
             )
)


