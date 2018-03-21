library(plotly)
library(enrichR)
library(DT)

source("filter_data_on_deltas.R")
source("find_optimal_match.R")
source("visualize_all_profiles.R")
#source("setPowerPointStyle.R")
load("IFN_TNF_1h")
load("GEOmatrix (1)")

#read profile definitions
codes = read.table('profile_codes_v2.txt', header = TRUE)

#will be immune-specific DB
GEOmatrix[,1] = make.names(GEOmatrix[,1], unique = T)
GEOmatrix = GEOmatrix[,-2]
row.names(GEOmatrix) = NULL

#genes to be visualized
genes = sort(as.character(degs$genes))


#DB for enrichment analysis
dbs <- c("GO_Biological_Process_2017b", "NCI-Nature_2016",
         "WikiPathways_2016", "Reactome_2016")


server = function(input, output, session) {
  
  options(shiny.maxRequestSize = 30 * 1024 ^ 2)
  
  input_data = reactiveValues()
  
  output$dynamic <- renderUI({
    
    inFile <- input$file
    
    if (is.null(inFile)) return(NULL)
    
    file_columns = dim(input_data$M)[2]
    
    file_columns_names = names(input_data$M)

    LL <- vector("list", file_columns)        
    
    for(i in 1:file_columns){
      
      LL[[i]] <- list(radioButtons(inputId = paste0("mVar",i), 
                                   label = file_columns_names[i], 
                                   choices = c("Info", "CTRL", "X", "Y", "X+Y"),
                                   inline = T),
                      br())
    }      
    
    return(LL)                      
    
  })
  

  
  output$summary_text = renderText({
    
    inFile <- input$file
    
    if (is.null(inFile)) return(NULL)
    
    paste("Uploaded file has", dim(input_data$M)[1], "rows", "and",
          dim(input_data$M)[2], "columns", "(NA=", sum(is.na(input_data)), ")")
    
    })
  
  
  output$test = renderText({
    
    inFile <- input$file
    
    if (is.null(inFile)) return(NULL)
    
    file_columns = dim(input_data$M)[2]
    
    unlist(reactiveValuesToList(input)[paste0("mVar",1:file_columns)])

  })
  
  
  output$define = renderText({
    
    inFile <- input$file
    
    if (is.null(inFile)) return(NULL)
    
    paste("Define samples")
    
  })
  
  
  filtered_data = reactive({
    
    expression_data = as.data.frame.matrix(input_data$M[,-(1:2)])
    
    if (max(expression_data)>25) expression_data = log2(expression_data)
    
    cof = apply(expression_data, 1, function(x) sd(x)/abs(mean(x)))
    
    cof_cutoff = 0.05
    
    cof_filter = which(cof > cof_cutoff)
    
    expression_data[cof_filter, ]

  })
  
  observeEvent(input$start, {
    
    my.pca <- prcomp(t(filtered_data()), center = TRUE, scale=TRUE)
    
    #we assume the same number of samples for each condition
    
    samples = ncol(filtered_data())/4

    cols = c(rep("ctrl", samples), rep("X", samples),
             rep("Y", samples), rep("X+Y", samples))
    
    pca_df = data.frame(PC1 = my.pca$x[, 1], PC2 = my.pca$x[, 2])
    
    row.names(pca_df) = names(input_data$M)[-(1:2)]
    
    pca_df$cols = cols
    
    output$pca_plot = renderPlotly(plot_ly(data = pca_df, x = ~PC1, y = ~PC2, 
            color = ~cols, marker = list(size = 15)) %>% layout(title = 'Principal Component Analysis'))
    
    updateTabsetPanel(session = session, inputId = "tabs", selected = "Quality Control")
     
    })
  
  observeEvent(input$start1, {
    
    updateTabsetPanel(session = session, inputId = "tabs", selected = "Run analysis")
    
  }
  )
  
  
  output$summary_table = renderTable({
    
    inFile <- input$file
    
    if (is.null(inFile))
      return(NULL)
    
    conditions = names(input_data$M)
    
    data.frame(ctrl = conditions[1:n_samples], 
               X = conditions[n_samples +  1:n_samples],
               Y = conditions[2*n_samples +  1:n_samples],
               'X+Y' = conditions[3*n_samples +  1:n_samples])
    
  })
  
  
  output$input_data_table <- renderDataTable({
    
    inFile <- input$file
    
    if (is.null(inFile)) return(NULL)
    
    input_data$M = read.table(inFile$datapath, header = TRUE, sep = "\t")
    
    head(input_data$M)
    
  },
  options=list(iDisplayLength=5,                    # initial number of records
               aLengthMenu=c(5,10),                  # records/page options
               bLengthChange=0,                       # show/hide records per page dropdown
               bFilter=0,                                    # global search box on/off
               bInfo=0,                                      # information on/off (how many records filtered, etc)
               bAutoWidth=0,                            # automatic column width calculation, disable if passing column width via aoColumnDefs
               aoColumnDefs = list(list(sWidth="300px", aTargets=c(list(0),list(1))))    # custom column size                       
  ))
  
  
  output$prof_txt <- renderText( {  paste(paste("www/case", input$case, sep = "_"),".png",sep = "") })
  
  output$case_img <- renderImage( {
    
    my_file = paste(paste("www/case", input$case, sep = "_"),".png",sep = "")
    
    list(
      src = my_file,
      contentType = 'image/png',
      width = 250,
      height = 240,
      alt = "Gene filtered out"
    )
    
  }, deleteFile = FALSE)
  
  
  degs_shown = reactive({
    
    degs_shown = filter(degs,
                  case == as.numeric(input$case),
                  type == input$interaction)
    
    #degs_shown = degs_shown[, c("genes", "adjusted_pvals", "bliss")]
    degs_shown = degs_shown[order(abs(degs_shown$bliss), decreasing = T), ]
    
    row.names(degs_shown) = NULL
    degs_shown$adjusted_pvals = round(degs_shown$adjusted_pvals, 4)
    degs_shown$bliss = round(degs_shown$bliss, 2)
    names(degs_shown)[20] = 'FDR'
    degs_shown
    
  })
  
  output$interactions_table <- renderDataTable({

    dat <- datatable(degs_shown()[1:10, c("genes", "FDR", "bliss")], 
                     options=list(iDisplayLength = 5,                    # initial number of records
                                  aLengthMenu=c(5,10),                  # records/page options
                                  bLengthChange=0,                       # show/hide records per page dropdown
                                  bFilter=0,                                    # global search box on/off
                                  bInfo=0,                                      # information on/off (how many records filtered, etc)
                                  bAutoWidth=0,                            # automatic column width calculation, disable if passing column width via aoColumnDefs
                                  aoColumnDefs = list(list(sWidth="300px", aTargets=c(list(0),list(1))))    # custom column size
                     )) %>% formatStyle('genes', color = 'white', backgroundColor = ifelse(input$interaction == 'P', 'blue', 'red'))
    
  })
  
  
  design = factor(c(rep("0", 3), rep("X", 3), rep("Y", 3), rep("Y+X", 3)))
  
  observeEvent(c(input$case,input$interaction, input$gene), {
    
    output$choose_gene <- renderUI({
      
      selectInput("gene", "Select gene:", as.character(degs_shown()[, "genes"]),
                  selected = input$gene)
      
    })
    
    log2expr = as.numeric(degs_shown()[degs_shown()$genes == input$gene, 3:14])
    
    gene_symbol = input$gene
    
      col = ifelse(input$interaction == 'P', 'blue', 'red')
      
      output$gene_plot = renderPlotly(
        ggplot(data.frame(design, log2expr), aes(x = design, y = log2expr)) +
          geom_boxplot(alpha = 0.80) +
          geom_point(colour = col, size = 2) +
          ylab('log2(expr)') + ggtitle(gene_symbol)
      )
         
      
    })
  # observeEvent(input$gene, {
  #   output$gene_plot = renderPlotly(
  #     genetab = as.numeric(degs_shown()[1, 3:14]),
  #     ggplot(data.frame(design, genetab), aes(x = design, y = genetab)) +
  #       geom_boxplot(alpha = 0.80) +
  #       geom_point(colour = 'blue', size = 2) +
  #       ylab('log2(expr)')
  #   )
  # })
  
  # observeEvent(input$gene,{
  #   
  #   #genetab = as.numeric(degs_shown_df[degs_shown_df$genes == input$gene, 3:14])
  #   genetab = as.numeric(degs_shown()[1, 3:14])
  #   
  #   names(genetab) = design
  #   
  #   output$gene_plot = renderPlotly(
  #     
  #       ggplot(data.frame(design, genetab), aes(x = design, y = genetab)) +
  #       geom_boxplot(alpha = 0.80) +
  #       geom_point(colour = 'blue', size = 2) +
  #         ylab('log2(expr)')
  #     )
  # })

  
  
  output$francois <- renderDataTable(as.data.frame(GEOmatrix))

}







#degs = degs[degs$type != 'A', ]
#degs$type = droplevels(degs$type)
#nodes = c(names(table(degs$type)), names(table(degs$case)))
#links = melt(table(degs$type, degs$case))
#links = links[links$value > 20, ]
#seqs = seq_along(nodes) - 1


# output$sankey = renderPlotly(
#  
#   #mapvalues(as.character(links$Var.1), nodes, seqs)
#   plot_ly(
#     type = "sankey",
#     orientation = "h",
# 
#     node = list(
#       label = c(names(table(degs$type)), names(table(degs$case))),
#       color = c("red", "blue"),
#       pad = 9,
#       thickness = 17,
#       line = list(
#         color = "black",
#         width = 0.1
#       )
#     ),
#     
#     link = list(
#       source = as.numeric(mapvalues(as.character(links$Var.1), nodes, seqs)),
#       target = as.numeric(mapvalues(as.character(links$Var.2), nodes, seqs)),
#       value =  links$value
#     )
#   ) 
# )