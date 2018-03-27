library(plotly)
library(enrichR)
library(DT)
library(memisc)
library(data.table)
#suppressWarnings
#suppressMessages

options(error = expression(NULL))
options(warn=-1) 

source("resolve_integration.R")
source("filter_data_on_deltas.R")
source("find_optimal_match.R")
source("visualize_all_profiles.R")
source("compute_bliss.R")
source("compute_minimum_delta.R")
source("filter_results.R")

#source("setPowerPointStyle.R")
load("IFN_TNF_1h")


immune_DB = read.csv('immune_selection/DB_immune_selection.txt', sep = '\t')
immune_DB = immune_DB[,1:12]
immune_DB = immune_DB[,-9]

#read profile definitions
PROFCODES = read.table("profile_codes_v2.txt", header = TRUE,sep = "\t") 

#DB for enrichment analysis
dbs <- c("GO_Biological_Process_2017b", "NCI-Nature_2016",
         "WikiPathways_2016", "Reactome_2016")


#degs = NULL
#suppressMessages()

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
    
    if (max(expression_data)>25) {
      
      expression_data = log2(expression_data)
    }
    
    cof = apply(expression_data, 1, function(x) sd(x)/abs(mean(x)))
    
    cof_cutoff = summary(cof)[2]
    
    cof_filter = which(cof > cof_cutoff)
    
    #expression_data[cof_filter, ]
    input_data$M[,-(1:2)] = expression_data
    
    input_data$M[cof_filter,]

  })
  
  
  observeEvent(input$start, {
    
    my.pca <- prcomp(t(as.data.frame.matrix(filtered_data()[,-(1:2)])), center = TRUE, scale=TRUE)
    
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
    
  })
  
  degs_shown = reactive({})
  
  observeEvent(input$start2, {
    
    samples = ncol(filtered_data()[,-(1:2)])/4
    design = factor(c(rep("CTRL", samples),
                      rep("X", samples),
                      rep("Y", samples),
                      rep("YX", samples)))

    results = suppressWarnings(resolve_integration(filtered_data(), design, PROFCODES))
    filtered_results = filter_results(results, adjusted_pval = 0.05)
    bliss = apply(filtered_results[[2]], 1, function(x) compute_bliss(x[2:5]))
    degs1 = filtered_results[[1]]
    degs1$bliss = bliss
    
    updateTabsetPanel(session = session, inputId = "tabs", selected = "Browse interactions")
    
  })
  
  
  
  observeEvent(input$go_to_DB, {
    
    updateTabsetPanel(session = session, inputId = "tabs", selected = "Immune X + Y Database")
    
  })
  
  
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
    
    #input_data$M = read.table(inFile$datapath, header = TRUE, sep = "\t", stringsAsFactors = F)
    input_data$M = read.csv(inFile$datapath, header = TRUE, sep = ",", stringsAsFactors = F)#, header = TRUE, sep = "\t", stringsAsFactors = F, row)
    
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
  
  
  observeEvent(input$start2, {
    
    updateTabsetPanel(session = session, inputId = "tabs", selected = "Browse interactions")
  
  })
  
  
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

    dat <- datatable(na.omit(degs_shown()[1:10, c("genes", "FDR", "bliss")]), 
                     options=list(iDisplayLength = 10,                    # initial number of records
                                  aLengthMenu=c(5,10),                  # records/page options
                                  bLengthChange=0,                       # show/hide records per page dropdown
                                  bFilter=0,                                    # global search box on/off
                                  bInfo=0,                                      # information on/off (how many records filtered, etc)
                                  scrollY = 180,
                                  paging = FALSE,
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
      
      output$gene_plot = renderPlotly({
        
        p = ggplot(data.frame(design, log2expr), aes(x = design, y = log2expr)) +
          geom_boxplot(alpha = 0.80) +
          geom_point(colour = col, size = 2) +
          ylab('log2(expr)') + ggtitle(gene_symbol)
        
        #shinyjs::delay(expr =({ 
        #  options(warn = storeWarn) 
        #}) ,ms = 10000)  
        p
      }
      )
         
      
    })
  
  
  observeEvent(input$explore1, {
    
    withProgress(message = 'Generating results', value = 0, {
      
      selected_genes = as.character(degs_shown()[, 'genes'])
      
      n <- 3
      
      for (i in 1:n) {
        
        enriched <- enrichr(selected_genes, dbs)
        
        enrich_tab = do.call("rbind", lapply(enriched, function(x) head(x)[,c("Term", "Overlap", "P.value", "Genes")]))
        
        enrich_tab$P.value = round(enrich_tab$P.value, 3)
        enrich_tab = enrich_tab[enrich_tab$P.value < 0.05, ]
        enrich_tab = enrich_tab[order(enrich_tab$P.value), ]
        

        output$enrich_tab = renderDataTable(datatable(enrich_tab,
                                            options = list(pageLength = 20)) %>% formatStyle('Term', color = 'white', backgroundColor = ifelse(input$interaction == 'P', 'blue', 'red')))
        
        incProgress(1/n, detail = paste(" ", i))
        
      }
      
    })
    
    updateTextInput(session, "explore1", value = "")
    updateTabsetPanel(session = session, inputId = "tabs", selected = "Functions/Pathways")
    
  })
  
  
  #output$francois <- renderDataTable(as.data.frame(GEOmatrix))
   
  output$immune_DB <- renderDataTable({as.data.frame(immune_DB)
  },

    options=list(iDisplayLength=50,                    # initial number of records
                 aLengthMenu=c(5,10),                  # records/page options
                 bLengthChange=0,                       # show/hide records per page dropdown
                 bFilter=0,                                    # global search box on/off
                 bInfo=0,                                      # information on/off (how many records filtered, etc)
                 bAutoWidth=0,                            # automatic column width calculation, disable if passing column width via aoColumnDefs
                 aoColumnDefs = list(list(sWidth="300px", aTargets=c(list(0),list(1))))    # custom column size
    )
    )
  
  
  output$imageGrid <- renderUI({
    images = list('case_1.png', 'case_2.png', 'case_3.png', 'case_4.png', 'case_5.png',
                  'case_6.png')
    
    #fluidRow(
      lapply(images, function(img) {
        #column(1, 
               #tags$img(src=paste0("www/", img), class="clickimg", 'data-value'=img)
               tags$img(src = img, class="clickimg", 'data-value' = img, width="150")
        #)
      })
    #)
  })
  
  
  
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