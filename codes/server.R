server = function(input, output, session) {
  
  #upload from file
#  observeEvent(input$go_to_upload, {
    
#    updateTabsetPanel(session = session, inputId = "tabs", selected = "Upload Data")}
    
#  )
  
  #upload from DB
#  observeEvent(input$go_to_DB, {
    
#    updateTabsetPanel(session = session, inputId = "tabs", selected = "Immune X + Y Database")
    
#  })
  
  #upload from GEO
#  observeEvent(input$go_to_GEO, {
    
#    updateTabsetPanel(session = session, inputId = "tabs", selected = "X")}
    
#  )
  
  #Go to demo
  observeEvent(input$go_to_demo, {
    
    updateTabsetPanel(session = session, inputId = "tabs", selected = "Demo")}
    
  )
  
  
  #import data from file
  input_data = reactiveValues()
  
  
  observeEvent(input$upload, {
    
    if (input$upload == 'go_to_upload'){
      
      output$main = renderUI({
        
        list(fileInput("file", label = ""),
        
             textOutput('summary_text'),
             
             dataTableOutput("input_data_table"),
             
             actionButton("start", "Quality control", icon("play")))})
      
      }
    
    
    if (input$upload == 'go_to_DB'){
      
      output$main = renderUI({
        
        dataTableOutput('immune_DB')
        
      })
      
      observeEvent(input$select_button, {
        
        selectedRow <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
        
        data_to_read <<- paste0('immune_selection/',immune_db_shown$data[selectedRow, 2], '.txt')
        
        data_from_file = fread(data_to_read, data.table = F, stringsAsFactors = F)
        #updateTabsetPanel(session = session, inputId = "tabs", selected = "Upload Data")
        names(data_from_file)[2] = 'genes'
        
        names(data_from_file) = make.names(names(data_from_file), unique = T)
        
        input_data$M = data_from_file
        
        my.pca <- prcomp(t(as.data.frame.matrix(filtered_data()[,-(1:2)])), center = TRUE, scale = TRUE)
        
        samples = ncol(filtered_data())/4
        
        cols = c(rep("ctrl", samples), rep("X", samples),
                 rep("Y", samples), rep("Y+X", samples))
        
        pca_df = data.frame(PC1 = my.pca$x[, 1], PC2 = my.pca$x[, 2])
        
        pca_df$cols = cols
        
        output$pca_plot = renderPlotly(plot_ly(data = pca_df, x = ~PC1, y = ~PC2, 
                                               color = ~cols, marker = list(size = 15)) %>% layout(title = 'Principal Component Analysis'))
        
        updateTabsetPanel(session = session, inputId = "tabs", selected = "2. Check Quality and Run")
        
        
        
      })
      
    
      
      }
    
    
   })
  
  
  
  output$input_data_table <- renderDataTable({
    
    inFile <- input$file
    
    if (is.null(inFile)) return(NULL)
    
    data_to_read = inFile$datapath

    #data_to_read = geo_dataset

    #data_from_file = fread(inFile$datapath, data.table = F, stringsAsFactors = F)
     
    #if (length(geo_dataset) > 0)  {
      
    data_from_file = fread(data_to_read, data.table = F, stringsAsFactors = F)
  
    names(data_from_file)[2] = 'genes'
    
    names(data_from_file) = make.names(names(data_from_file), unique = T)
    
    input_data$M = data_from_file
    
    head(input_data$M, 10)
    
  },
  
  options=list(iDisplayLength=5,                    # initial number of records
               aLengthMenu=c(5,10),                  # records/page options
               bLengthChange=0,                       # show/hide records per page dropdown
               bFilter=0,                                    # global search box on/off
               bInfo=0,                                      # information on/off (how many records filtered, etc)
               bAutoWidth=0,                            # automatic column width calculation, disable if passing column width via aoColumnDefs
               aoColumnDefs = list(list(sWidth="300px", aTargets=c(list(0),list(1))))    # custom column size                       
  ))
  
  
  output$summary_text = renderText({
    
    inFile <- input$file
    
    if (is.null(inFile)) return(NULL)
    
    paste("Uploaded file has", nrow(input_data$M), "rows", "and",
          
          ncol(input_data$M), "columns", "(NA=", sum(is.na(input_data$M)), ")")
    
  })
  
  
  output$define = renderText({
    
    inFile <- input$file
    
    if (is.null(inFile)) return(NULL)
    
    paste("Define samples")
    
  })
  
  
  output$set_columns <- renderUI({
    
    inFile <- input$file
    
    if (is.null(inFile)) return(NULL)
    
    file_columns = dim(input_data$M)[2]
    
    file_columns_names = names(input_data$M)

    column_list <- vector("list", file_columns)        
    
    for(i in 1:file_columns){
      
      column_list[[i]] <- list(radioButtons(inputId = paste0("mVar",i), 
                                   
                                   label = file_columns_names[i], 
                                   
                                   choices = c("Info", "CTRL", "X", "Y", "X+Y"),
                                   
                                   inline = T),
                      
                      br()
                      
      )}      
    
    return(column_list)                      
    
  })
  
  
  output$test = renderText({
    
    inFile <- input$file
    
    if (is.null(inFile)) return(NULL)
    
    file_columns = ncol(input_data$M)
    
    unlist(reactiveValuesToList(input)[paste0("mVar", 1:file_columns)])

  })
  
  
  observeEvent(input$start1, {
    
    updateTabsetPanel(session = session, inputId = "tabs", selected = "2. Check Quality and Run")
    
  })
  
  
  filtered_data = reactive({
    
    filter_data(input_data$M)
    
  })
  
  
  observeEvent(input$start, {
    
    my.pca <- prcomp(t(as.data.frame.matrix(filtered_data()[,-(1:2)])), center = TRUE, scale = TRUE)
    
    samples = ncol(filtered_data())/4

    cols = c(rep("ctrl", samples), rep("X", samples),
             rep("Y", samples), rep("Y+X", samples))
    
    pca_df = data.frame(PC1 = my.pca$x[, 1], PC2 = my.pca$x[, 2])
    
    pca_df$cols = cols
    
    output$pca_plot = renderPlotly(plot_ly(data = pca_df, x = ~PC1, y = ~PC2, 
                                           color = ~cols, marker = list(size = 15)) %>% layout(title = 'Principal Component Analysis'))
    
    updateTabsetPanel(session = session, inputId = "tabs", selected = "2. Check Quality and Run")
     
    })
  
  
  observeEvent(input$start2, {
    
    samples = ncol(filtered_data()[,-(1:2)])/4
    
    design = factor(c(rep("CTRL", samples),
                      rep("X", samples),
                      rep("Y", samples),
                      rep("YX", samples)))

    n = 3
    
    withProgress(message = 'Generating results', value = 0, {
      
      for (i in 1:n) {
        
        results = suppressWarnings(resolve_integration(filtered_data(), design, PROFCODES))
        
        incProgress(1/n, detail = paste(" ", i))
        
        }
      
      })
    
    filtered_results = filter_results(results, adjusted_pval = 0.05)
    
    bliss = apply(filtered_results[[2]], 1, function(x) compute_bliss(x[2:5]))
    
    degs = filtered_results[[1]]
    
    visualize_all_profiles(degs)
    
    degs$bliss = bliss
    
    updateTextInput(session, "case", value = "1")
    
    output$case_img <- renderImage( {
      
      my_file = paste(paste("www/case", input$case, sep = "_"),".png",sep = "")
      
      list(
        src = my_file,
        contentType = 'image/png',
        width = 250,
        height = 240,
        alt = "No data"
      )
      
    }, deleteFile = FALSE)

    output$imageGrid <- renderUI({
      
      images = list('case_1.png', 'case_2.png', 'case_3.png', 'case_4.png', 'case_5.png',
                    'case_6.png', 'case_7.png', 'case_8.png', 'case_9.png', 'case_10.png',
                    'case_11.png', 'case_12.png', 'case_13.png', 'case_14.png',
                    'case_15.png', 'case_16.png', 'case_17.png')
      
      #fluidRow(
      lapply(images, function(img) {
        #column(1, 
        #tags$img(src=paste0("www/", img), class="clickimg", 'data-value'=img)
        #tags$img(src = img, class="clickimg", 'data-value' = img, width="150")
        
        tags$button(
          #id = "web_button",
          id = strsplit(img, split = '\\.')[[1]][1],
          #class = "btn action_button",
          class="btn action-button btn-large btn-primary",
          img(src = img,
              width="120"),
          tags$style(HTML('color: #4d3a7d;'))
        )
        
        
        #)
      })
      #)
    })
    
    answer <- reactiveValues()
    
    observeEvent(input$case_1, {
      
      answer$text = paste0("www/", paste("case_", input$case_1,".png",sep = ""))
      
      output$case_img_new <- renderImage( {
        
        my_file = paste0("www/", paste("case_", input$case_1,".png",sep = ""))
        
        list(
          src = my_file,
          contentType = 'image/png',
          width = 250,
          height = 240,
          alt = "No data"
        )
        
      }, deleteFile = FALSE)
    
    })
    
    observeEvent(input$case_2, {
      
      answer$text = my_file = paste0("www/", paste("case_", input$case_2,".png",sep = ""))
      
      output$case_img_new <- renderImage( {
        
        my_file = paste0("www/", paste("case_", input$case_2,".png",sep = ""))
        
        list(
          src = my_file,
          contentType = 'image/png',
          width = 250,
          height = 240,
          alt = "No data"
        )
        
      }, deleteFile = FALSE)
    })
    
    observeEvent(input$case_3, {

      answer$text = my_file = paste0("www/", paste("case_", input$case_3,".png",sep = ""))
      
      output$case_img_new <- renderImage( {
        
        my_file = paste0("www/", paste("case_", input$case_3,".png",sep = ""))
        
        list(
          src = my_file,
          contentType = 'image/png',
          width = 250,
          height = 240,
          alt = "No data"
        )
        
      }, deleteFile = FALSE)
      
    })
    
    observeEvent(input$case_4, {
      
      answer$text = "case 4"
      
    })
    
    observeEvent(input$case_5, {
      #answer$data <- get.focus.spec(input=input, premise=premise, 
      #                            itemname=input$dropdown.itemname, spec.info=spec.info)
      answer$text = "case 5"
      
    })
    


    output$click_case = renderText(answer$text)
      

    # observeEvent(input$case_2, {
    #   
    #   output$click_case2 = renderText("case 2!")
    #   
    # })
    
    output$imageGrid1 <- renderUI({
      images = list('case_7.png', 'case_8.png', 'case_9.png', 'case_10.png', 'case_11.png',
                    'case_12.png')
      
      #fluidRow(
      lapply(images, function(img) {
        #column(1, 
        #tags$img(src=paste0("www/", img), class="clickimg", 'data-value'=img)
        tags$img(src = img, class="clickimg", 'data-value' = img, width="150")
        #)
      })
      #)
    })
    
    
    observeEvent(c(input$start2, input$case, input$interaction, input$gene), {
      
      
      degs_shown = filter(degs, case == as.numeric(input$case), type == input$interaction)
      
      if (nrow(degs_shown)== 0){
        
        degs_shown = NULL
        
        output$interactions_table = renderDataTable({})
        
        output$choose_gene = renderUI({})
        
        output$gene_plot = renderPlotly(plotly_empty())
        
      } 
      
      else{
        
        degs_shown = degs_shown[order(abs(degs_shown$bliss), decreasing = T), ]
        
        row.names(degs_shown) = NULL
        
        degs_shown$adjusted_pvals = round(degs_shown$adjusted_pvals, 4)
        
        degs_shown$bliss = round(degs_shown$bliss, 2)
        
        names(degs_shown)[20] = 'FDR'
        
        urls = paste0('http://www.genecards.org/cgi-bin/carddisp.pl?gene=', degs_shown$genes)
        
        degs_shown$urls = paste0("<a href='", urls,"' target='_blank'>", degs_shown$genes," </a>")
        
        output$interactions_table <- renderDataTable({
          
          dat <- datatable(na.omit(degs_shown[1:10, c("urls", "FDR", "bliss")]), escape = FALSE,
                           
                           options=list(iDisplayLength = 10,                    # initial number of records
                                        aLengthMenu=c(5,10),                  # records/page options
                                        bLengthChange=0,                       # show/hide records per page dropdown
                                        bFilter=0,                                    # global search box on/off
                                        bInfo=0,                                      # information on/off (how many records filtered, etc)
                                        scrollY = 180,
                                        paging = FALSE,
                                        bAutoWidth=0,                            # automatic column width calculation, disable if passing column width via aoColumnDefs
                                        aoColumnDefs = list(list(sWidth="300px", aTargets=c(list(0),list(1))))    # custom column size
                           )) #%>% formatStyle('urls', color = 'white', backgroundColor = ifelse(input$interaction == 'P', 'blue', 'red'))
          
        })
        
        output$choose_gene <- renderUI({
          
          selectInput("gene", "Select gene:", choices = degs_shown[, "genes"], selected = input$gene)
          
        })
        
        gene_symbol = input$gene
        
        samples = ncol(filtered_data()[,-(1:2)])/4
        
        design = factor(c(rep("CTRL", samples),
                          
                          rep("X", samples),
                          
                          rep("Y", samples),
                          
                          rep("YX", samples)))
        
        log2expr = as.numeric(degs_shown[(degs_shown$genes == gene_symbol), 3:(2 + length(design))])
        
        col = ifelse(input$interaction == 'P', 'blue', 'red')
        
        output$gene_plot = renderPlotly({
          
          p = ggplot(data.frame(design, log2expr), aes(x = design, y = log2expr)) + 
            geom_boxplot(alpha = 0.80) +
            geom_point(colour = col, size = 2) +
            ylab('log2(expr)') + ggtitle(gene_symbol)
          
          p
        })
        
      } 
      
    updateTextInput(session, "select_button", value = "")

    updateTabsetPanel(session = session, inputId = "tabs", selected = "3. Browse Results")
    
  })
     
  
     observeEvent(input$explore1, {
       
       degs_shown = filter(degs,
                           case == as.numeric(input$case),
                           type == input$interaction)
       
       withProgress(message = 'Generating results (see below when finished)', value = 0, {
         
         selected_genes = as.character(degs_shown[, 'genes'])
         
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
       #updateTabsetPanel(session = session, inputId = "tabs", selected = "Functions/Pathways")
       
     })
     
  })   

  
  myValue <- reactiveValues(geo_dataset = '')
  
       
  shinyInput <- function(FUN, len, id, ...) {
    
    inputs <- character(len)
    
    for (i in seq_len(len)) {
      
      inputs[i] <- as.character(FUN(paste0(id, i), ...))
    }
    
    inputs
  }
  
  
  immune_db_shown <- reactiveValues(data = data.frame(
    
    Actions = shinyInput(actionButton, nrow(immune_DB), 'button_', 
                         label = "Fire", 
                         onclick = 'Shiny.onInputChange(\"select_button\",  this.id)'),
    immune_DB, row.names = 1:nrow(immune_DB))
    
  )
  
  output$immune_DB <- DT::renderDataTable(
    
    immune_db_shown$data, server = FALSE, escape = FALSE, selection = 'none'
    
    )
  
  
  geo_dataset = reactive({NULL})
  
#  observeEvent(input$select_button, {

#    selectedRow <- as.numeric(strsplit(input$select_button, "_")[[1]][2])

#    geo_dataset <<- paste0('immune_selection/',immune_db_shown$data[selectedRow, 2], '.txt')

#    updateTabsetPanel(session = session, inputId = "tabs", selected = "Upload Data")
  
#    })

  # output$myText <- renderText({
  # 
  #   myValue$geo_dataset
  # 
  # })
  
  
  # degs_shown = reactive({
  # 
  #   #load("results")
  # 
  #   degs_shown = filter(degs,
  #                       case == as.numeric(input$case),
  #                       type == input$interaction)
  # 
  #   #degs_shown = degs_shown[, c("genes", "adjusted_pvals", "bliss")]
  #   degs_shown = degs_shown[order(abs(degs_shown$bliss), decreasing = T), ]
  # 
  #   row.names(degs_shown) = NULL
  #   degs_shown$adjusted_pvals = round(degs_shown$adjusted_pvals, 4)
  #   degs_shown$bliss = round(degs_shown$bliss, 2)
  #   names(degs_shown)[20] = 'FDR'
  #   degs_shown
  # 
  # })
  

  
  # output$summary_table = renderTable({
  #   
  #   inFile <- input$file
  #   
  #   if (is.null(inFile))
  #     return(NULL)
  #   
  #   conditions = names(input_data)
  #   
  #   data.frame(ctrl = conditions[1:n_samples], 
  #              X = conditions[n_samples +  1:n_samples],
  #              Y = conditions[2*n_samples +  1:n_samples],
  #              'X+Y' = conditions[3*n_samples +  1:n_samples])
  #   
  # })
  
  
  
  
  #output$prof_txt <- renderText( {  paste(paste("www/case", input$case, sep = "_"),".png",sep = "") })
  
  
  
  #design = factor(c(rep("0", 3), rep("X", 3), rep("Y", 3), rep("Y+X", 3)))

  
  
  
  
  # output$immune_DB <- renderDataTable({immune_db_shown
  # },
  # 
  #   options=list(iDisplayLength=50,                    # initial number of records
  #                aLengthMenu=c(5,10),                  # records/page options
  #                bLengthChange=0,                       # show/hide records per page dropdown
  #                #bFilter=0,                                    # global search box on/off
  #                bInfo=0,                                      # information on/off (how many records filtered, etc)
  #                bAutoWidth=0,                            # automatic column width calculation, disable if passing column width via aoColumnDefs
  #                aoColumnDefs = list(list(sWidth="300px", aTargets=c(list(0),list(1))))    # custom column size
  #   )
  #   )
  
  
   
  # output$immune_DB <- renderDataTable({as.data.frame(immune_DB)
  # },
  # 
  #   options=list(iDisplayLength=50,                    # initial number of records
  #                aLengthMenu=c(5,10),                  # records/page options
  #                bLengthChange=0,                       # show/hide records per page dropdown
  #                #bFilter=0,                                    # global search box on/off
  #                bInfo=0,                                      # information on/off (how many records filtered, etc)
  #                bAutoWidth=0,                            # automatic column width calculation, disable if passing column width via aoColumnDefs
  #                aoColumnDefs = list(list(sWidth="300px", aTargets=c(list(0),list(1))))    # custom column size
  #   )
  #   )
  
  
  # output$imageGrid <- renderUI({
  #   images = list('case_1.png', 'case_2.png', 'case_3.png', 'case_4.png', 'case_5.png',
  #                 'case_6.png')
  #   
  #   #fluidRow(
  #     lapply(images, function(img) {
  #       #column(1, 
  #              #tags$img(src=paste0("www/", img), class="clickimg", 'data-value'=img)
  #              tags$img(src = img, class="clickimg", 'data-value' = img, width="150")
  #       #)
  #     })
  #   #)
  # })
  
  # output$imageGrid1 <- renderUI({
  #   images = list('case_7.png', 'case_8.png', 'case_9.png', 'case_10.png', 'case_11.png',
  #                 'case_12.png')
  #   
  #   #fluidRow(
  #   lapply(images, function(img) {
  #     #column(1, 
  #     #tags$img(src=paste0("www/", img), class="clickimg", 'data-value'=img)
  #     tags$img(src = img, class="clickimg", 'data-value' = img, width="150")
  #     #)
  #   })
  #   #)
  # })
  
  
  
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