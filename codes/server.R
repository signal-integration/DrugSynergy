server = function(input, output, session) {
  
  input_data = reactiveValues()
  
  observeEvent(input$upload, {
    
    if (input$upload == 'go_to_upload') {
      
      output$main = renderUI({
        
        list(fileInput("file", label = ""),
             
             textOutput('summary_text'),
             
             dataTableOutput("input_data_table"),
             
             h5(textOutput("annotation_samples")),
          
          #h5(textOutput("control_samples")),
          
          #h5(textOutput("X_samples")),
          
          #h5(textOutput("Y_samples")),
          
          #h5(textOutput("X+Y samples")),
          
          dropdownButton(
            
            uiOutput("set_columns"),
            
            circle = TRUE,
            
            status = "danger",
            
            icon = icon("gear"),
            
            width = "300px",
            
            tooltip = tooltipOptions(title = "Define Samples")
            
          ),
          
          actionButton("start", "Quality control", icon("play"))
          
        )
      })
      
    }
    
    
    if (input$upload == 'go_to_DB') {
      
      output$main = renderUI({dataTableOutput('immune_DB')
        
        })
      
      observeEvent(input$select_button, {
        
        selectedRow = as.numeric(strsplit(input$select_button, "_")[[1]][2])
        
        data_to_read = paste0('immune_selection/', immune_db_shown$data[selectedRow, 2], '.txt')
        
        data_from_file = fread(data_to_read,
                               data.table = F,
                               stringsAsFactors = F)
        
        #updateTabsetPanel(session = session, inputId = "tabs", selected = "Upload Data")
        names(data_from_file)[2] = 'genes'
        
        names(data_from_file) = make.names(names(data_from_file), unique = T)
        
        input_data$M = data_from_file
        
        my_pca <-
          prcomp(t(as.data.frame.matrix(filtered_data()[,-(1:2)])),
                 center = TRUE,
                 scale = TRUE)
        
        samples = ncol(filtered_data()) / 4
        
        cols = c(rep("ctrl", samples),
                 rep("X", samples),
                 rep("Y", samples),
                 rep("Y+X", samples))
        
        pca_df = data.frame(PC1 = my_pca$x[, 1], PC2 = my_pca$x[, 2])
        
        pca_df$cols = cols
        
        output$pca_plot = renderPlotly(
          
          plot_ly(data = pca_df,
                  x = ~ PC1,
                  y = ~ PC2,
                  color = ~ cols,
                  marker = list(size = 15)) %>% layout(title = 'Principal Component Analysis')
          
          )
        
        updateTabsetPanel(session = session,
                          
                          inputId = "tabs",
                          
                          selected = "2. Check Quality and Run")
        
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
  
  options = list(
    iDisplayLength = 5,
    # initial number of records
    aLengthMenu = c(5, 10),
    # records/page options
    bLengthChange = 0,
    # show/hide records per page dropdown
    bFilter = 0,
    # global search box on/off
    bInfo = 0,
    # information on/off (how many records filtered, etc)
    bAutoWidth = 0,
    # automatic column width calculation, disable if passing column width via aoColumnDefs
    aoColumnDefs = list(list(
      sWidth = "300px", aTargets = c(list(0), list(1))
    ))    # custom column size
  ))
  
  
  output$summary_text = renderText({
    inFile <- input$file
    
    if (is.null(inFile))
      return(NULL)
    
    paste(
      "Uploaded file has",
      nrow(input_data$M),
      "rows",
      "and",
      
      ncol(input_data$M),
      "columns",
      "(NA=",
      sum(is.na(input_data$M)),
      ")"
    )
    
  })
  
  
  output$define = renderText({
    inFile <- input$file
    
    if (is.null(inFile))
      return(NULL)
    
    paste("Define samples")
    
  })
  
  
  output$set_columns <- renderUI({
    
    inFile <- input$file
    
    if (is.null(inFile)) return(NULL)
    
    file_columns = dim(input_data$M)[2]
    
    file_columns_names = names(input_data$M)
    
    column_list <- vector("list", file_columns)
    
    samples = ncol(filtered_data()[,-(1:2)]) / 4
    
    column_names = c(rep("Info", 2),
                     rep("CTRL", samples),
                     rep("X", samples),
                     rep("Y", samples),
                     rep("X+Y", samples)
    )
    
    for (i in 1:file_columns) {
      
      column_list[[i]] <- list(
        
        radioButtons(
          
          inputId = paste0("mVar", i),
          
          label = file_columns_names[i],
          
          choices = c("Info", "CTRL", "X", "Y", "X+Y"),
          
          inline = T,
          
          selected = column_names[i]
          
        ),
        
        br()
        
      )
    }
    
    return(column_list)
    
  })
  
  output$annotation_samples = renderText({
    inFile <- input$file
    
    if (is.null(inFile))
      return(NULL)
    
    file_columns = ncol(input_data$M)
    
    #col_names = unlist(reactiveValuesToList(input)[paste0("mVar", 1:file_columns)])[-(1:2)]
    
    #col_names = names(input_data$M)
    #samples = (ncol(input_data$M) - 2) / 4
    #paste(c("X samples:", col_names[(2 + (samples + 1)):(2 + 2 * (samples))]))
    unlist(reactiveValuesToList(input)[paste0("mVar", 3:file_columns)])
    
  })
  
  output$y_samples = renderText({
    inFile <- input$file
    
    if (is.null(inFile))
      return(NULL)
    
    #file_columns = ncol(input_data$M)
    
    #col_names = unlist(reactiveValuesToList(input)[paste0("mVar", 1:file_columns)])[-(1:2)]
    
    col_names = names(input_data$M)
    
    samples = (ncol(input_data$M) - 2) / 4
    
    paste(c("Y samples:", col_names[(2 + (2 * samples + 1)):(2 + 3 * (samples))]))
    #unlist(reactiveValuesToList(input)[paste0("mVar", 1:file_columns)])
    
  })
  
  
  output$x_samples_1 = renderText({
    inFile <- input$file
    
    if (is.null(inFile))
      return(NULL)
    
    #file_columns = ncol(input_data$M)
    
    #col_names = unlist(reactiveValuesToList(input)[paste0("mVar", 1:file_columns)])[-(1:2)]
    
    col_names = names(input_data$M)
    samples = (ncol(input_data$M) - 2) / 4
    paste(c("X samples:", col_names[(2 + (samples + 1)):(2 + 2 * (samples))]))
    #unlist(reactiveValuesToList(input)[paste0("mVar", 1:file_columns)])
    
  })
  
  output$y_samples_1 = renderText({
    inFile <- input$file
    
    if (is.null(inFile))
      return(NULL)
    
    #file_columns = ncol(input_data$M)
    
    #col_names = unlist(reactiveValuesToList(input)[paste0("mVar", 1:file_columns)])[-(1:2)]
    
    col_names = names(input_data$M)
    samples = (ncol(input_data$M) - 2) / 4
    paste(c("Y samples:", col_names[(2 + (2 * samples + 1)):(2 + 3 * (samples))]))
    #unlist(reactiveValuesToList(input)[paste0("mVar", 1:file_columns)])
    
  })
  
  
  output$x_samples_2 = renderText({
    inFile <- input$file
    
    if (is.null(inFile))
      return(NULL)
    
    #file_columns = ncol(input_data$M)
    
    #col_names = unlist(reactiveValuesToList(input)[paste0("mVar", 1:file_columns)])[-(1:2)]
    
    col_names = names(input_data$M)
    samples = (ncol(input_data$M) - 2) / 4
    paste(c("X samples:", col_names[(2 + (samples + 1)):(2 + 2 * (samples))]))
    #unlist(reactiveValuesToList(input)[paste0("mVar", 1:file_columns)])
    
  })
  
  output$y_samples_2 = renderText({
    inFile <- input$file
    
    if (is.null(inFile))
      return(NULL)
    
    #file_columns = ncol(input_data$M)
    
    #col_names = unlist(reactiveValuesToList(input)[paste0("mVar", 1:file_columns)])[-(1:2)]
    
    col_names = names(input_data$M)
    samples = (ncol(input_data$M) - 2) / 4
    paste(c("Y samples:", col_names[(2 + (2 * samples + 1)):(2 + 3 * (samples))]))
    #unlist(reactiveValuesToList(input)[paste0("mVar", 1:file_columns)])
    
  })
  
  
  observeEvent(input$start1, {
    updateTabsetPanel(session = session,
                      inputId = "tabs",
                      selected = "2. Check Quality and Run")
    
  })
  
  
  filtered_data = reactive({
    
    filter_data(input_data$M)
    
  })
  
  
  observeEvent(input$start, {
    my.pca <-
      prcomp(t(as.data.frame.matrix(filtered_data()[,-(1:2)])),
             center = TRUE,
             scale = TRUE)
    
    samples = ncol(filtered_data()) / 4
    
    cols = c(rep("ctrl", samples),
             rep("X", samples),
             rep("Y", samples),
             rep("Y+X", samples))
    
    pca_df = data.frame(PC1 = my.pca$x[, 1], PC2 = my.pca$x[, 2])
    
    pca_df$cols = cols
    
    output$pca_plot = renderPlotly(
      plot_ly(
        data = pca_df,
        x = ~ PC1,
        y = ~ PC2,
        color = ~ cols,
        marker = list(size = 15)
      ) %>% layout(title = 'Principal Component Analysis')
    )
    
    updateTabsetPanel(session = session,
                      inputId = "tabs",
                      selected = "2. Check Quality and Run")
    
  })
  
  
  observeEvent(input$start2, {
    
    samples = ncol(filtered_data()[,-(1:2)]) / 4

     design = factor(c(
       rep("CTRL", samples),
       rep("X", samples),
       rep("Y", samples),
       rep("YX", samples)
     ))
    
    combinatorial_data = filtered_data()
    
    if (input$platform == 'rnaseq') combinatorial_data = pre_process_RNA_seq(combinatorial_data, design)
    
    #file_columns = ncol(input_data$M)
    #design = factor(unlist(reactiveValuesToList(input)[paste0("mVar", 1:file_columns)])[-(1:2)])
    
    n = 4
    
    withProgress(message = 'Generating results', value = 0, {
      
      for (i in 1:n) {
        
        incProgress(1 / n, detail = paste(" ", i))
        
      }
      
      results = suppressWarnings(resolve_integration(combinatorial_data, design, PROFCODES))
      
    })
    

    #a = reactive({filter_results(results, adjusted_pval = as.numeric(input$pval))})
    
    filtered_results = filter_results(results, adjusted_pval = input$pval)

    
    bliss = apply(filtered_results[[2]], 1, function(x)
      compute_bliss(x[2:5]))
    
    degs = filtered_results[[1]]
    
    visualize_all_profiles(degs)
    
    degs$bliss = bliss
    
    #updateTextInput(session, "explore1", value = NULL)
    updateTextInput(session, "case", value = "1")
    updateTextInput(session, "interaction", value = "P")
    
    output$case_img <- renderImage({
      my_file = paste(paste("www/case", input$case, sep = "_"), ".png", sep = "")
      
      list(
        src = my_file,
        contentType = 'image/png',
        width = 250,
        height = 240,
        alt = "No data"
      )
      
    }, deleteFile = FALSE)
    
    
    output$imageGrid <- renderUI({
      
      images = as.list(paste0("case_", 1:17, '.png'))

      A = lapply(images, function(img) {

        tags$button(
          id = strsplit(img, split = '\\.')[[1]][1],
          class = "btn action-button btn-large btn-default",
          img(src = img,
              width = "120")
        )})
      
      # list(radioButtons("rb", "Choose case:",
      #                     choiceNames = list(
      #                     
      #                       tags$img(src= 'case_1.png', width="125"),
      #                       tags$img(src= 'case_2.png', width="125"),
      #                       tags$img(src= 'case_3.png', width="125"),
      #                       tags$img(src= 'case_4.png', width="125"),
      #                       tags$img(src= 'case_5.png', width="125"),
      #                       tags$img(src= 'case_6.png', width="125"),
      #                       tags$img(src= 'case_7.png', width="125"),
      #                       tags$img(src= 'case_8.png', width="125"),
      #                       tags$img(src= 'case_9.png', width="125"),
      #                       tags$img(src= 'case_10.png', width="125"),
      #                       tags$img(src= 'case_11.png', width="125"),
      #                       tags$img(src= 'case_12.png', width="125")
      #                       
      #                       ),
      #                     choiceValues = paste0("case_", 1:12), inline = T))
             
             #A)
      # list(
      #   
      #   div(style="width: 80;",
      #       sliderInput("pval",
      #                   "Significance threshold:",
      #                   min = 0.01,
      #                   max = 0.1,
      #                   value = 0.05,
      #                   step = 0.01)),
      # 
      #   div(style="width: 80;",
      #          sliderInput("minfc",
      #                      "Minimum fold-change (log2):",
      #                       min = 0.1,
      #                       max = 1,
      #                       value = 0.5,
      #                       step = 0.1)),
      #      
      #                 A)
      
      
           #div(style="width: 25px;",
#)
      # column_list = list()
      # 
      # img = paste0("case_", 1:17, '.png')
      # 
      # for (i in 1:length(img)) {
      # 
      #   column_list[[i]] = tags$button(
      #     #id = strsplit(img, split = '\\.')[[1]][1],
      #     id = paste0("case_", i),
      #         type="submit",
      #         class = "btn btn-success",
      #         img(src = img[i],
      #             width = "120")
      #   )
      #   }
      # 
      # return(column_list)
      
      })
      
    #})
    
    #answer <- reactiveValues()
    #answer$text = unlist(reactiveValuesToList(input)["case_1"])
    #output$grid_text = renderText(answer$text)
    

    observeEvent(input$case_1, {
      
      output$case_img_new <- renderImage({
        
        my_file = paste0("www/", paste("case_", "1", ".png", sep = ""))
        
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
      #answer$text = my_file = paste0("www/", paste("case_", "2", ".png", sep = ""))

      output$case_img_new <- renderImage({
        my_file = paste0("www/", paste("case_", "2", ".png", sep = ""))

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
      
      my_file = paste0("www/", paste("case_", "3", ".png", sep = ""))

      output$case_img_new <- renderImage({
        my_file = paste0("www/", paste("case_", "3", ".png", sep = ""))

        list(
          src = my_file,
          contentType = 'image/png',
          width = 250,
          height = 240,
          alt = "No data"
        )

      }, deleteFile = FALSE)

    })
    # 
    # observeEvent(input$case_4, {
    #   answer$text = "case 4"
    #   
    # })
    # 
    # observeEvent(input$case_5, {
    #   #answer$data <- get.focus.spec(input=input, premise=premise,
    #   #                            itemname=input$dropdown.itemname, spec.info=spec.info)
    #   answer$text = "case 5"
    #   
    # })
    
    
    
    #output$click_case = renderText(answer$text)
    
    
    # observeEvent(input$case_2, {
    #
    #   output$click_case2 = renderText("case 2!")
    #
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
    
    
    #observeEvent(c(input$start2, input$case, input$interaction, input$gene),
    observeEvent(c(input$start2, input$case, input$interaction, input$gene),
                 {
                   #degs_shown = filter(degs, case == as.numeric(input$case), type == input$interaction)
                   
                   save(degs, file = paste0('results', input$select_button))
                   
                   degs_shown = filter(degs,
                                       case == as.numeric(input$case),
                                       type == input$interaction)
                   
                   if (nrow(degs_shown) == 0) {
                     degs_shown = NULL
                     
                     output$interactions_table = renderDataTable({
                       
                     })
                     
                     output$choose_gene = renderUI({
                       
                     })
                     
                     output$gene_plot = renderPlotly(plotly_empty())
                     
                   }
                   
                   else{
                     
                     degs_shown = degs_shown[order(abs(degs_shown$bliss), decreasing = T), ]
                     
                     row.names(degs_shown) = NULL
                     
                     degs_shown$adjusted_pvals = round(degs_shown$adjusted_pvals, 4)
                     
                     degs_shown$bliss = round(degs_shown$bliss, 2)
                     
                     names(degs_shown)[20] = 'FDR'
                     
                     
                     urls = paste0('http://www.genecards.org/cgi-bin/carddisp.pl?gene=',
                                   degs_shown$genes)
                     
                     degs_shown$gene = paste0("<a href='",
                                              urls,
                                              "' target='_blank'>",
                                              degs_shown$genes,
                                              " </a>")
                     
                     output$interactions_table <- renderDataTable({
                       dat <-
                         datatable(
                           na.omit(degs_shown[, c("gene", "FDR", "bliss")]),
                           escape = FALSE,
                           
                           options = list(
                             iDisplayLength = 20,
                             # initial number of records
                             aLengthMenu = c(5, 10),
                             # records/page options
                             bLengthChange = 0,
                             # show/hide records per page dropdown
                             bFilter = 0,
                             # global search box on/off
                             bInfo = 0,
                             # information on/off (how many records filtered, etc)
                             scrollY = 180,
                             paging = FALSE,
                             bAutoWidth = 0,
                             # automatic column width calculation, disable if passing column width via aoColumnDefs
                             aoColumnDefs = list(list(
                               sWidth = "300px",
                               aTargets = c(list(0), list(1))
                             ))    # custom column size
                           )
                         ) #%>% formatStyle('urls', color = 'white', backgroundColor = ifelse(input$interaction == 'P', 'blue', 'red'))
                       
                     })
                     
                     output$choose_gene <- renderUI({
                       selectInput("gene",
                                   "Select gene:",
                                   choices = degs_shown[, "genes"],
                                   selected = input$gene)
                       
                     })
                     
                     gene_symbol = input$gene
                     
                     samples = ncol(filtered_data()[,-(1:2)]) / 4
                     
                     design = factor(c(
                       rep("CTRL", samples),
                       
                       rep("X", samples),
                       
                       rep("Y", samples),
                       
                       rep("YX", samples)
                     ))
                     
                     log2expr = as.numeric(degs_shown[(degs_shown$genes == gene_symbol), 3:(2 + length(design))])
                     
                     col = ifelse(input$interaction == 'P', 'blue', 'red')
                     
                     output$gene_plot = renderPlotly({
                       p = ggplot(data.frame(design, log2expr),
                                  aes(x = design, y = log2expr)) +
                         geom_boxplot(alpha = 0.80) +
                         geom_point(colour = col, size = 2) +
                         ylab('log2(expr)') + ggtitle(gene_symbol)  + theme_bw()            
                       ggplotly(p) %>% config(displayModeBar = F)
                       
                     })
                     
                   }
                   
                   
                   updateTabsetPanel(session = session,
                                     inputId = "tabs",
                                     selected = "3. Browse Results")
                   
                   
                   
                   #updateTextInput(session, "select_button", value = "")
                   
                   #                   updateTabsetPanel(session = session,
                   #                                     inputId = "tabs",
                   #                                     selected = "3. Browse Results")
                   
                   
                 })
    
    observeEvent(input$explore1, {
      degs_shown = filter(degs,
                         case == as.numeric(input$case),
                         type == input$interaction)
      
      selected_genes = as.character(degs_shown[, 'genes'])
      
      
      withProgress(message = 'Generating results (results below)', value = 0, {
        
        for (i in 1:n) {
      
      incProgress(1 / n, detail = paste(" ", i))
      
        }
        enriched <- enrichr(selected_genes, dbs)
      })
      
      enrich_tab = do.call("rbind", lapply(enriched, function(x)
        head(x)[, c("Term", "Overlap", "P.value", "Genes")]))
      
      enrich_tab$P.value = round(enrich_tab$P.value, 3)
      
      enrich_tab = enrich_tab[enrich_tab$P.value < 0.05, ]
      
      enrich_tab = enrich_tab[order(enrich_tab$P.value), ]
      output$enrich_tab = renderDataTable(
        datatable(enrich_tab,
                  options = list(pageLength = 20)) %>% formatStyle(
                    'Term',
                    color = 'white',
                    backgroundColor = ifelse(input$interaction == 'P', 'blue', 'red')
                  )
      )
      
      
      #updateTextInput(session, "explore1", value = NULL)

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
    Actions = shinyInput(
      actionButton,
      nrow(immune_DB),
      'button_',
      label = "Analyze",
      onclick = 'Shiny.onInputChange(\"select_button\",  this.id)'
    ),
    immune_DB,
    row.names = 1:nrow(immune_DB)
  ))
  
  output$immune_DB <- DT::renderDataTable(
    immune_db_shown$data,
    server = FALSE,
    escape = FALSE,
    selection = 'none',
    options = list(
      pageLength = 35)
    #pageLength = 35
    
  )
  
  
  geo_dataset = reactive({
    NULL
  })
  
}





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