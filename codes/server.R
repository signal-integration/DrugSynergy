library(plotly)
library(enrichR)

source("filter_data_on_deltas.R")
source("find_optimal_match.R")
source("visualize_all_profiles.R")
source("setPowerPointStyle.R")

load("GEOmatrix (1)")
GEOmatrix[,1] = make.names(GEOmatrix[,1], unique = T)
GEOmatrix = GEOmatrix[,-2]
row.names(GEOmatrix) = NULL

#load("my_data_filtered_matched_2.5h")
#load("GSE39292_results")
#my_results_2.5h = GSE39292_results

load("GSE43452_filtered_matched")
my_results_2.5h = my_data_filtered_matched

genes = sort(as.character(my_results_2.5h$gene))


dbs <- c("GO_Biological_Process_2017b", "NCI-Nature_2016",
         "WikiPathways_2016", "Reactome_2016")


codes = read.table('profile_codes_v2.txt', header = TRUE)

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
    
    paste(
      
      "Uploaded file has", dim(input_data$M)[1], "rows", "and",
      
      dim(input_data$M)[2], "columns", "(NA=", sum(is.na(input_data)), ")"
      
    )
    
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
    
    #expression_data = as.data.frame.matrix(input_data$M[,-(1:2)])
    
    #if (max(expression_data)>25) expression_data = log2(expression_data)
    
    #removing uninformative probes (very small coefficient of variation)
    #cof_cutoff = 0.05
    
    #cof = apply(expression_data, 1, function(x) sd(x)/abs(mean(x)))
    
    #cof_filter = which(cof > cof_cutoff)
    
    #my.pca <- prcomp(t(expression_data[cof_filter, ]), center = TRUE, scale=TRUE)
    
    my.pca <- prcomp(t(filtered_data()), center = TRUE, scale=TRUE)
    
    #we assume the same number of samples for each condition
    #samples = ncol(expression_data)/4
    
    samples = ncol(filtered_data_data())/4
    

    cols = c(rep("ctrl", samples), rep("X", samples),
             rep("Y", samples), rep("X+Y", samples))
    
    
    pca_df = data.frame(PC1 = my.pca$x[, 1], PC2 = my.pca$x[, 2])
    
    row.names(pca_df) = names(input_data$M)[-(1:2)]
    
    pca_df$cols = cols
    
    output$pca_plot = renderPlotly(plot_ly(data = pca_df, x = ~PC1, y = ~PC2, 
            color = ~cols, marker = list(size = 15)) %>% layout(title = 'Principal Component Analysis'))
    
    updateTabsetPanel(session = session, inputId = "tabs", selected = "Quality Control")
     }
    )
  
  observeEvent(input$start1, {
    
    updateTabsetPanel(session = session, inputId = "tabs", selected = "Run analysis")
    
  }
  )
  
  
  
  #output$pca_plot = renderPlotly({
    
  #  pca_plot
    #inFile <- input$file
    
    #if (is.null(inFile)) return(NULL)
    
    #if input$start 
  #})
  
  
  output$summary_table = renderTable({
    
    inFile <- input$file
    
    if (is.null(inFile))
      return(NULL)
    
    #n_samples = (dim(input_data$M)[2] - 2)/4
    
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
      width = 260,
      height = 250,
      alt = "Gene filtered out"
    )
    
  }, deleteFile = FALSE)
  
  
  output$sankey = renderPlotly(
    
    plot_ly(
      type = "sankey",
      orientation = "h",
      
      node = list(
        label = c("synergistic", "antagonistic", "case 1", "case 2", "case 3"),
        color = c("blue", "red", "white", "white", "white"),
        pad = 15,
        thickness = 20,
        line = list(
          color = "black",
          width = 0.5
        )
      ),
      
      link = list(
        source = c(0,0,0, 1,1,1),
        target = c(2,3,4, 2,3,4),
        value =  c(8,1,4, 8,2, 2)
      )
    ) %>% 
      layout(
        title = "Interactions",
        font = list(
          size = 10
        )
      )
  )
  
 # output$interaction_table = renderDataTable({
    
#    my_data_filtered_matched[my_data_filtered_matched$]
#  })
  
  
  output$prof_1_img <- renderImage( {
    
    index = which(codes$case == input$case)
    
    my_file = paste(paste("www/profile", index[1], sep = "_"),".png",sep = "")
    
    list(
      src = my_file,
      contentType = 'image/png',
      width = 200,
      height = 200,
      alt = "Gene filtered out"
    )
  }, deleteFile = FALSE)
  
  
  output$prof_2_img <- renderImage( {
    
    index = which(codes$case == input$case)
    
    my_file = paste(paste("www/profile", index[2], sep = "_"),".png",sep = "")
    
    list(
      src = my_file,
      contentType = 'image/png',
      width = 200,
      height = 200,
      alt = "Gene filtered out"
    )
  }, deleteFile = FALSE)
  
  
  output$prof_3_img <- renderImage( {
    
    index = which(codes$case == input$case)
    
    my_file = paste(paste("www/profile", index[3], sep = "_"),".png",sep = "")
    
    list(
      src = my_file,
      contentType = 'image/png',
      width = 200,
      height = 200,
      alt = " "
    )
  }, deleteFile = FALSE)
  
  
  output$prof_4_img <- renderImage( {
    
    index = which(codes$case == input$case)
    
    my_file = paste(paste("www/profile", index[4], sep = "_"),".png",sep = "")
    
    list(
      src = my_file,
      contentType = 'image/png',
      width = 200,
      height = 200,
      alt = " "
    )
  }, deleteFile = FALSE)
  
  
  output$prof_5_img <- renderImage( {
    
    index = which(codes$case == input$case)
    
    my_file = paste(paste("www/profile", index[5], sep = "_"),".png",sep = "")
    
    list(
      src = my_file,
      contentType = 'image/png',
      width = 200,
      height = 200,
      alt = " "
    )
  }, deleteFile = FALSE)
  
  
  output$prof_6_img <- renderImage( {
    
    index = which(codes$case == input$case)
    
    my_file = paste(paste("www/profile", index[6], sep = "_"),".png",sep = "")
    
    list(
      src = my_file,
      contentType = 'image/png',
      width = 200,
      height = 200,
      alt = " "
    )
  }, deleteFile = FALSE)
  
  
  
  output$prof_7_img <- renderImage( {
    
    index = which(codes$case == input$case)
    
    my_file = paste(paste("www/profile", index[7], sep = "_"),".png",sep = "")
    
    list(
      src = my_file,
      contentType = 'image/png',
      width = 200,
      height = 200,
      alt = " "
    )
  }, deleteFile = FALSE)
  
  
  
  output$prof_8_img <- renderImage( {
    
    index = which(codes$case == input$case)
    
    my_file = paste(paste("www/profile", index[8], sep = "_"),".png",sep = "")
    
    list(
      src = my_file,
      contentType = 'image/png',
      width = 200,
      height = 200,
      alt = " "
    )
  }, deleteFile = FALSE)
  
  
  output$prof_9_img <- renderImage( {
    
    index = which(codes$case == input$case)
    
    my_file = paste(paste("www/profile", index[9], sep = "_"),".png",sep = "")
    
    list(
      src = my_file,
      contentType = 'image/png',
      width = 200,
      height = 200,
      alt = " "
    )
  }, deleteFile = FALSE)
  
  
  observeEvent(input$explore1, {
    
    index = which(codes$case == input$case & codes$outcome == 1)
    
    tab = my_results_2.5h[my_results_2.5h$prof_index == index, ]

    withProgress(message = 'Generating results', value = 0, {
        
        n <- 3
        
        for (i in 1:n) {
          
            enriched <- enrichr(as.character(tab$gene), dbs)
            
            enrich_tab = do.call("rbind", lapply(enriched, function(x) head(x)[,c("Term", "Overlap", "P.value", "Genes")]))
            
            enrich_tab = enrich_tab[order(enrich_tab$P.value), ]
            
            output$enrich_tab = renderDataTable(enrich_tab,
                                           options = list(pageLength = 100))
            
            incProgress(1/n, detail = paste(" ", i))
            
          }
          
        })
      
      updateTextInput(session, "explore1", value = "")
      updateTabsetPanel(session = session, inputId = "tabs", selected = "Function/Pathways")
 
      })
  
  
  
  observeEvent(input$explore2, {
    
    index = which(codes$case == input$case & codes$outcome == 2)
    
    tab = my_results_2.5h[my_results_2.5h$prof_index == index, ]
    
    withProgress(message = 'Generating results', value = 0, {
      
      n <- 3
      
      for (i in 1:n) {
        
        enriched <- enrichr(as.character(tab$gene), dbs)
        
        enrich_tab = do.call("rbind", lapply(enriched, function(x) head(x)[,c("Term", "Overlap", "P.value", "Genes")]))
        
        enrich_tab = enrich_tab[order(enrich_tab$P.value), ]
        
        output$enrich_tab = renderDataTable(enrich_tab,
                                            options = list(pageLength = 100))
        
        incProgress(1/n, detail = paste(" ", i))
        
        }
      
      
    })
    
    updateTextInput(session, "explore2", value = "")
    updateTabsetPanel(session = session, inputId = "tabs", selected = "Function/Pathways")
    
  })
  
  
  
  observeEvent(input$explore3, {
    
    index = which(codes$case == input$case & codes$outcome == 3)
    
    tab = my_results_2.5h[my_results_2.5h$prof_index == index, ]
    
    withProgress(message = 'Generating results', value = 0, {
      
      n <- 3
      
      for (i in 1:n) {
        
        enriched <- enrichr(as.character(tab$gene), dbs)
        
        enrich_tab = do.call("rbind", lapply(enriched, function(x) head(x)[,c("Term", "Overlap", "P.value", "Genes")]))
        
        enrich_tab = enrich_tab[order(enrich_tab$P.value), ]
        
        output$enrich_tab = renderDataTable(enrich_tab,
                                            options = list(pageLength = 100))
        
        incProgress(1/n, detail = paste(" ", i))
        
      }
      
      
    })
    
    updateTextInput(session, "explore3", value = "")
    updateTabsetPanel(session = session, inputId = "tabs", selected = "Function/Pathways")
    
  })
  
  
  
  observeEvent(input$explore4, {
    
    index = which(codes$case == input$case & codes$outcome == 4)
    
    tab = my_results_2.5h[my_results_2.5h$prof_index == index, ]
    
    withProgress(message = 'Generating results', value = 0, {
      
      n <- 3
      
      for (i in 1:n) {
        
        enriched <- enrichr(as.character(tab$gene), dbs)
        
        enrich_tab = do.call("rbind", lapply(enriched, function(x) head(x)[,c("Term", "Overlap", "P.value", "Genes")]))
        
        enrich_tab = enrich_tab[order(enrich_tab$P.value), ]
        
        output$enrich_tab = renderDataTable(enrich_tab,
                                            options = list(pageLength = 100))
        
        incProgress(1/n, detail = paste(" ", i))
        
      }
      
      
    })
    
    updateTextInput(session, "explore4", value = "")
    updateTabsetPanel(session = session, inputId = "tabs", selected = "Function/Pathways")
    
  })
  
  
  
  observeEvent(input$explore5, {
    
    index = which(codes$case == input$case & codes$outcome == 5)
    
    tab = my_results_2.5h[my_results_2.5h$prof_index == index, ]
    
    withProgress(message = 'Generating results', value = 0, {
      
      n <- 3
      
      for (i in 1:n) {
        
        enriched <- enrichr(as.character(tab$gene), dbs)
        
        enrich_tab = do.call("rbind", lapply(enriched, function(x) head(x)[,c("Term", "Overlap", "P.value", "Genes")]))
        
        enrich_tab = enrich_tab[order(enrich_tab$P.value), ]
        
        output$enrich_tab = renderDataTable(enrich_tab,
                                            options = list(pageLength = 100))
        
        incProgress(1/n, detail = paste(" ", i))
        
      }
      
      
    })
    
    updateTextInput(session, "explore5", value = "")
    updateTabsetPanel(session = session, inputId = "tabs", selected = "Function/Pathways")
    
  })
  
  
  
  observeEvent(input$explore6, {
    
    index = which(codes$case == input$case & codes$outcome == 6)
    
    tab = my_results_2.5h[my_results_2.5h$prof_index == index, ]
    
    withProgress(message = 'Generating results', value = 0, {
      
      n <- 3
      
      for (i in 1:n) {
        
        enriched <- enrichr(as.character(tab$gene), dbs)
        
        enrich_tab = do.call("rbind", lapply(enriched, function(x) head(x)[,c("Term", "Overlap", "P.value", "Genes")]))
        
        enrich_tab = enrich_tab[order(enrich_tab$P.value), ]
        
        output$enrich_tab = renderDataTable(enrich_tab,
                                            options = list(pageLength = 100))
        
        incProgress(1/n, detail = paste(" ", i))
        
      }
      
      
    })
    
    updateTextInput(session, "explore6", value = "")
    updateTabsetPanel(session = session, inputId = "tabs", selected = "Function/Pathways")
    
  })
  
  
  
  observeEvent(input$explore7, {
    
    index = which(codes$case == input$case & codes$outcome == 7)
    
    tab = my_results_2.5h[my_results_2.5h$prof_index == index, ]
    
    withProgress(message = 'Generating results', value = 0, {
      
      n <- 3
      
      for (i in 1:n) {
        
        enriched <- enrichr(as.character(tab$gene), dbs)
        
        enrich_tab = do.call("rbind", lapply(enriched, function(x) head(x)[,c("Term", "Overlap", "P.value", "Genes")]))
        
        enrich_tab = enrich_tab[order(enrich_tab$P.value), ]
        
        output$enrich_tab = renderDataTable(enrich_tab,
                                            options = list(pageLength = 100))
        
        incProgress(1/n, detail = paste(" ", i))
        
        
      }
      
      
    })
    
    updateTextInput(session, "explore7", value = "")
    updateTabsetPanel(session = session, inputId = "tabs", selected = "Function/Pathways")
    
  })
  
  
  
  observeEvent(input$explore8, {
    
    index = which(codes$case == input$case & codes$outcome == 8)
    
    tab = my_results_2.5h[my_results_2.5h$prof_index == index, ]
    
    withProgress(message = 'Generating results', value = 0, {
      
      n <- 3
      
      for (i in 1:n) {
        
        enriched <- enrichr(as.character(tab$gene), dbs)
        
        enrich_tab = do.call("rbind", lapply(enriched, function(x) head(x)[,c("Term", "Overlap", "P.value", "Genes")]))
        
        enrich_tab = enrich_tab[order(enrich_tab$P.value), ]
        
        output$enrich_tab = renderDataTable(enrich_tab,
                                            options = list(pageLength = 100))
        
        incProgress(1/n, detail = paste(" ", i))
        
      }
      
    })
    
    updateTextInput(session, "explore8", value = "")
    updateTabsetPanel(session = session, inputId = "tabs", selected = "Function/Pathways")
    
  })
  
  
  
  observeEvent(input$explore9, {
    
    index = which(codes$case == input$case & codes$outcome == 9)
    
    tab = my_results_2.5h[my_results_2.5h$prof_index == index, ]
    
    withProgress(message = 'Generating results', value = 0, {
      
      n <- 3
      
      for (i in 1:n) {
        
        enriched <- enrichr(as.character(tab$gene), dbs)
        
        enrich_tab = do.call("rbind", lapply(enriched, function(x) head(x)[,c("Term", "Overlap", "P.value", "Genes")]))
        
        enrich_tab = enrich_tab[order(enrich_tab$P.value), ]
        
        output$enrich_tab = renderDataTable(enrich_tab,
                                            options = list(pageLength = 100))
        
        incProgress(1/n, detail = paste(" ", i))
        
      }
      
    })
    
    updateTextInput(session, "explore9", value = "")
    updateTabsetPanel(session = session, inputId = "tabs", selected = "Function/Pathways")
    
  })
    
  n_sample = 2
  
  design = factor(c(rep("0", n_sample), rep("X", n_sample), rep("Y", n_sample), rep("Y + X", n_sample)))
  
  output$gene_plot = renderPlot({
    
    incProgress(1/n, detail = paste(" ", i))
    
    genetab = my_results_2.5h[which(my_results_2.5h$gene == input$gene)[1], 3:(3+length(design)-1)]
    
    colnames(genetab) = design

    boxplot(as.numeric(genetab) ~ design, ylab = 'log2(expression)', col = 'blue', medcol = "yellow")

    })

  output$francois <- renderDataTable(as.data.frame(GEOmatrix))

}

