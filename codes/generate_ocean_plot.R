generate_ocean_plot = function(deg){
  
  library(alluvial)
  
  #remove additive
  deg = deg[deg$type != 'A', ]
  
  deg$type = droplevels(deg$type)
  
  alluvial_cases = melt(table(deg$type, deg$case))
  
  col = as.character(alluvial_cases$Var.1)
  col[col == 'N'] = 'red'
  col[col == 'P'] = 'blue'
  
  alluvial(alluvial_cases[,1:2], freq = alluvial_cases$value, 
           col = col, cex = 0.01, alpha=0.8, 
           hide = alluvial_cases$value < 10)
  }