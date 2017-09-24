filter_data_on_deltas = function(my_data, design, min_delta = 0.5){
  
  source("compute_minimum_delta_analysis.R")
  
  expression_data_means = t(apply(my_data[,-(1:2)], 1, function(x) tapply(x, design, mean)))
  
  minimum_delta_filter = apply(expression_data_means, 1, function(x) compute_minimum_delta_analysis(x, min_delta))
  
  my_data_filtered = my_data[which(minimum_delta_filter), ]
  
  return(my_data_filtered)
  
}

