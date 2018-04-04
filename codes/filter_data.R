filter_data = function(combinatorial_data){
  
  if (max(combinatorial_data[,-(1:2)])>25) {
    
    combinatorial_data[,-(1:2)] = log2(combinatorial_data[,-(1:2)])
  
    }
  
  #remove rows without gene symbol
  combinatorial_data = combinatorial_data[!is.na(combinatorial_data$genes), ]
  
  #compute the coefficient of variation for each row
  combinatorial_data$cof = apply(combinatorial_data[,-(1:2)], 1, function(x) sd(x)/abs(mean(x)))
  
  #select gene symbols with largest cof
  #combinatorial_data = aggregate(cof ~ genes + ., combinatorial_data, function(x) x[,which.max(x)])
  
  combinatorial_data = data.frame(combinatorial_data %>% group_by(genes) %>% slice(which.max(cof)))
  
  #filter on cof
  cof_cutoff = summary(combinatorial_data$cof)[3]
  
  cof_filter = which(combinatorial_data$cof >= cof_cutoff)
  
  filtered_data = combinatorial_data[cof_filter, ]
  
  filtered_data$cof = NULL

  return(filtered_data)
}
