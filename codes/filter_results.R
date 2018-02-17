filter_results = function(results, adjusted_pval = 0.05){
  
  #filter by adj_pval
  filter_deg = results[[1]]$adjusted_pvals < adjusted_pval
  
  results_f = results[[1]][filter_deg, ]
  features_f = results[[2]][filter_deg, ]
  
  #filter by median fc
  median_fc = median(results_f$min_fc)
  
  results_ff = results_f[results_f$min_fc >= median_fc, ]
  feautures_ff = features_f[results_f$min_fc >= median_fc, ]
  
  return(list(results_ff, feautures_ff))
}