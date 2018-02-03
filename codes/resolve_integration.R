resolve_integration = function(my_data, design, PROFCODES, model = rf_model){
  
  source("extract_limma_features_v1.R")
  
  stat_features = extract_limma_features_v1(my_data[,-(1:2)], design)
  
  predicted_classes = predict(rf_model, newdata = stat_features, type = 'prob')
  
  max_match = t(apply(predicted_classes, 1, function(x) c(max(x), which.max(x))))
  
  optimal_match = data.frame(my_data, max_match, PROFCODES[max_match[,2], c('case', 'outcome', 'type')])
  
  names(optimal_match)[ncol(optimal_match)-4] = "score"
  names(optimal_match)[ncol(optimal_match)-3] = "prof_index"
  names(optimal_match)[ncol(optimal_match)-2] = "case"
  names(optimal_match)[ncol(optimal_match)-1] = "outcome"
  names(optimal_match)[ncol(optimal_match)] = "type"
  
  #set to constant profiles genes with p.adj > 0.05
  adjusted_pvals = p.adjust(stat_features[,1], method = 'BH')
  
  optimal_match[adjusted_pvals>0.05, "prof_index"] = 2
  optimal_match[adjusted_pvals>0.05, "case"] = 1
  optimal_match[adjusted_pvals>0.05, "outcome"] = 2
  optimal_match[adjusted_pvals>0.05, "type"] = 'A'
  
  #compute minimum relevant fold-change
  min_fc = vector()
  for (s in 1:nrow(optimal_match)){
    
    min_fc[s] = compute_minimum_delta(as.numeric(stat_features[s,2:5]),
                                      PROFCODES,
                                      optimal_match$prof_index[s])
  }
  
  optimal_match$min_fc = min_fc
  
  optimal_match$adjusted_pvals = adjusted_pvals
  
  return(optimal_match)
  }