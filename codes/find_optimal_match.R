find_optimal_match = function(my_data_filtered){
  
  library(caret)
  library(randomForest)
  source("extract_stat_features.R")
  source("setPowerPointStyle.R")
  setPowerPointStyle()

  load("rf_model")

  #this small perturbations fixes the statistical errors related to null group variances
  noise_mat = replicate(12, rnorm(nrow(my_data_filtered[, -(1:2)]), sd = 0.0001))
  
  profile_features = apply(my_data_filtered[, -(1:2)] + noise_mat, 1, extract_stat_features)
  
  predicted_classes = predict(rf_model, newdata = data.frame(t(profile_features)), type = 'prob')

  max_match = t(apply(predicted_classes, 1, function(x) c(max(x), which.max(x))) )
  
  optimal_match = data.frame(my_data_filtered, max_match)
  
  names(optimal_match)[ncol(optimal_match)-1] = "score"
  names(optimal_match)[ncol(optimal_match)] = "prof_index"
  
  return(optimal_match)
}
