find_optimal_match = function(my_data_filtered){
  
  library(caret)
  library(randomForest)
  source("match11.R")
  source("setPowerPointStyle.R")
  setPowerPointStyle()

  load("rf_model")

  profile_features = apply(my_data_filtered[, -(1:2)] + replicate(12, rnorm(nrow(my_data_filtered[, -(1:2)]), sd = 0.0001)), 1, match11)
  
  predicted_classes=predict(rf_model, newdata = data.frame(t(profile_features)), type = 'prob')

  max_match=t(apply(predicted_classes, 1, function(x) c(max(x), which.max(x))) )
  
  optimal_match = data.frame(my_data_filtered, max_match)
  
  names(optimal_match)[ncol(optimal_match)-1] = "score"
  names(optimal_match)[ncol(optimal_match)] = "prof_index"
  
  return(optimal_match)
}
