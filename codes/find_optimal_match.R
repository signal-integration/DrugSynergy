find_optimal_match = function(my_data_filtered, design, PROFCODES){
  
  library(caret)
  library(randomForest)
  #source("extract_stat_features.R")
  source("extract_limma_features.R")
  source("setPowerPointStyle.R")
  setPowerPointStyle()

  #load("rf_model")
  load("classifiers")
  rf_model = classifiers[[3]]

  #this small perturbations fixes the statistical errors related to null group variances
  #noise_mat = replicate((ncol(my_data_filtered)-2), runif(nrow(my_data_filtered[, -(1:2)]), 0.001, 0.005))
  
  #profile_features = apply(my_data_filtered[, -(1:2)] + noise_mat, 1, 
  #                         function(x) extract_stat_features(x, design))
  
  profile_features = extract_limma_features(my_data_filtered, design)
  
  predicted_classes = predict(rf_model, newdata = data.frame(profile_features), type = 'prob')

  max_match = apply(predicted_classes, 1, function(x) c(max(x), which.max(x)))
                                                         
  
  optimal_match = data.frame(my_data_filtered, max_match, PROFCODES[max_match[,2], c('case', 'outcome', 'type')])
  
  names(optimal_match)[ncol(optimal_match)-4] = "score"
  names(optimal_match)[ncol(optimal_match)-3] = "prof_index"
  names(optimal_match)[ncol(optimal_match)-2] = "case"
  names(optimal_match)[ncol(optimal_match)-1] = "outcome"
  names(optimal_match)[ncol(optimal_match)] = "type"
  
  
  return(optimal_match)
}
