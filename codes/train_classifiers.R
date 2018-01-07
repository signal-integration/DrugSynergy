train_classifiers = function(training_set){
  
  library(caret)
  library(randomForest)
  
  lda_model = train(training_set[[1]][,-1], 
                    as.factor(training_set[[1]][,1]), 
                    method="lda", preProcess = "pca")
  
  knn_model = train(training_set[[1]][,-1], 
                    as.factor(training_set[[1]][,1]),
                    method = "knn",
                    preProcess = "pca")
  
  rf_model = randomForest(training_set[[1]][,-1],
                          as.factor(training_set[[1]][,1]))
  
  return(list(lda_model, knn_model, rf_model))
  
  }