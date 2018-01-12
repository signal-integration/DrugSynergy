test_accuracy = function(classifiers, test_set){
  
  #test LDA
  lda_test = diag(table(test_set[[1]][,1], predict(classifiers[[1]],
                                                   test_set[[1]][,-1])))
  lda_test = data.frame(as.numeric(lda_test), "LDA")
  names(lda_test) = c("accuracy", "model")
  
  
  #test KNN
  knn_test = diag(table(test_set[[1]][,1], predict(classifiers[[2]],
                                                   test_set[[1]][,-1])))
  
  knn_test = data.frame(as.numeric(knn_test), "KNN")
  names(knn_test) = c("accuracy", "model")
  
  
  #test RF
  rf_test = diag(table(test_set[[1]][,1], predict(classifiers[[3]],
                                                  test_set[[1]][,-1])))
  
  rf_test = data.frame(as.numeric(rf_test), "RF")
  names(rf_test) = c("accuracy", "model")
  
  
  test = rbind(knn_test, lda_test, rf_test)
  
  density2 <- ggplot(data = test, aes(x = accuracy, fill = model)) + 
    geom_density(stat = "density", alpha=I(0.7)) + xlim(0, 100)

  density2 = density2 + xlab("% accuracy on test set") + theme_bw() + theme(axis.text.y = element_blank(),axis.ticks = element_blank(), 
                             panel.border = element_blank(), 
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank()) + scale_fill_manual( values = c("red","blue", 'green'))

  density2
#  density2 = density2 + facet_grid(model ~ .) + xlab("% accuracy per profile") +  theme_bw() + theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank()) 
#  density2 = density2 + theme(
#    axis.text.y = element_blank(),
#    axis.ticks = element_blank())
#  density2
}
