mypvals = function(expr.b, myfactor, test.type, alternative) {
  
  if (test.type == "t.test") {
    test = tryCatch(
      pairwise.t.test(
        expr.b,
        myfactor,
        alternative = alternative,
        pool.sd = TRUE
      ),
      error = function(e)
        NA
    )
  }
  
  else if (test.type == "wilcox.test") {
    test = tryCatch(
      pairwise.wilcox.test(
        expr.b,
        myfactor,
        alternative = alternative,
        pool.sd = TRUE
      ),
      error = function(e)
        NA
    )
    
  }
  
  #initialize vector of p-values
  testP = rep(NA, 10)
  
  #pull out all p-values
  if (!is.na(test[1]))
    
    testP = c(test[[3]][1:4, 1], test[[3]][2:4, 2], test[[3]][3:4, 3], test[[3]][4, 4])
  
  names1 = c("ADD.0",
             "X.0",
             "Y.0",
             "X+Y.0",
             "X.ADD",
             "Y.ADD",
             "X+Y.ADD",
             "Y.X",
             "X+Y.X",
             "X+Y.Y")
  
  names(testP) = paste(test.type, alternative, names1, sep = "_")
  
  return(testP)
  
}
