#this function generates qualitative plots for each profile

generate_all_profiles = function(){
  
  source("compute_profile_means.R")
  source("compute_minimum_delta.R")
  source("setPowerPointStyle.R")
  load("constraints_vector")
  
  #read profile codes
  PROFCODES = read.table("profile_codes_v2.txt", header = TRUE, sep = "\t")
  
  ntimes = 10 #n. of simulations
  exp_min = 2 #min range of expression value
  exp_max = 16 #max range of expression value
  min_delta = 0.5 #signal (minimum difference in expression between any two comparisons)
  mean_expression_vector = vector(length = 4)
  names(mean_expression_vector) = c("0", "X", "Y", "X+Y")
  
  #runs through all profiles
  for (k in 1:123) {
    
    prof_index = k
    
    mean_expression_vector = compute_profile_means(PROFCODES, prof_index, ntimes, exp_min, exp_max,
                                                   constraints_vector, min_delta)[, 1:4]
    
    additive_level = mean_expression_vector[1, 1] + 
      (mean_expression_vector[1, 2] - mean_expression_vector[1, 1]) +    
      (mean_expression_vector[1, 3] - mean_expression_vector[1, 1])
    
    #fix boundaries
    if (additive_level < exp_min) additive_level = 0.1
    
    if (additive_level > exp_max) additive_level = max(mean_expression_vector[1, ]) + 1.5
    
    setPowerPointStyle()
    barplot(mean_expression_vector[1, ], ylab = '', yaxt = 'n', ann = FALSE,
            col = 'black', main = paste('Profile', prof_index), ylim = c(0, max(mean_expression_vector[1, ]) + 1.6))
    
    abline(h = additive_level, col = "magenta", lwd = 3)
    
  }
  
}