generate_training_set = function(samples, ntimes, expr_range = c(2,16), min_delta = 0.5, noise_range = c(0.2, 0.4, 0.6, 0.8, 1.0)){

  source("compute_profile_means.R")
  source("simulate_from_means.R")
  source("extract_stat_features.R")
  load("constraints_vector")
  
  noise_range=c(0.2, 0.4, 0.6, 0.8, 1.0)
  
  #n. of vars in the dataframe (75 features + class label)
  NVAR=75 + 1
  
  #initialize dataframe containing training set
  big_simulation = data.frame()
  
  for (h in 1:length(noise_range)){
    
    for (k in 1:123){
      
      temp = data.frame(matrix(ncol = NVAR, nrow = ntimes))
      temp[,1] = k
      simulated_means = compute_profile_means(PROFCODES, k, ntimes, min(expr_range),
                                              max(expr_range), constraints_vector, min_delta)[,-5]
      
      simulated_data = t(apply(simulated_means, 1, function(x) simulate_from_means(x, k, samples,noise_range[h],
                                                                                   min(expr_range), max(expr_range))))
      
      features = t(apply(simulated_data, 1, extract_stat_features))
      
      temp[,2:NVAR] = features
      
      big_simulation = rbind(big_simulation, temp)
      
    }
    
  }
  
  names(big_simulation) = c("TCIND", colnames(features))
  
  return(big_simulation)
  
  }