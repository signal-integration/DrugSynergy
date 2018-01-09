generate_synthetic_data = function(expr_range, ntimes, 
                                   signal_to_noise_range, PROFCODES, 
                                   samples = 4, min_delta = 0.5){
  
  source("compute_profile_means.R")
  source("simulate_from_means.R")
  source("extract_stat_features.R")
  load("constraints_vector")
  
  #n. of vars in the dataframe (n. features + class label)
  NVAR = 50 + 1
  NVAR1 = 16
  
  #design factor
  design = factor(c(rep("0", samples), rep("X", samples), rep("Y", samples), rep("Y+X", samples)))
  
  #initialize dataframe containing training set
  big_simulation = data.frame()
  big_simulated_data = data.frame()
  
  for (h in 1:length(signal_to_noise_range)) {
    
    for (k in 1:123) {
      
      temp = data.frame(matrix(ncol = NVAR, nrow = ntimes))
      temp[, 1] = k
      
      temp1 = data.frame(matrix(ncol = NVAR1, nrow = ntimes))

      simulated_means = compute_profile_means(PROFCODES, k, ntimes, 
                                              min(expr_range), max(expr_range),
                                              constraints_vector, min_delta)
      
      signal_to_noise = signal_to_noise_range[h]
      
      simulated_data = t(apply(simulated_means, 
                               1, function(x) simulate_from_means(x, samples,
                                                                  signal_to_noise,
                                                                  min(expr_range),
                                                                  max(expr_range))))
      temp1 = simulated_data
      big_simulated_data = rbind(big_simulated_data, temp1)
      

      features = t(apply(simulated_data, 1, 
                         function(x) extract_stat_features(x, design)))
      
      temp[, 2:NVAR] = features
      big_simulation = rbind(big_simulation, temp)
      
    }
    
  }
  
  names(big_simulation) = c("TCIND", colnames(features))
  
  return(list(big_simulation, big_simulated_data))
}