simulate_beta_from_means = function(profile_means, replicates, sigma, exp_min = 2, exp_max = 16){
  
  #REQUIRES
  #profile_means: a vector of profile means (e0, eX, eY, eX+Y) as computed by "compute_profile_means" 
  #replicates = how many replicates should be simulated for each condition
  #noise: sigma in the normal distribution
  #exp_min = minimum value of the expression range
  #exp_max = maximum value of the expression range
  
  #OUTPUT
  #a simulated profile with replicated instances of e0, eX, eY, eX+Y
  
  e_0 = profile_means[1]
  e_X = profile_means[2]
  e_Y = profile_means[3]
  e_XY = profile_means[4]
  
  signal = profile_means[5]
  
  #for constant profile
  if (signal == Inf){
    
    simulated_profile = c(
      
      rnorm(replicates, e_0, sd = sigma), 
      rnorm(replicates, e_X, sd = sigma),
      rnorm(replicates, e_Y, sd = sigma), 
      rnorm(replicates, e_XY, sd = sigma)
      
    )
  }
  else{
    
    
    #avoid excessive noise

    simulated_profile = c(
      
      simulate_beta_on_interval(replicates, e_0, sigma),
      simulate_beta_on_interval(replicates, e_X, sigma),
      simulate_beta_on_interval(replicates, e_Y, sigma),
      simulate_beta_on_interval(replicates, e_XY, sigma)

    )
  }
  

  names(simulated_profile) = c(rep("CTRL", replicates),
                               rep("X", replicates),
                               rep("Y", replicates),
                               rep("Y+X", replicates))
  return(simulated_profile)
  
}