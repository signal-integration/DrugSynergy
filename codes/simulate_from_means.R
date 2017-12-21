simulate_from_means = function(profile_means, replicates, noise, exp_min, exp_max){
  
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
  
  simulated_profile = c(
    
    rnorm(replicates, e_0, sd = noise), 
    rnorm(replicates, e_X, sd = noise),
    rnorm(replicates, e_Y, sd = noise), 
    rnorm(replicates, e_XY, sd = noise)
  
    )
  
  #keep simulated_profileession values within simulated_profileession range
  simulated_profile[simulated_profile < exp_min] = exp_min
  simulated_profile[simulated_profile > exp_max] = exp_max
  
  names(simulated_profile) = c(rep("0", replicates),
                               rep("X", replicates),
                               rep("Y", replicates),
                               rep("Y+X", replicates))
  return(simulated_profile)
  
}