simulate_from_means = function(profile_means, replicates, min_noise, max_noise) {
  
  #this function simulates a noisy profile 
  #INPUTS
  #profile_means: a vector of mean expression values as computed with the function compute_profile_means)
  #replicates: number of replicates for each condition 
  #min_noise: minumum level of noise
  #max_noise: maximum level of noise 

  Mu = profile_means
  
  ADD = Mu[1] + Mu[2] - Mu[1] + Mu[3] - Mu[1]
  
  V1 = c(
    ADD - Mu[1],
    Mu[2] - Mu[1],
    Mu[3] - Mu[1],
    Mu[4] - Mu[1],
    Mu[2] - ADD,
    Mu[3] - ADD,
    Mu[4] - ADD,
    Mu[3] - Mu[2],
    Mu[4] - Mu[2],
    Mu[4] - Mu[3]
  )
  
  
  #if constant profile, minimum shift is zero: so noise is in range 0.2-1
  if (var(Mu) == 0) {
    
    sd = runif(1, min = min_noise, max = max_noise)
    
  } else {
    #if not constant, noise is 20%-50% of minimum shift
    
    D = min(abs(V1[ind]))
    
    sd = runif(1, min = min_noise, max = max_noise) * D
    
  }
  
  
  #here simulate the distributions
  expr = c(
    rnorm(replicates, Mu[1], sd = sd),
    rnorm(replicates, Mu[2], sd = sd),
    
    rnorm(replicates, Mu[3], sd = sd),
    rnorm(replicates, Mu[4], sd = sd)
  )
  
  return(expr)
  
}