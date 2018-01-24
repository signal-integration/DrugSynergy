simulate_beta_on_interval = function(replicates, mu, sigma, expr_min = 2, expr_max = 16){
  
  lambda = (mu - a)*(b - mu)/sigma^2 - 1

  #parameters of beta distribution
  a1 = lambda*(mu - a)/(b - a)
  a2 = lambda*(b - mu)/(b - a)
  
  simulated_beta = a + (b - a)*rbeta(replicates, a1, a2)
  
  return(simulated_beta)
  
  }