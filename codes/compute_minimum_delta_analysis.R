compute_minimum_delta_analysis = function(Mu, min_delta){
  
  #tests if there is at least one pair of conditions with 
  #a difference larger than min_delta (in absolute value)
  #inputs: Mu = vector of means in e0, eX, eY, eX+Y, min_delta
  #output: logical variable which answers the test

  #expression corresponding to additivity
  #e0 + (eX - e0) + (eY - e0)
  additive = Mu[1] + Mu[2] - Mu[1] + Mu[3] - Mu[1]
  
  all_deltas = c(additive - Mu[1], 
                 Mu[2] - Mu[1], 
                 Mu[3] - Mu[1],
                 Mu[4] - Mu[1], 
                 Mu[2] - additive, 
                 Mu[3] - additive,
                 Mu[4] - additive, 
                 Mu[3] - Mu[2],
                 Mu[4] - Mu[2], 
                 Mu[4] - Mu[3])
  
  test_all_deltas = abs( all_deltas[all_deltas!=0]) > min_delta
  
  return(sum( test_all_deltas) >= 1)
}