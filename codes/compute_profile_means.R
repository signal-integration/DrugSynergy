compute_profile_means = function(PROFCODES,
                                 prof_index,
                                 ntimes,
                                 exp_min,
                                 exp_max,
                                 constraints_vector,
                                 min_delta) {
  
  #this function computes a vector of mean values of (e0, eX, eY, eX+Y) for a profile of interest.
  
  #INPUTS
  #PROFCODES (definitions of all profiles) 
  #prof_index: index of profile of interest
  #ntimes: how many different vectors of means is computed;
  #exp_min: min of expression range 
  #exp_max: max of expression range;
  #constraints_vector: expresses the inequalities of the profiles in the format needed for linsolve
  #min_delta: minimum non-zero difference of all pairwise comparisons
  
  #OUTPUT
  #a vector of mean values of (e0, eX, eY, eX+Y) for a profile of interest
  
  source("compute_minimum_delta.R")
  
  library(limSolve)
  
  options(warn = -1)
  
  b = cbind(as.numeric(PROFCODES[prof_index, 1:10]))
  
  E = matrix(constraints_vector[b == 0, ], ncol = 4)
  
  F = cbind(rep(0, dim(E)[1]))
  
  G1 = matrix(constraints_vector[b == 1, ], ncol = 4)
  H1 = cbind(rep(0, dim(G1)[1]))
  
  G2 = -matrix(constraints_vector[b == -1, ], ncol = 4)
  H2 = cbind(rep(0, dim(G2)[1]))
  
  Gvar = rbind(diag(4), -diag(4))
  
  Hvar = cbind(c(
    exp_min,
    exp_min,
    exp_min,
    exp_min,-exp_max,
    -exp_max,
    -exp_max,
    -exp_max
  ))
  
  G = rbind(G1, G2, Gvar)
  H = rbind(H1, H2, Hvar)
  
  if (dim(E)[1] == 0 & dim(G)[1] > 0) {
    synth = xsample(G = G,
                    H = H,
                    iter = (5 * ntimes + 1))[[1]]
    
  } else {
    synth = xsample(
      E = E,
      F = F,
      G = G,
      H = H,
      iter = (5 * ntimes + 1)
    )[[1]]
    
  }
  
  synth = synth[-1, ]
  
  min_shift = apply(synth, 1, function(x) compute_minimum_delta(x, PROFCODES, prof_index))
  
  output = cbind(synth, min_shift)
  
  output = output[which(min_shift > min_delta), ]
  
  output = output[sample(1:ntimes), ]
  
  return(output)
  
}
