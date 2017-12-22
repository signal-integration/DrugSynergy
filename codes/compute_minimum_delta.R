compute_minimum_delta = function(profile_means, PROFCODES, prof_index){

  #REQUIRES
  #profile_means: a vector of profile means (e0, eX, eY, eX+Y) as computed by "compute_profile_means" 
  #PROFCODES: definitions of all profiles in terms of ternary vectors
  #prof_index = index of a profile of interest (1 to 123)

  #OUTPUT
  #the minimum non-zero difference of all pairwise comparisons among (additive_level, e0, eX, eY, eX+Y) 

  e_0 = profile_means[1]
  e_X = profile_means[2]
  e_Y = profile_means[3]
  e_XY = profile_means[4]
  
  additive_level = e_0 + (e_X - e_0) + (e_Y - e_0)
  
  pairwise_deltas = c(additive_level - e_0,
                      e_X - e_0,
                      e_Y - e_0,
                      e_XY - e_0,
                      e_X - additive_level,
                      e_Y - additive_level,
                      e_XY - additive_level,
                      e_Y - e_X,
                      e_XY - e_X,
                      e_XY - e_Y)

  non_zero_pairwise_deltas_indeces = which(PROFCODES[prof_index,1:10]!=0)
  
  non_zero_pairwise_deltas = pairwise_deltas[non_zero_pairwise_deltas_indeces]
  
  min_delta = min(abs(non_zero_pairwise_deltas))
  
  return(min_delta)
  
  }