compute_pairwise_deltas  = function(profile_means){
  
  mean_e_0 = profile_means['0']
  mean_e_X = profile_means['X']
  mean_e_Y = profile_means['Y']
  mean_e_XY = profile_means['Y+X']
  mean_additive_level = profile_means['additivity']
  
  pairwise_deltas = c(mean_additive_level - mean_e_0,
                      mean_e_X - mean_e_0,
                      mean_e_Y - mean_e_0,
                      mean_e_XY - mean_e_0,
                      mean_e_X - mean_additive_level,
                      mean_e_Y - mean_additive_level,
                      mean_e_XY - mean_additive_level,
                      mean_e_Y - mean_e_X,
                      mean_e_XY - mean_e_X,
                      mean_e_XY - mean_e_Y)
  
  names_pairwise_deltas = c("ADD.0", "X.0", "Y.0", "X+Y.0", "X.ADD", "Y.ADD",
                            "X+Y.ADD", "Y.X", "X+Y.X", "X+Y.Y")
  
  names(pairwise_deltas) = paste("Delta", names_pairwise_deltas, sep="_")
  
  return(pairwise_deltas)
  }