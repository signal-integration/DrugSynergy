extract_stat_features = function(integration_profile) {
  
  source("mypvals.R")
  
  #REQUIRES
  #integration_profile a vector of expression values with n replicates of 0, X, Y, X+Y
  #(order matters)
  
  #OUTPUT
  #   - the Bliss index
  #   - the average integration_profileession values for each condition
  #   - the average deltas in integration_profileession for all pairwise comparisons   
  #   - all the one-tailed and two-tailed t-test p-values for all pairwise comparisons 
  #   - all the one-tailed and two-tailed Wilcoxon p-values for all pairwise comparisons 
  #   The total number of features is xxx
  #   These features will be used to train a classifier with synthetic profiles and machine learning.
  
  
  #number of replicates (assuming same n. of donors/condition)
  replicates = length(integration_profile)/4

  design_factor = as.factor(c(names(integration_profile), rep('additivity', replicates)))
  
  e_0 = integration_profile[which(design_factor == "0")]
  e_X = integration_profile[which(design_factor == "X")]
  e_Y = integration_profile[which(design_factor == "Y")]

  additivity = as.vector(e_0 + (e_X - e_0) + (e_Y - e_0))

  integration_profile_with_additivity = c(integration_profile, additivity)
  names(integration_profile_with_additivity) = design_factor

  #calculate vector of means
  profile_means = tapply(integration_profile_with_additivity, design_factor, mean)
  
  #mean bliss index
  bliss = as.numeric(profile_means["Y+X"] - profile_means["additivity"])
  
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
  
  
  #all t-tests p-values
  EQP = mypvals(integration_profile_with_additivity, design_factor, "t.test", "two.sided")
  UPP = mypvals(integration_profile_with_additivity, design_factor, "t.test", "greater")
  DOWNP = mypvals(integration_profile_with_additivity, design_factor, "t.test", "less")
  
  #all wicoxon p-values
  EQPw = mypvals(integration_profile_with_additivity, design_factor, "wilcox.test", "two.sided")
  UPPw = mypvals(integration_profile_with_additivity, design_factor, "wilcox.test", "greater")
  DOWNPw = mypvals(integration_profile_with_additivity, design_factor, "wilcox.test", "less")
  
  statistical_features = c(bliss, profile_means[-2], pairwise_deltas, 1-EQP, DOWNP, UPP, 1-EQPw, DOWNPw, UPPw)
  
  names(statistical_features)[1] = "Bliss"
  
  names(statistical_features) = make.names(names(statistical_features))
  
  return(statistical_features)
  
  }





