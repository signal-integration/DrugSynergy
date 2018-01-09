extract_stat_features = function(integration_profile, design) {
  
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
  
  source("compute_pairwise_deltas.R")
  source("get_p_values.R")
  
  #number of replicates (assuming same n. of donors/condition)
  replicates = length(integration_profile)/4

  integration_profile = as.numeric(integration_profile)
  names(integration_profile) = design
  
  design_factor = as.factor(c(names(integration_profile), rep('additivity', replicates)))
  
  e_0 = integration_profile[which(design_factor == "0")]
  e_X = integration_profile[which(design_factor == "X")]
  e_Y = integration_profile[which(design_factor == "Y")]

  additivity = as.vector(e_0 + (e_X - e_0) + (e_Y - e_0))

  integration_profile_with_additivity = c(integration_profile, additivity)
  names(integration_profile_with_additivity) = design_factor

  #calculate vector of means
  profile_means = tapply(integration_profile_with_additivity, design_factor, mean)
  profile_sds = tapply(integration_profile_with_additivity, design_factor, sd)
  names(profile_sds) = c("X0_sd", "add_sd", "X_sd", "Y_sd", "Y.X_sd")
  #mean bliss index
  #bliss = as.numeric(profile_means["Y+X"] - profile_means["additivity"])
  
  #getting all pairwise deltas 
  pairwise_deltas = compute_pairwise_deltas(profile_means)
  
  #all t-tests p-values
  EQP = get_p_values(integration_profile_with_additivity, design_factor, "t.test", "two.sided")
  UPP = get_p_values(integration_profile_with_additivity, design_factor, "t.test", "greater")
  DOWNP = get_p_values(integration_profile_with_additivity, design_factor, "t.test", "less")
  
  #all wicoxon p-values
  #EQPw = get_p_values(integration_profile_with_additivity, design_factor, "wilcox.test", "two.sided")
  #UPPw = get_p_values(integration_profile_with_additivity, design_factor, "wilcox.test", "greater")
  #DOWNPw = get_p_values(integration_profile_with_additivity, design_factor, "wilcox.test", "less")
  
  #statistical_features = c(bliss, profile_means[-2], pairwise_deltas, EQP, DOWNP, UPP)#, 1-EQPw, DOWNPw, UPPw)
  
  #names(statistical_features)[1] = "Bliss"
  
  statistical_features = c(profile_means, profile_sds, pairwise_deltas, EQP, DOWNP, UPP)
  
  names(statistical_features) = make.names(names(statistical_features))
  names(statistical_features)[2] = 'add_level'
  return(statistical_features)
  
  }





