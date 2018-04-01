library(limma)
library(Biobase)

extract_limma_features_v1 = function(integration_profile, design){
  
  names(integration_profile) = make.names(as.character(design), unique = T)
  
  my_eset = new("ExpressionSet", exprs = as.matrix(integration_profile))
  
  design_matrix <- model.matrix(~ 0 + design)
  
  fit <- lmFit(my_eset, design_matrix)
  
  contrast_matrix <- makeContrasts(designX - designCTRL,
                                   designY - designCTRL,
                                   designY - designX,
                                   designYX - designCTRL,
                                   designYX - designX,
                                   designYX - designY,
                                   levels = design_matrix)
  
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)

  dt = classifyTestsP(fit2)
  
  stat_features = data.frame(fit2$F.p.value, fit$coefficients, fit2$coefficients, fit2$p.value, dt[,1:6], fit2$sigma)
  
  return(stat_features)
}





#compute expected additivity line
#  e_0 = integration_profile[,which(design == "CTRL")]
#  e_X = integration_profile[,which(design == "X")]
#  e_Y = integration_profile[,which(design == "Y")]
#  additivity = e_0 + (e_X - e_0) + (e_Y - e_0)
#  names(additivity) = rep('add', replicates)
#extend design vector and data matrix to include additivity line
#  design_extended = as.factor(c(colnames(integration_profile), rep('add', replicates)))
#  levels(design_extended) = c("CTRL", "X", "Y", "combo", "add")

# integration_profile_extended = data.frame(integration_profile, additivity)  
# names(integration_profile_extended) = make.names(as.character(design_extended), unique = T)
# 
# #convert to eset
# my_eset = new("ExpressionSet", exprs = as.matrix(integration_profile_extended))
# 
# design_matrix <- model.matrix(~ 0 + design_extended)
# 
# #colnames(design_matrix) <- make.names(c("CTRL", "X", "Y", "Y+X", "add"))
# 
# fit <- lmFit(my_eset, design_matrix)
# 
# 
# contrast_matrix <- makeContrasts(design_extendedadd - design_extendedCTRL, 
#                                  design_extendedX - design_extendedCTRL, 
#                                  design_extendedY - design_extendedCTRL, 
#                                  design_extendedcombo - design_extendedCTRL,
#                                  design_extendedX - design_extendedadd,
#                                  design_extendedY - design_extendedadd,
#                                  design_extendedcombo - design_extendedadd,
#                                  design_extendedY - design_extendedX,
#                                  design_extendedcombo - design_extendedX,
#                                  design_extendedcombo - design_extendedY,
#                                  levels = design_matrix)
# 
# fit2 <- contrasts.fit(fit, contrast_matrix)
# fit2 <- eBayes(fit2)
# 
# dt = classifyTestsP(fit2)
# 
# stat_features = data.frame(fit2$F.p.value, 
#                            dt[,1:10], 
#                            fit2$coefficients, 
#                            fit2$sigma)
