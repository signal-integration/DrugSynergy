library(limma)
library(Biobase)

extract_limma_features = function(integration_profile, design){
  
  replicates = ncol(integration_profile)/4
  
  names(integration_profile) = design
  
  #compute expected additivity line
  e_0 = integration_profile[,which(design == "CTRL")]
  e_X = integration_profile[,which(design == "X")]
  e_Y = integration_profile[,which(design == "Y")]
  additivity = e_0 + (e_X - e_0) + (e_Y - e_0)
  
  #extend design vector and data matrix to include additivity line
  design_extended = as.factor(c(colnames(integration_profile), rep('add', replicates)))
  
  my_data = data.frame(integration_profile, additivity)  
  names(my_data) = make.names(as.character(design_extended), unique = T)
  
  #convert to eset
  my_eset = new("ExpressionSet", exprs = as.matrix(my_data))
  
  design <- model.matrix(~ 0 + design_extended)

  colnames(design) <- make.names(c("CTRL", "X", "Y", "Y+X", "add"))
  
  fit <- lmFit(my_eset, design)

  contrast.matrix <- makeContrasts(add - CTRL, 
                                   X - CTRL, 
                                   Y - CTRL, 
                                   Y.X - CTRL,
                                   X - add,
                                   Y - add,
                                   Y.X - add,
                                   Y - X,
                                   Y.X - X,
                                   Y.X - Y,
                                   levels = design)
  
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  #features to add
  dt = decideTests(fit2)
  stat_features = data.frame(p.adjust(fit2$F.p.value, method="fdr"), 
                             dt[,1:10], 
                             fit2$coefficients, 
                             fit2$sigma)
  return(stat_features)
}


#data_file = "tnf_ifn_1_v3"
#my_data = read.csv(data_file)
#additivity = my_data[,2+(1:3)] + (my_data[,2+(4:6)]- my_data[,2+(1:3)]) + (my_data[,2+(7:9)]- my_data[,2+(1:3)])
#my_data = data.frame(my_data, additivity)  
#my_eset = new("ExpressionSet", exprs = as.matrix(my_data[,-(1:2)]))
#featureNames(my_eset) = my_data[,1]

#design <- model.matrix(~ 0+factor(c(1,1,1,2,2,2,3,3,3, 4, 4, 4, 5, 5, 5)))
#fit <- lmFit(my_eset, design)

# colnames(design) <- make.names(c("0", "X", "Y", "Y+X", "add"))
# 
# contrast.matrix <- makeContrasts(add - X0, 
#                                  X - X0, 
#                                  Y - X0, 
#                                  Y.X - X0,
#                                  X - add,
#                                  Y - add,
#                                  Y.X - add,
#                                  Y - X,
#                                  Y.X - X,
#                                  Y.X - Y,
#                                  levels = design)
# 
# fit2 <- contrasts.fit(fit, contrast.matrix)
# fit2 <- eBayes(fit2)
# 
# #features to add
# dt = decideTests(fit2)
# stat_features = data.frame(dt[,1:10], 
#                 fit2$coefficients, 
#                 fit2$sigma,
#                 fit2$p.value)

#sum(p.adjust(fit2$F.p.value, method="fdr")<0.05)
