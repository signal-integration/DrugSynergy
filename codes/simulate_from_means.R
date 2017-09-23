simulate_from_means=function(profile_means,prof_code,replicates,noise,expr_min,expr_max){
  
  Mu=profile_means

  ADD = Mu[1] + Mu[2] - Mu[1] + Mu[3] - Mu[1]
  
  V1 = c(ADD - Mu[1],Mu[2] - Mu[1],Mu[3] - Mu[1],Mu[4] - Mu[1],
         Mu[2] - ADD,Mu[3] - ADD,Mu[4] - ADD,Mu[3] - Mu[2],
         Mu[4] - Mu[2],Mu[4] - Mu[3])
  
  sd = noise
    
  expr = c(
    
    rnorm(replicates,Mu[1],sd = sd),rnorm(replicates,Mu[2],sd = sd),
    
    rnorm(replicates,Mu[3],sd = sd),rnorm(replicates,Mu[4],sd = sd)
  )
  
  #keep expression values within expression range
  expr[expr<expr_min]=expr_min
  expr[expr>expr_max]=expr_max
  
  return(expr)
  
}