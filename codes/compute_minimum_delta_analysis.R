compute_minimum_delta_analysis = function(Mu,delta){
  
  ADD = Mu[1] + Mu[2] - Mu[1] + Mu[3] - Mu[1]
  
  V1 = c(ADD - Mu[1],Mu[2] - Mu[1],Mu[3] - Mu[1],Mu[4] - Mu[1],
         Mu[2] - ADD,Mu[3] - ADD,Mu[4] - ADD,Mu[3] - Mu[2],
         Mu[4] - Mu[2],Mu[4] - Mu[3])
  
  return(sum( abs( V1[V1!=0])>delta)>=1)
}