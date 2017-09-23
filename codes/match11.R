match11 = function(expr) {
  
  source("mypvals.R")
  #   this function takes as input expr, a vector of expression    
  #   values in 0,X,Y,Y+X. Same n. of replicates per condition is assumed.
  #   It returns:
  #   - the Bliss index
  #   - the average expression values for each condition
  #   - the average deltas in expression for all pairwise comparisons   
  #   - all the one-tailed and two-tailed t-test p-values for all pairwise comparisons 
  #   - all the one-tailed and two-tailed Wilcoxon p-values for all pairwise comparisons 
  #   The total number of features is xxx
  #   These features will be used to train a classifier with synthetic profiles and machine learning.
  
  #   to run the function on a synthetic profile (emergent positive synergy):
  #   source("match11.R")
#     N=4
#     sd=0.5
#     mean.values=c(2,2,2,4)
#     expr = c(rnorm(N,mean.values[1],sd = sd),rnorm(N,mean.values[2],sd = sd),
#             rnorm(N,mean.values[3],sd = sd),rnorm(N,mean.values[4],sd = sd))
  #   profile.features=match11(expr)
  
  
  #number of replicatates (assuming same n. of donors/condition)
  N = length(expr)/4


  #this is a factor whose levels correspond to expression in: 0, X, Y, X+Y  
  myfactor = factor(c(rep("0",N),rep("X",N),rep("Y",N),
                      rep("Y+X",N),rep("ADD_id",N) ))
  
  #compute absolute level of additivity for each replicate 

    ADD_id = expr[which(myfactor == "0")] +
    (expr[which(myfactor == "Y")] - expr[which(myfactor == "0")]) +
    (expr[which(myfactor == "X")] - expr[which(myfactor == "0")])

  #extend factor F and expression vector to include absolute level of additivity  
  expr.b = c(expr,ADD_id)

  #calculate vector of means
  means=tapply(expr.b,myfactor,mean)
  
  #bliss index
  Bliss = means["Y+X"]-means["ADD_id"]
  
  Mu=means[-2]
  
  
  #initialize vector with all possible pairwise fold-changes
  V1=vector(mode="numeric",length=10)
  V1 = c(means["ADD_id"] - Mu[1], Mu[2] - Mu[1],Mu[3] - Mu[1], Mu[4] - Mu[1],
         Mu[2] - means["ADD_id"], Mu[3] - means["ADD_id"], Mu[4] - means["ADD_id"], Mu[3] - Mu[2],
         Mu[4] - Mu[2], Mu[4] - Mu[3])
  
  names1=c("ADD.0","X.0","Y.0","X+Y.0","X.ADD","Y.ADD","X+Y.ADD","Y.X","X+Y.X","X+Y.Y")
  names(V1)=paste("Delta",names1,sep="_")
  
  
  #all t-tests p-values
  EQP=mypvals(expr.b,myfactor,"t.test","two.sided")
  UPP=mypvals(expr.b,myfactor,"t.test","greater")
  DOWNP=mypvals(expr.b,myfactor,"t.test","less")
  
  #all wicoxon p-values
  EQPw=mypvals(expr.b,myfactor,"wilcox.test","two.sided")
  UPPw=mypvals(expr.b,myfactor,"wilcox.test","greater")
  DOWNPw=mypvals(expr.b,myfactor,"wilcox.test","less")
  

  Q=c(Bliss,Mu,V1,1-EQP,DOWNP,UPP,1-EQPw,DOWNPw,UPPw)
  
  names(Q)[1]="Bliss"
  
  names(Q)=make.names(names(Q))
  
  return(Q)
}





