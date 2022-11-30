# obtain a list of initial values for parameters delta, beta_star, sigmasq, tausq, pi, ind and z
get_init <- function(N,ncomp,k,m,V){
  
  # logistic parameters
  ndelta=ncol(V)
  delta <- matrix(0,nrow = ndelta,ncol = ncomp)
  if (ncomp==1){
    delta[,1]=rep(0,ndelta)
  } else{
    delta[,1:(ncomp-1)] <- runif((ndelta)*(ncomp-1),0,0.1)
  }
  
  # intercept, slope and basis function coefficients
  beta_star <- array(rep(0,(m+2)*ncomp*k),dim = c(m+2,ncomp,k))
  # error variance
  sigmasq <- matrix(1,nrow = ncomp,ncol = k)
  # smoothing parameters
  tausq <- matrix(1,nrow = ncomp,ncol = k)
  
  # mixing weights
  exb = exp(V %*% delta)
  pi=exb/apply(exb,1,sum)
  ind=rep(NA,N)
  for (i in 1:N){
    ind[i] <- sample(1:ncomp,1,replace = T,prob = pi[i,])
  }
  # table(ind)
  # table(ind,ind_true_all[,1])
  
  #binary indicators
  z=matrix(0,N,ncomp)
  for (j in 1:ncomp) z[ind==j,j]=1
  
  # return a list of parameters
  return(list(delta=delta,beta_star=beta_star,sigmasq=sigmasq,tausq=tausq,pi=pi,ind=ind,z=z))
}
