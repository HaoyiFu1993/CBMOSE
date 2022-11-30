# update smoothing parameters using Gibbs sampling 
get_tausq <- function(nu_tau,tausq,A_tau,m,beta){
  a_tau <- rinvgamma(1,(nu_tau+1)/2,nu_tau/tausq + 1/A_tau^2)
  rinvgamma(1,(nu_tau+m)/2,t(beta)%*%beta/2 + nu_tau/a_tau)
}
