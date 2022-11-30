# update variance of random intercepts using Gibbs sampling 
get_zetasq <- function(nu_zeta,zetasq,A_zeta,m,delta,N){
  a_zeta <- rinvgamma(1,(nu_zeta+1)/2,nu_zeta/zetasq + 1/A_zeta^2)
  zetasq <- rinvgamma(1,(nu_zeta+N)/2,t(delta)%*%delta/2 + nu_zeta/a_zeta)
  if (zetasq>10000) zetasq=10000
  if (zetasq<0.000001) zetasq=0.000001
  zetasq
}
