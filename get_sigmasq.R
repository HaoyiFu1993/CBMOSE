# update error variance using Gibbs sampling 
get_sigmasq <- function(nu_sigma,sigmasq,A_sigma,Y,Zstar,beta_star,n){
  a_sigma <- rinvgamma(1,(nu_sigma+1)/2,nu_sigma/sigmasq + 1/A_sigma^2)
  epsilon <- t(apply(Y,1,function(s) s-Zstar%*%beta_star))
  epsilon2 =epsilon^2
  epsilonprod = apply(epsilon2,1,sum)
  epsilonprod2 = sum(epsilonprod)
  rinvgamma(1,(nu_sigma+n*nrow(Y))/2,epsilonprod2/2 + nu_sigma/a_sigma)
}
