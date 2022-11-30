# update intercept, slope and basis function coefficients using Gibbs sampling 
get_betastar <- function(V,Zstar,sigmasq,Y,m,PrVar_alpha,tausq){
  Sigma_beta <- matrix(0,m+2,m+2)
  Sigma_beta <- diag(c(PrVar_alpha,rep(tausq,m)))
  Gamma <- solve(dim(V)[1]*t(Zstar) %*% Zstar + sigmasq*solve(Sigma_beta))
  mu <- Gamma%*%(t(Zstar)%*%apply(Y,2,sum))
  rmvnorm(1,mu,sigmasq*Gamma)
}
