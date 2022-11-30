# obtain log density of p(y|Theta)
get_log_den <- function(N,ncomp,Y,Zstar,beta_star,sigmasq,k,n){
  log_den=matrix(0,N,ncomp)
  for (r in 1:N){
    for (j in 1:ncomp){
      log_den[r,j] <- sum(dnorm(Y[r,],as.matrix(as.vector(apply(as.matrix(beta_star[,j,1:k]),2,function(l) Zstar %*% l))),
                                sqrt(rep(sigmasq[j,],each=n)),log = T))
    }
  }
  log_den
}
