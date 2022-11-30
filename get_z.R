# update latent indicators z using Gibbs sampling
get_z <- function(N,ncomp,pi,Y,Zstar,beta_star,sigmasq,k,n){
  # compute log pi times density for each subject and each component
  prob_post=matrix(0,N,ncomp)
  for (r in 1:N){
    for (j in 1:ncomp){
      prob_post[r,j] <- log(pi[r,j]) + sum(dnorm(Y[r,],as.matrix(as.vector(apply(as.matrix(beta_star[,j,1:k]),2,function(l) Zstar %*% l))),
                                                 sqrt(rep(sigmasq[j,],each=n)),log = T))
    }
  }
  # obtain the max density of component for each subject
  pif_max <- rep(NA,N)
  for (r in 1:N){
    pif_max[r] <- which.max(prob_post[r,])
  }
  # compute probability of z using the LogExpSum trick
  for (r in 1:N){
    for (j in 1:ncomp){
      prob_post[r,j] <- exp(log(pi[r,j])-log(pi[r,pif_max[r]])+
                              sum(dnorm(Y[r,],as.matrix(as.vector(apply(as.matrix(beta_star[,j,1:k]),2,function(l) Zstar %*% l))),
                                        sqrt(rep(sigmasq[j,],each=n)),log = T))-
                              sum(dnorm(Y[r,],as.matrix(as.vector(apply(as.matrix(beta_star[,pif_max[r],1:k]),2,function(l) Zstar %*% l))),
                                        sqrt(rep(sigmasq[pif_max[r],],each=n)),log = T)))
    }
  }
  # reset the threshold for very small number
  for (r in 1:N){
    for (j in 1:ncomp){
      if (prob_post[r,j] < 2.220446e-7) {prob_post[r,j]=2.220446e-7}
    }
  }
  # compute probability of z
  prob_post/apply(prob_post,1,sum)
}
