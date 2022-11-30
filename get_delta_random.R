# update logistic parameters and random intercepts using Gibbs sampling and Polya-Gamma augmentation method
get_delta_random <- function(N,N_split,V_2,c,delta,z,PrVar_delta,ndelta,zetasq){
  # Update latents
  Sigma_delta <- matrix(0,N+ndelta,N+ndelta)
  if (zetasq==0){
    Sigma_delta <- diag(c(rep(PrVar_delta[1,1],ndelta),rep(0.00001,N)))
  } else{
    Sigma_delta <- diag(c(rep(PrVar_delta[1,1],ndelta),rep(zetasq,N)))
  }
  w = rpg(N_split, rep(1,N_split), V_2 %*% delta- c)
  nu = w*c + z - .5
  # Update regression coefficients
  Vd = solve(t(V_2) %*% diag(w) %*% V_2 + solve(Sigma_delta))
  mud = Vd %*% (t(V_2) %*% nu)
  t(rmvnorm(1,mud,Vd))
}
