# update logistic parameters using Gibbs sampling and Polya-Gamma augmentation method
get_delta <- function(N,V,c,delta,z,PrVar_delta){
  # Update latents
  w = rpg(N, rep(1,N), V %*% delta- c)
  nu = w*c + z - .5
  # Update regression coefficients
  Vd = solve(t(V) %*% diag(w) %*% V + solve(PrVar_delta))
  mud = Vd %*% (t(V) %*% nu)
  t(rmvnorm(1,mud,Vd))
}
