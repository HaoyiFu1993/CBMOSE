# Obtain basis functions using eigen decomposition
get_Zstar <- function(n,x,m){
  Phi <- matrix(NA,nrow = n,ncol = n)
  for (i in 1:n){
    for (j in i:n){
      Phi[i,j] <- 0.5*(x[i]^2)*(x[j]-x[i]/3)
      Phi[j,i] <- Phi[i,j]
    }
  }
  X=matrix(0,n,2)
  X[,1]=rep(1,n)
  X[,2]=x
  decomp <- eigen(Phi,symmetric = T)
  m=10
  Z=decomp$vectors[,1:m] %*% (diag(decomp$values[1:m]))^0.5
  Zstar <- cbind(X,Z)
  Zstar
}
