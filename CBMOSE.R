#####  Main function for the covariate-guided Bayesian mixture of spline experts model (CBMOSE)

#####  This function performs a Bayesian mixture of splines model and allow mixing weights to depend on values of selected covariates via
#####  multinomial logits. The splines model part creates basis functions based on the eigen decomposition. For the multinomial logits, there
#####  is an option to include random intercepts for each subject into the multinomial logits model.

#####  For more details, please refer to the article and supplemental materials of "Covariate-Guided Bayesian Mixture of Spline
#####  Experts for the Analysis of Multivariate Time Series"

#####  Arguments of the function
# @param Y: a (N*(k*n)) matrix of time series data, each row is a multivariate time series with dimension k, the length of each variate is n
# @param V: a (N*(p+1)) matrix of covariate matrix (the first column is a vector of 1 as the intercept)
# @param x: a vector of length n, indicating the time
# @param ncomp: an integer > 1, the number of components
# @param k: an integer >= 1, the number of dimensions/entries
# @param m: an integer >=1, the number of spline basis functions
# @param nloop: an integer, the number of iterations for the Gibbs sampling
# @param nburnin: an integer, the number of iterations that should drop out from the beginning
# @param init: a list including initial values for each set of parameters
# @param include.random: a logical value indicating if the random intercepts are included in the multinomial logits part
# @param hyperparam: a list including all hyperparameter values of priors

##### Other arguments
# @param N: an integer, the number of multivariate time series
# @param n: an integer, the number of time points for each univariate time series
# @param p: an integer, the number of covariates included in the multinomial logits

##### functions need to source
# source("get_betastar.R")
# source("get_delta.R")
# source("get_delta_random.R")
# source("get_init.R")
# source("get_init_random.R")
# source("get_log_den.R")
# source("get_pi.R")
# source("get_sigmasq.R")
# source("get_tausq.R")
# source("get_z.R")
# source("get_zetasq.R")
# source("get_Zstar.R")

##### packages need to import
# library(BayesLogit)
# library(mvtnorm)
# library(LaplacesDemon)
# library(dplyr)

##### return a list of posterior estimates of each parameter for each iteration
# \item pi_save: a list of length (nloop-nburnin), each with a (N*ncomp) matrix of mixing weights for each subject and component
# \item ind_save: a (N*(nloop-nburnin)) matrix of latent indicators for each subject and each iteration
# \item beta_star_save: a list of length (nloop-nburnin), each with a ((m+2)*ncomp*k) array of model coefficients (intercept, slope and basis
# function coefficients for each component and each dimension/entry
# \item sigmasq_save: a list of length (nloop-nburnin), each with a (ncomp*k) matrix of error variance for each component and dimension/entry
# \item tausq_save: a list of length (nloop-nburnin), each with a (ncomp*k) matrix of smoothing parameters 
# for each component and dimension/entry
# \item delta_save: a list of length (nloop-nburnin), if include.random=FALSE, each of delta_save is a ((p+1)*ncomp) matrix of logistic 
# parameters for each covariate and intercept, and each component. If include.random=TRUE, each of delta_save is a ((p+1+N)*ncomp) matrix of
# logistic parameters and random intercepts. The last component is used as the reference component so the last column of each element of 
# delta_save are all zero
# \item zetasq_save: If include.random=TRUE, zetasq_save is a list of length (nloop-nburnin), each with a vector of length ncomp. Since
# the last component is used as the reference component, the last value of each element of zetasq_save is set to 1
# \item prob_post_save: a list of length (nloop-nburnin), each with a (N*ncomp) matrix of probabilities of latent indicators for each subject and
# each component
# \item log_den_save: a list of length (nloop-nburnin), each with a (N*ncomp) matrix of log density of each subject given parameter estimates
# from each component p(y_i | theta_g), y_i is the multivariate time series and theta_g includes all parameters for the gth component
# \item Z_star: m basis functions using eigen decomposition



CBMOSE <- function(Y,V,x,ncomp,k,m,nloop=10000,nburnin=2000,init=NULL,include.random=F,hyperparam=NULL){
  
  ###################################
  # Error checking                  #
  ###################################
  
  ncomp=as.integer(ncomp)
  n=length(x)
  N=nrow(V)
  ndelta=ncol(V)
  if(ncomp<2){
    stop ("The number of component must be an integer greater than 1")}
  if(nloop <= 0){
    stop("The number of iterations has to be positive.")}
  if(nburnin <= 0){
    stop("The number of iterations has to be positive.")}
  if(nrow(V) != nrow(Y)){
    stop("The number of rows of V should equal to the number of rows of Y")}
  if (k*n != ncol(Y)){
    stop("Dimension times time points should equal to the number of columns of Y")}
  
  
  ###################################
  # Store posterior values          #
  ###################################
  
  pi_save=vector("list",nloop)
  ind_save=matrix(0,N,nloop)
  # prob_post_save=vector("list",nloop)
  delta_save <- vector("list",nloop)
  beta_star_save <- vector("list",nloop)
  sigmasq_save <- vector("list",nloop)
  tausq_save <- vector("list",nloop)
  zetasq_save <- vector("list",nloop)
  # probabilities of latent indicators
  prob_post_save=vector("list",nloop)
  # obtain density of p(y|Theta)
  log_den_save=vector("list",nloop)
  
  ######################################################################
  # Obtain m basis functions using eigen decomposition                 #
  ######################################################################
  
  Zstar <- get_Zstar(n,x,m)
  
  ######################################################################
  # not include random intercept in multinomial logits                 #
  ######################################################################
  
  if(include.random==F){
    
    # if init is null, the default is to generate random values for model coefficients, and set all variance term to 1
    if(is.null(init)){
      # set up initial values
      init <- get_init(N,ncomp,k,m,V)
      delta <- init$delta
      beta_star <- init$beta_star
      sigmasq <- init$sigmasq
      tausq <- init$tausq
      pi <- init$pi
      ind <- init$ind
      z <- init$z
    } else{
      # check the dimension of all initial values
      init_check <- get_init(N,ncomp,k,m,V)
      if (identical(sapply(init_check, function(l) dim(l)),sapply(init, function(l) dim(l)))==F){
        stop("Dimensions of initial values are incorrect")}
    }
    
    # if hyperparam is null, use default hyperparameter values
    if(is.null(hyperparam)){
      # set up default hyperparameter value
      nu_tau=nu_sigma=3; A_tau=A_sigma=10  #hyperparameters of prior on sigmasq and tausq, assuming same across components and dimensions
      PrVar_alpha=rep(100,2)   # prior variance for the intercept and slope,
      PrVar_delta=diag(rep(10,ndelta))   # prior variance for each of the logistic parameter
      hyperparam <- list(nu_tau=nu_tau,nu_sigma=nu_sigma,A_tau=A_tau,A_sigma=A_sigma,PrVar_alpha=PrVar_alpha,PrVar_delta=PrVar_delta)
    } else{
      # extract each hyperparameter in an order
      nu_tau=hyperparam[[1]]
      nu_sigma=hyperparam[[2]]
      A_tau=hyperparam[[3]]
      A_sigma=hyperparam[[4]]
      PrVar_alpha=hyperparam[[5]]
      PrVar_delta=hyperparam[[6]]
      # if only one value is given, use it for all parameters
      if (length(PrVar_alpha)==1){
        PrVar_alpha <- rep(PrVar_alpha,2)
      }
      if (length(PrVar_delta)==1){
        PrVar_delta <- diag(rep(PrVar_delta,ndelta))
      }
    }
    
    ######################################################################
    # run Gibbs sampling                                                 #
    ######################################################################
    
    for (p in 1:nloop){
      cat("p:", p, " \n")
      
      # fit splines model for each component and each dimension/entry
      for (j in 1:ncomp){
        Y_exp=Y[ind==j,]
        V_exp=V[ind==j,]
        if (length(V_exp)==ndelta) V_exp=matrix(V_exp,ncol=ndelta)
        if (length(Y_exp)==n*k) Y_exp=matrix(Y_exp,ncol=n*k)
        
        # within each dimension/entry, sample intercept, slope and basis function coefficients beta_star
        for (i in 1:k){
          beta_star[,j,i] <- get_betastar(V_exp,Zstar,sigmasq[j,i],Y_exp[,(n*(i-1)+1):(n*i),drop=F],m,PrVar_alpha,tausq[j,i])
        }
        
        # within each dimension/entry, sample error variance sigmasq
        for (i in 1:k){
          sigmasq[j,i] <- get_sigmasq(nu_sigma,sigmasq[j,i],A_sigma,Y_exp[,(n*(i-1)+1):(n*i),drop=F],Zstar,beta_star[,j,i],n)
        }
        
        # within each dimension/entry, sample smoothing parameter tausq
        for (i in 1:k){
          tausq[j,i] <- get_tausq(nu_tau=3,tausq[j,i],A_tau,m,beta=beta_star[3:(m+2),j,i])
        }
      }
      
      # sample logistic parameter delta
      if (ncomp==1){
        # for number of components equal to 1
        delta <- rep(0,ndelta)
        pi <- rep(1,N)
        ind <- rep(1,N)
        z <- rep(1,N)
        
      } else{
        # for number of components larger than 1
        z=matrix(0,N,ncomp)
        for (j in 1:ncomp){
          z[ind==j,j]=1
        }
        c=matrix(0,N,ncomp)
        exb=exp(V %*% delta)
        exb[which(is.infinite(exb))]=.Machine$double.xmax
        for (j in 1:(ncomp-1)){
          if(ncomp==2){c[,j]=rep(0,N)
          } else{c[,j]=log(apply(exb[,-j],1,sum))}
        }
        c[which(is.infinite(c))]=.Machine$double.xmax
        for (j in 1:(ncomp-1)){
          delta[,j] <- get_delta(N,V,c[,j],delta[,j],z[,j],PrVar_delta)
        }
        
        
        # compute mixing weights pi
        pi <- get_pi(V,delta)
        for (j in 1:ncomp){
          pi[,j][which(pi[,j]==0)] <- .Machine$double.xmin
        }
        
        # compute log density of p(y|Theta)
        log_den <- get_log_den(N,ncomp,Y,Zstar,beta_star,sigmasq,k,n)
        
        # computing probability of indicators z
        prob_post <- get_z(N,ncomp,pi,Y,Zstar,beta_star,sigmasq,k,n)
        
        # sample z based on the probabilities of latent indicators
        ind=rep(NA,N)
        for (i in 1:N){
          ind[i] <- sample(1:ncomp,1,replace = T,prob = prob_post[i,])
        }
        # obtain the matrix of latent indicators
        z=matrix(0,N,ncomp)
        for (j in 1:ncomp){
          z[ind==j,j]=1
        }
      }
      
      ######################################################################
      # store values for each iteration                                    #
      ######################################################################
      
      pi_save[[p]] <- pi
      ind_save[,p] <- ind
      beta_star_save[[p]] <- beta_star
      sigmasq_save[[p]] <- sigmasq
      tausq_save[[p]] <- tausq
      delta_save[[p]] <- delta
      prob_post_save[[p]] <- prob_post
      log_den_save[[p]] <- log_den
    }
    
    ######################################################################
    # remove burnin results                                              #
    ######################################################################
    
    result <- list(pi_save=pi_save[(nburnin+1):nloop],ind_save=ind_save[,((nburnin+1):nloop)],
                   beta_star_save=beta_star_save[(nburnin+1):nloop],sigmasq_save=sigmasq_save[(nburnin+1):nloop],
                   tausq_save=tausq_save[(nburnin+1):nloop],delta_save=delta_save[(nburnin+1):nloop],
                   prob_post_save=prob_post_save[(nburnin+1):nloop], log_den_save=log_den_save[(nburnin+1):nloop],Zstar=Zstar)
    
    
  }
  
  ######################################################################
  # include random intercept in multinomial logits                     #
  ######################################################################
  
  if(include.random==T){
    
    V_2 <- cbind(V,diag(1,N,N))
    # if init is null, the default is to generate random values for model coefficients, and set all variance term to 1
    if(is.null(init)){
      # set up initial values
      init <- get_init_random(N,ncomp,k,m,V)
      delta <- init$delta   # delta also includes random intercept for each subject
      beta_star <- init$beta_star
      sigmasq <- init$sigmasq
      tausq <- init$tausq
      zetasq <- init$zetasq
      pi <- init$pi
      ind <- init$ind
      z <- init$z
    } else{
      # check the dimension of all initial values
      init_check <- get_init_random(N,ncomp,k,m,V)
      if (identical(sapply(init_check, function(l) dim(l)),sapply(init, function(l) dim(l)))==F){
        stop("Dimensions of initial values are incorrect")}
    }
    
    if(is.null(hyperparam)){
      # set up default hyperparameter value
      
      #hyperparameters of prior on sigmasq,tausq and zetasq, assuming same across components and dimensions
      nu_tau=nu_sigma=nu_zeta=3; A_tau=A_sigma=A_zeta=10  
      PrVar_alpha=rep(100,2)   # prior variance for the intercept and slope,
      PrVar_delta=diag(rep(10,ndelta))   # prior variance for each of the logistic parameter
      hyperparam <- list(nu_tau=nu_tau,nu_sigma=nu_sigma,nu_zeta=nu_zeta,A_tau=A_tau,A_sigma=A_sigma,A_zeta=A_zeta,
                         PrVar_alpha=PrVar_alpha,PrVar_delta=PrVar_delta)
    } else{
      # extract each hyperparameter in an order
      nu_tau=hyperparam[[1]]
      nu_sigma=hyperparam[[2]]
      nu_zeta=hyperparam[[3]]
      A_tau=hyperparam[[4]]
      A_sigma=hyperparam[[5]]
      A_zeta=hyperparam[[6]]
      PrVar_alpha=hyperparam[[7]]
      PrVar_delta=hyperparam[[8]]
      # if only one value is given, use it for all parameters
      if (length(PrVar_alpha)==1){
        PrVar_alpha <- rep(PrVar_alpha,2)
      }
      if (length(PrVar_delta)==1){
        PrVar_delta <- diag(rep(PrVar_delta,ndelta))
      }
    }
    
    ######################################################################
    # run Gibbs sampling                                                 #
    ######################################################################
    
    for (p in 1:nloop){
      cat("p:", p, " \n")
      
      # fit splines model for each component and each dimension/entry
      for (j in 1:ncomp){
        Y_exp=Y[ind==j,]
        V_exp=V_2[ind==j,]
        if (length(V_exp)==ndelta) V_exp=matrix(V_exp,ncol=ndelta)
        if (length(Y_exp)==n*k) Y_exp=matrix(Y_exp,ncol=n*k)
        
        # within each dimension/entry, sample intercept, slope and basis function coefficients beta_star
        for (i in 1:k){
          beta_star[,j,i] <- get_betastar(V_exp,Zstar,sigmasq[j,i],Y_exp[,(n*(i-1)+1):(n*i),drop=F],m,PrVar_alpha,tausq[j,i])
        }
        
        # within each dimension/entry, sample error variance sigmasq
        for (i in 1:k){
          sigmasq[j,i] <- get_sigmasq(nu_sigma,sigmasq[j,i],A_sigma,Y_exp[,(n*(i-1)+1):(n*i),drop=F],Zstar,beta_star[,j,i],n)
        }
        
        # within each dimension/entry, sample smoothing parameter tausq
        for (i in 1:k){
          tausq[j,i] <- get_tausq(nu_tau=3,tausq[j,i],A_tau,m,beta=beta_star[3:(m+2),j,i])
        }
      }
      
      # sample logistic parameter delta
      if (ncomp==1){
        # for number of components equal to 1
        delta <- rep(0,ndelta)
        pi <- rep(1,N)
        ind <- rep(1,N)
        z <- rep(1,N)
        
      } else{
        # for number of components larger than 1
        z=matrix(0,N,ncomp)
        for (j in 1:ncomp){
          z[ind==j,j]=1
        }
        c=matrix(0,N,ncomp)
        exb=exp(V_2 %*% delta)
        exb[which(is.infinite(exb))]=.Machine$double.xmax
        for (j in 1:(ncomp-1)){
          if(ncomp==2){c[,j]=rep(0,N)
          } else{c[,j]=log(apply(exb[,-j],1,sum))}
        }
        c[which(is.infinite(c))]=.Machine$double.xmax
        for (j in 1:(ncomp-1)){
          delta[,j] <- get_delta_random(N,N,V_2,c[,j],delta[,j],z[,j],PrVar_delta,ndelta,zetasq[j])
          # sample variance of random intercept zetasq
          zetasq[j] <- get_zetasq(nu_zeta,zetasq[j],A_zeta,m,delta[(ndelta+1):(N+ndelta),j],N)
        }
        
        
        # compute mixing weights pi
        pi <- get_pi(V_2,delta)
        for (j in 1:ncomp){
          pi[,j][which(pi[,j]==0)] <- .Machine$double.xmin
        }
        
        # compute log density of p(y|Theta)
        log_den <- get_log_den(N,ncomp,Y,Zstar,beta_star,sigmasq,k,n)
        
        # computing probability of indicators z
        prob_post <- get_z(N,ncomp,pi,Y,Zstar,beta_star,sigmasq,k,n)
        
        # sample z based on the probabilities of latent indicators
        ind=rep(NA,N)
        for (i in 1:N){
          ind[i] <- sample(1:ncomp,1,replace = T,prob = prob_post[i,])
        }
        # obtain the matrix of latent indicators
        z=matrix(0,N,ncomp)
        for (j in 1:ncomp){
          z[ind==j,j]=1
        }
      }
      
      ######################################################################
      # store values for each iteration                                    #
      ######################################################################
      pi_save[[p]] <- pi
      ind_save[,p] <- ind
      beta_star_save[[p]] <- beta_star
      sigmasq_save[[p]] <- sigmasq
      tausq_save[[p]] <- tausq
      zetasq_save[[p]] <- zetasq
      delta_save[[p]] <- delta
      prob_post_save[[p]] <- prob_post
      log_den_save[[p]] <- log_den
    }
    
    ######################################################################
    # remove burnin results                                              #
    ######################################################################
    
    result <- list(pi_save=pi_save[(nburnin+1):nloop],ind_save=ind_save[,((nburnin+1):nloop)],
                   beta_star_save=beta_star_save[(nburnin+1):nloop],sigmasq_save=sigmasq_save[(nburnin+1):nloop],
                   tausq_save=tausq_save[(nburnin+1):nloop],delta_save=delta_save[(nburnin+1):nloop],
                   zetasq_save=zetasq_save[(nburnin+1):nloop],
                   prob_post_save=prob_post_save[(nburnin+1):nloop], log_den_save=log_den_save[(nburnin+1):nloop], Zstar=Zstar)
    
  }
  return(result)
}
