##### Obtain mean and selected quantiles for each parameter
##### Always extract posterior mean from the MCMC result output
##### The default is to extract posterior mean, 2.5%, 50% (median) and 97.5% quantiles of the MCMC result output

##### Arguments of the function
# @param MCMC_output: the return list from the function CBMOSE
# @param quantile: a vector indicating selected quantiles

##### return a list of selected quantiles for each set of parameters
# \item betastar_quantile: mean and selected quantiles of intercepts, slopes and basis function coefficients for each component and
# dimension/entry
# \item sigmasq_quantile: mean and selected quantiles of error variances for each component and dimension/entry
# \item tausq_quantile: mean and selected quantiles of smoothing parameters for each component and dimension/entry
# \item delta_quantile: mean and selected quantiles of logistic parameters in multinomial logits for each component
# \item zeta_quantile: mean and selected quantiles of variances of random intercepts in multinomial logits for each component (if including
# random intercepts)


MCMC_quantile <- function(MCMC_output,quantile=c(0.025,0.5,0.975)){
  
  ##############################################################
  # extract mean and selected quantiles for beta_star          #
  ##############################################################
  
  # error check
  if (any(as.numeric(quantile)<0) | any(as.numeric(quantile)>1)){
    stop("Quantiles must be between 0 and 1")
  }
  
  # create the list for selected quantiles (always extract posterior mean)
  betastar_quantile <- vector("list",length = length(quantile)+1)
  names(betastar_quantile)[1] <- "beta_q_mean"
  for (i in 2:(length(quantile)+1)){
    names(betastar_quantile)[i] <- paste0("beta_q_",quantile[i-1])
  }
  
  # compute posterior mean
  d1 <- dim(MCMC_output$beta_star_save[[1]])[1]
  d2 <- dim(MCMC_output$beta_star_save[[1]])[2]
  d3 <- dim(MCMC_output$beta_star_save[[1]])[3]
  betastar_quantile[[1]] <- array(rep(NA,d1*d2*d3),dim = c(d1,d2,d3))
  for (i in 1:d1){
    for (j in 1:d2){
      for (q in 1:d3){
        betastar_quantile[[1]][i,j,q] <- mean(do.call("c",lapply(MCMC_output$beta_star_save,function(l) l[i,j,q])))
      }
    }
  }
  
  # compute selected quantiles
  for (r in 2:(length(quantile)+1)){
    betastar_quantile[[r]] <- array(rep(NA,d1*d2*d3),dim = c(d1,d2,d3))
  }
  for (r in 2:(length(quantile)+1)){
    for (i in 1:d1){
      for (j in 1:d2){
        for (q in 1:d3){
          betastar_quantile[[r]][i,j,q] <- quantile(do.call("c",lapply(MCMC_output$beta_star_save,function(l) l[i,j,q])),quantile[r-1])
        }
      }
    }
  }
  
  ##############################################################
  # extract mean and selected quantiles for sigmasq            #
  ##############################################################
  
  # create the list for selected quantiles (always extract posterior mean)
  sigmasq_quantile <- vector("list",length = length(quantile)+1)
  names(sigmasq_quantile)[1] <- "sigmasq_q_mean"
  for (i in 2:(length(quantile)+1)){
    names(sigmasq_quantile)[i] <- paste0("sigmasq_q_",quantile[i-1])
  }
  
  # compute posterior mean
  d4 <- dim(MCMC_output$sigmasq_save[[1]])[1]
  d5 <- dim(MCMC_output$sigmasq_save[[1]])[2]
  sigmasq_quantile[[1]] <- matrix(NA,nrow = d4,ncol = d5)
  for (i in 1:d4){
    for (j in 1:d5){
        sigmasq_quantile[[1]][i,j] <- mean(do.call("c",lapply(MCMC_output$sigmasq_save,function(l) l[i,j])))
    }
  }
  
  # compute selected quantiles
  for (r in 2:(length(quantile)+1)){
    sigmasq_quantile[[r]] <- matrix(NA,nrow = d4,ncol = d5)
  }
  for (r in 2:(length(quantile)+1)){
    for (i in 1:d4){
      for (j in 1:d5){
          sigmasq_quantile[[r]][i,j] <- quantile(do.call("c",lapply(MCMC_output$sigmasq_save,function(l) l[i,j])),quantile[r-1])
      }
    }
  }
  
  ##############################################################
  # extract mean and selected quantiles for tausq            #
  ##############################################################
  
  # create the list for selected quantiles (always extract posterior mean)
  tausq_quantile <- vector("list",length = length(quantile)+1)
  names(tausq_quantile)[1] <- "tausq_q_mean"
  for (i in 2:(length(quantile)+1)){
    names(tausq_quantile)[i] <- paste0("tausq_q_",quantile[i-1])
  }
  
  # compute posterior mean
  d4 <- dim(MCMC_output$tausq_save[[1]])[1]
  d5 <- dim(MCMC_output$tausq_save[[1]])[2]
  tausq_quantile[[1]] <- matrix(NA,nrow = d4,ncol = d5)
  for (i in 1:d4){
    for (j in 1:d5){
      tausq_quantile[[1]][i,j] <- mean(do.call("c",lapply(MCMC_output$tausq_save,function(l) l[i,j])))
    }
  }
  
  # compute selected quantiles
  for (r in 2:(length(quantile)+1)){
    tausq_quantile[[r]] <- matrix(NA,nrow = d4,ncol = d5)
  }
  for (r in 2:(length(quantile)+1)){
    for (i in 1:d4){
      for (j in 1:d5){
        tausq_quantile[[r]][i,j] <- quantile(do.call("c",lapply(MCMC_output$tausq_save,function(l) l[i,j])),quantile[r-1])
      }
    }
  }
  
  ##############################################################
  # extract mean and selected quantiles for delta              #
  ##############################################################
  
  # create the list for selected quantiles (always extract posterior mean)
  delta_quantile <- vector("list",length = length(quantile)+1)
  names(delta_quantile)[1] <- "delta_q_mean"
  for (i in 2:(length(quantile)+1)){
    names(delta_quantile)[i] <- paste0("delta_q_",quantile[i-1])
  }
  
  # compute posterior mean
  d6 <- dim(MCMC_output$delta_save[[1]])[1]
  d7 <- dim(MCMC_output$delta_save[[1]])[2]
  delta_quantile[[1]] <- matrix(NA,nrow = d6,ncol = d7)
  for (i in 1:d6){
    for (j in 1:d7){
      delta_quantile[[1]][i,j] <- mean(do.call("c",lapply(MCMC_output$delta_save,function(l) l[i,j])))
    }
  }
  
  # compute selected quantiles
  for (r in 2:(length(quantile)+1)){
    delta_quantile[[r]] <- matrix(NA,nrow = d6,ncol = d7)
  }
  for (r in 2:(length(quantile)+1)){
    for (i in 1:d6){
      for (j in 1:d7){
        delta_quantile[[r]][i,j] <- quantile(do.call("c",lapply(MCMC_output$delta_save,function(l) l[i,j])),quantile[r-1])
      }
    }
  }
  
  ################################################################################################
  # extract mean and selected quantiles for zetasq if including random intercepts                #
  ################################################################################################
  
  # if including random intercepts, give outputs of zetasq
  if("zetasq_save" %in% names(MCMC_output)){
    
    # create the list for selected quantiles (always extract posterior mean)
    zetasq_quantile <- vector("list",length = length(quantile)+1)
    names(zetasq_quantile)[1] <- "zetasq_q_mean"
    for (i in 2:(length(quantile)+1)){
      names(zetasq_quantile)[i] <- paste0("zetasq_q_",quantile[i-1])
    }
    
    # compute posterior mean
    d8 <- length(MCMC_output$zetasq_save[[1]])
    zetasq_quantile[[1]] <- rep(NA,d8)
    for (i in 1:d8){
        zetasq_quantile[[1]][i] <- mean(do.call("c",lapply(MCMC_output$zetasq_save,function(l) l[i])))
    }
    
    # compute selected quantiles
    for (r in 2:(length(quantile)+1)){
      zetasq_quantile[[r]] <- rep(NA,d8)
    }
    for (r in 2:(length(quantile)+1)){
      for (i in 1:d8){
          zetasq_quantile[[r]][i] <- quantile(do.call("c",lapply(MCMC_output$zetasq_save,function(l) l[i])),quantile[r-1])
      }
    }
    
  } else{
    # if not including random intercepts, give NULL for zetasq
    zetasq_quantile <- NULL
  }
  
  # make a list of statistics of all parameters
  param_quantile <- list(betastar_quantile=betastar_quantile,sigmasq_quantile=sigmasq_quantile,
                         tausq_quantile=tausq_quantile,delta_quantile=delta_quantile,zetasq_quantile=zetasq_quantile)
  param_quantile
  
  
}
