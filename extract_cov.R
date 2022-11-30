##### extract covariates with mean and selected quantiles

#### Arguments of function
# @param MCMC_output: the return list from the function CBMOSE
# @oaram CI_level: a vector of length 2 indicating the quantiles for credible intervals

##### return a data frame of covariates with mean and selected quantiles


extract_cov <- function(MCMC_output,CI_level=c(0.025,0.975),V){
  
  # error check
  if (any(as.numeric(CI_level)<0) | any(as.numeric(CI_level)>1)){
    stop("Quantiles must be between 0 and 1")
  }
  
  ##############################################################
  # extract mean and selected quantiles for delta              #
  ##############################################################
  
  # create the list for selected quantiles (always extract posterior mean)
  delta_quantile <- vector("list",length = length(CI_level)+2)
  names(delta_quantile)[1] <- "cov_mean"
  names(delta_quantile)[2] <- paste0("cov_q_",CI_level[1])
  names(delta_quantile)[3] <- "cov_median"
  names(delta_quantile)[4] <- paste0("cov_q_",CI_level[2])
  
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
  for (r in 2:(length(CI_level)+2)){
    delta_quantile[[r]] <- matrix(NA,nrow = d6,ncol = d7)
  }
  quantile_2 <- c(CI_level[1],0.5,CI_level[2])
  for (r in 2:(length(quantile_2)+1)){
    for (i in 1:d6){
      for (j in 1:d7){
        delta_quantile[[r]][i,j] <- quantile(do.call("c",lapply(MCMC_output$delta_save,function(l) l[i,j])),quantile_2[r-1])
      }
    }
  }
  
  # convert covariates of selected quantiles to the dataframe
  # whether include random intercepts
  if (nrow(delta_quantile[[1]])==(nrow(V)+ncol(V))){
    cov_df <- data.frame("term"= paste0("cov_",1:(nrow(V)+ncol(V))),
                         "mean"=as.vector(delta_quantile$cov_mean[,(1:(ncol(delta_quantile[[1]])-1))]),
                         "median"=as.vector(delta_quantile$cov_median[,(1:(ncol(delta_quantile[[1]])-1))]),
                         "CI.low"=as.vector(delta_quantile[[2]][,(1:(ncol(delta_quantile[[1]])-1))]),
                         "CI.high"=as.vector(delta_quantile[[4]][,(1:(ncol(delta_quantile[[1]])-1))]))
  } else{
    cov_df <- data.frame("term"= paste0("cov_",1:(ncol(V))),
                         "mean"=as.vector(delta_quantile$cov_mean[,(1:(ncol(delta_quantile[[1]])-1))]),
                         "median"=as.vector(delta_quantile$cov_median[,(1:(ncol(delta_quantile[[1]])-1))]),
                         "CI.low"=as.vector(delta_quantile[[2]][,(1:(ncol(delta_quantile[[1]])-1))]),
                         "CI.high"=as.vector(delta_quantile[[4]][,(1:(ncol(delta_quantile[[1]])-1))]))
  }
  
  cov_df
}
