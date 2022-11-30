##### Trajectory plot for selected component and dimension/entry
##### Credible intervals for trajectory are computed based on selected quantiles of trajectory estimates across all iterations

#### Arguments of function
# @param MCMC_output: the return list from the function CBMOSE
# @param use.median: a logical value indicating if posterior median is used. The default is false and use posterior mean
# @param include.CI: a logical value indicating if credible intervals are drawn for the trajectory plot. The default is false
# @oaram CI_level: a vector of length 2 indicating the quantiles for credible intervals
# @param ncomp: an integer > 1, the number of components
# @param k: an integer >= 1, the number of dimensions/entries
# @param x: a vector of length n, indicating the time

##### return a list of length two
# \item post_est: a ((nloop-nburnin)*n) matrix of trajectory estimates of all iterations for selected component and dimension/entry
# \item traj_quantile: a list of length 4 with mean, median and two credible intervals of trajectory estimates across all iterations for
# selected component and dimension/entry

plot_traj <- function(MCMC_output, use.median=F, include.CI=F, CI_level=c(0.025,0.975), ncomp, k, x, ...){
  
  # extract basis functions 
  Zstar <- MCMC_output$Zstar
  
  # error check
  if (length(CI_level) != 2){
    stop("The level of credible interval must equal to 2")
  }
  
  # compute posterior estimates of trajectory for selected component and dimension/entry for all iterations
  post_est <- matrix(NA,ncol = nrow(Zstar),nrow = length(MCMC_output$beta_star_save))
  for (i in 1:length(MCMC_output$beta_star_save)){
    post_est[i,] <- Zstar %*% MCMC_output$beta_star_save[[i]][,ncomp,k]
  }
  
  # compute posterior mean, median and selected credible intervals based on trajectory estimates of all iterations
  traj_quantile <- vector("list",length = length(CI_level)+2)
  names(traj_quantile)[1] <- "traj_mean"
  names(traj_quantile)[2] <- paste0("traj_q_",CI_level[1])
  names(traj_quantile)[3] <- "traj_median"
  names(traj_quantile)[4] <- paste0("traj_q_",CI_level[2])
  traj_quantile[[1]] <- apply(post_est,2,mean)
  traj_quantile[[2]] <- apply(post_est,2,quantile,CI_level[1])
  traj_quantile[[3]] <- apply(post_est,2,median)
  traj_quantile[[4]] <- apply(post_est,2,quantile,CI_level[2])
  
  
  # if use.median is false, posterior mean will be used to draw trajectory plot
  if (use.median==F){
    # if include.CI is false, draw trajectory plot without credible intervals 
    if (include.CI==F){
      plot(x,traj_quantile$traj_mean,type = "l",...)
    } else{
      # if include.CI is TRUE, draw trajectory plot with credible intervals of selected levels
      plot(x,traj_quantile$traj_mean,type = "l",ylim = 
             c(min(sapply(traj_quantile,min)),max(sapply(traj_quantile,max))),...)
      lines(x,traj_quantile[[2]],col=2,lwd=1.5,lty=2)
      lines(x,traj_quantile[[4]],col=2,lwd=1.5,lty=2)
    }
  }
  
  # if use.median is true, posterior median will be used to draw trajectory plot
  if (use.median==T){
    # if include.CI is false, draw trajectory plot without credible intervals 
    if (include.CI==F){
      plot(x,traj_quantile$traj_median,type = "l",...)
    } else{
      # if include.CI is TRUE, draw trajectory plot with credible intervals of selected levels
      plot(x,traj_quantile$traj_median,type = "l",ylim = 
             c(min(sapply(traj_quantile,min)),max(sapply(traj_quantile,max))),...)
      lines(x,traj_quantile[[2]],col=2,lwd=1.5,lty=2)
      lines(x,traj_quantile[[4]],col=2,lwd=1.5,lty=2)
    }
  }
  
  post_est_result <- list(post_est=post_est,traj_quantile=traj_quantile)
  post_est_result
}
