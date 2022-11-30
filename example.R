##### an example of the proposed CBMOSE algorithm using simulated data
##### the 'example.RData' is simulated from a true model with two components (ncomp=2) and three dimensions/entries (k=3)

# import helper functions
source("CBMOSE.R")
source("extract_cov.R")
source("get_betastar.R")
source("get_delta.R")
source("get_delta_random.R")
source("get_init.R")
source("get_init_random.R")
source("get_log_den.R")
source("get_pi.R")
source("get_sigmasq.R")
source("get_tausq.R")
source("get_z.R")
source("get_zetasq.R")
source("get_Zstar.R")
source("MCMC_quantile.R")
source("plot_traj.R")

# import required library
library(BayesLogit)
library(mvtnorm)
library(LaplacesDemon)
library(dplyr)
library(ggplot2)

# load the example data
load("/home/haf48/R_functions/example.RData", verbose = T)

# run CBMOSE with time series data Y, covariate matrix V and time vector x
Y=example_data$Y
V=example_data$V
x=example_data$x

# set two components, three dimensions/entries, 10 basis functions, run 1000 iterations with 200 burnin and include random intercepts
set.seed(07069)
result <- CBMOSE(Y,V,x,ncomp=2,k=3,m=10,nloop = 1000,nburnin = 200,include.random = T)

# extract posterior mean, 2.5%, 50% (median) and 97.5$ quantiles of estimated parameters
param_quantile <- MCMC_quantile (MCMC_output = result,quantile=c(0.025,0.5,0.975))

# draw trajectory plots for selected component and dimension/entry
traj_res <- plot_traj(MCMC_output = result, use.median=T, include.CI=T, CI_level=c(0.025,0.975), ncomp=1, k=2, x=example_data$x,
                      xlab="Time",ylab="Posterior estimate",lwd=2)

# extract covariates based on mean and selected quantiles
cov_df <- extract_cov(result,CI_level = c(0.025,0.975),V=example_data$V)
cov_df
