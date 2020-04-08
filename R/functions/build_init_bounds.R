## Functions.R
## M S KAPUR adapted from J Sullivan seak_sablefish & Grant Adams
## preset data gen, post-processing and plotting functions for use in *_Master.R scripts

# Build parameter bounds
# Original code by Grant Adams, adapted for use with the sablefish model

build_bounds <- function(param_list = NULL, data_list){
  
  # Debug function
  # param_list <- parameters
  # data_list <- data
  
  upper_bnd <- param_list
  lower_bnd <- param_list
  
  # General bounds
  for(i in 1:length(param_list)){
    upper_bnd[[i]] <- replace(upper_bnd[[i]], values = rep(Inf, length(upper_bnd[[i]])))
    lower_bnd[[i]] <- replace(lower_bnd[[i]], values = rep(-Inf, length(lower_bnd[[i]])))
  }
  
  # Natural mortality
  lower_bnd$log_M <- replace(lower_bnd$log_M, values = -3) 
  upper_bnd$log_M <- replace(upper_bnd$log_M, values = -1) 
  
  # Fishery selectivity
  lower_bnd$log_fsh_slx_pars[,,] <- replace(lower_bnd$log_fsh_slx_pars[,,], values = -2) 
  upper_bnd$log_fsh_slx_pars[,,] <- replace(upper_bnd$log_fsh_slx_pars[,,], values = 2) 
  
  # Survey selectivity
  lower_bnd$log_srv_slx_pars[,,] <- replace(lower_bnd$log_srv_slx_pars[,,], values = -2) 
  upper_bnd$log_srv_slx_pars[,,] <- replace(upper_bnd$log_srv_slx_pars[,,], values = 2) 
  
  # Fishery catchability
  lower_bnd$fsh_logq <- replace(lower_bnd$fsh_logq, values = rep(-30, length(lower_bnd$fsh_logq)))
  upper_bnd$fsh_logq <- replace(upper_bnd$fsh_logq, values = rep(5, length(upper_bnd$fsh_logq)))
  
  # Survey catchability
  lower_bnd$srv_logq <- replace(lower_bnd$srv_logq, values = rep(-30, length(lower_bnd$srv_logq)))
  upper_bnd$srv_logq <- replace(upper_bnd$srv_logq, values = rep(5, length(upper_bnd$srv_logq)))
  
  # Mark-recapture catchability
  lower_bnd$mr_logq <- replace(lower_bnd$mr_logq, values = rep(-20, length(lower_bnd$mr_logq)))
  upper_bnd$mr_logq <- replace(upper_bnd$mr_logq, values = rep(1, length(upper_bnd$mr_logq)))
  
  # Recruitment devs
  lower_bnd$log_rec_devs <- replace(lower_bnd$log_rec_devs, values = rep(-10, length(lower_bnd$log_rec_devs)))
  upper_bnd$log_rec_devs <- replace(upper_bnd$log_rec_devs, values = rep(10, length(upper_bnd$log_rec_devs)))
  
  # Initial numbers-at-age devs
  lower_bnd$log_rinit_devs <- replace(lower_bnd$log_rinit_devs, values = rep(-10, length(lower_bnd$log_rinit_devs)))
  upper_bnd$log_rinit_devs <- replace(upper_bnd$log_rinit_devs, values = rep(10, length(upper_bnd$log_rinit_devs)))
  
  # F devs
  lower_bnd$log_F_devs <- replace(lower_bnd$log_F_devs, values = rep(-15, length(lower_bnd$log_F_devs)))
  upper_bnd$log_F_devs <- replace(upper_bnd$log_F_devs, values = rep(15, length(upper_bnd$log_F_devs)))
  
  # SPR F rates
  lower_bnd$log_spr_Fxx <- replace(lower_bnd$log_spr_Fxx, values = rep(-5, length(lower_bnd$log_spr_Fxx)))
  upper_bnd$log_spr_Fxx <- replace(upper_bnd$log_spr_Fxx, values = rep(5, length(upper_bnd$log_spr_Fxx)))
  
  # Fishery age comp Dirichlet-multinomial theta
  lower_bnd$log_fsh_theta <- replace(lower_bnd$log_fsh_theta, values = rep(-5, length(lower_bnd$log_fsh_theta)))
  upper_bnd$log_fsh_theta <- replace(upper_bnd$log_fsh_theta, values = rep(15, length(upper_bnd$log_fsh_theta)))
  
  # Fishery age comp Dirichlet-multinomial theta
  lower_bnd$log_srv_theta <- replace(lower_bnd$log_srv_theta, values = rep(-5, length(lower_bnd$log_srv_theta)))
  upper_bnd$log_srv_theta <- replace(upper_bnd$log_srv_theta, values = rep(15, length(upper_bnd$log_srv_theta)))
  
  # Put bounds together
  bounds <- list(upper= upper_bnd, lower = lower_bnd)
  
  # Make sure inits are within bounds
  if( sum(unlist(bounds$upper) < as.numeric(unlist(param_list))) > 0 | sum(as.numeric(unlist(param_list)) < unlist(bounds$lower)) > 0 ){
    lower_check <- param_list
    upper_check <- param_list
    param_check <- data.frame(matrix(NA, nrow = length(param_list), ncol = 3))
    colnames(param_check) <- c("Parameter", "Lower", "Upper")
    param_check$Parameter <- names(param_list)
    
    for(i in 1:length(param_list)){
      lower_check[[i]] <- param_list[[i]] < lower_bnd[[i]]
      upper_check[[i]] <- param_list[[i]] > upper_bnd[[i]]
      param_check$Lower[i] <- sum(lower_check[[i]])
      param_check$Upper[i] <- sum(upper_check[[i]])
    }
    
    print("Non-zero value indicates error in initial value")
    print(param_check)
    stop("Initial parameter values are not within bounds")
  }
  
  return(bounds)
}