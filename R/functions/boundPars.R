boundPars <- function(obj, r0_lower = 10, boundSlx = FALSE){
  
  ## bounds on repro pars ----
  
  lower <- obj$par-Inf
  upper <- obj$par+Inf
  
  lower[names(lower) == 'logh_k'] <- log(0.2) ## duh
  upper[names(upper) == 'logh_k'] <- log(0.99)
  # lower[names(lower) == 'epsilon_tau'] <- 0
  lower[names(lower) == 'logR_0k'] <- r0_lower
  # 
  # lower[names(lower) == 'logSDR'] <- log(0.0001)
  # lower[names(lower) == 'omega_0ij'] <- log(0.0001)
  # upper[names(upper) == 'omega_0ij'] <- log(0.999)
  # 
  # lower[names(lower) == 'b'] <- 0 
  # upper[names(upper) == 'b'] <- 1
  ## bounds on fsh slx ----
  ## first check if slx was fixed at all
  # if(length(grep("log_fsh_slx_pars", mappy)) != 0){
  #   ## make a master array with everything
  #   array(1:length(obj$par[names(obj$par) == "log_fsh_slx_pars"]), dim = c(7,2,2))
  #   ## identify which fleets were NA
  #   which(is.na(mappy[[grep("log_fsh_slx_pars", mappy)]]))
  #   
  # }
  
  ## brute force par locations; if mappy mirrors first ak fleets,
  ## assume that fleets 1 and 2 are now AK FIX mirror and AK TWL mirror
  # array(1:length(obj$par[names(obj$par) == "log_fsh_slx_pars"]), dim = c(4,2,2))
  
  
  if(boundSlx == TRUE){
    ## if no fleets to fix, go with the normal slx bounds
    p1_logistic_idx <- c(1:2,15:16) #c(1:4,19:22)
    p2_logistic_idx <- p1_logistic_idx+df$nfleets_fish
    p1_norm_idx <- c(3,4,6,7,17,18,20,21)  #c(5,6,8,9,23,24,26,27)
    p2_norm_idx <- p1_norm_idx+df$nfleets_fish
    p1_gamma_idx <- c(5,19)#c(7,25)
    p2_gamma_idx <- p1_gamma_idx+df$nfleets_fish
    # #* fsh slx lower bounds ----
    # ## logistic p1 (a50)
    # array(1:length(lower[names(lower) == 'log_fsh_slx_pars']), dim = c(7,2,2))
    
    lower[names(lower) == 'log_fsh_slx_pars'][p1_logistic_idx] <- log(35)
    ## logistic p2 (a95)
    lower[names(lower) == 'log_fsh_slx_pars'][p2_logistic_idx] <- log(60)
    ## normal p1 (mean)
    lower[names(lower) == 'log_fsh_slx_pars'][p1_norm_idx] <- log(15)
    ## normal p2 (sd)
    lower[names(lower) == 'log_fsh_slx_pars'][p2_norm_idx] <- log(1)
    ## gamma shape (k*theta equals mean)
    lower[names(lower) == 'log_fsh_slx_pars'][p1_gamma_idx] <- log(15)
    ## gamma rate
    lower[names(lower) == 'log_fsh_slx_pars'][p2_gamma_idx] <- log(2)
    # #* fsh slx upper bounds ----
    # ## logistic p1 (a50)
    upper[names(upper) == 'log_fsh_slx_pars'][p1_logistic_idx] <- log(60)
    ## logistic p2 (a95)
    upper[names(upper) == 'log_fsh_slx_pars'][p2_logistic_idx] <- log(70)
    ## normal p1 (mean)
    upper[names(upper) == 'log_fsh_slx_pars'][p1_norm_idx] <- log(65)
    ## normal p2 (sd)
    upper[names(upper) == 'log_fsh_slx_pars'][p2_norm_idx] <- log(15)
    ## gamma shape (k*theta equals mean)
    upper[names(upper) == 'log_fsh_slx_pars'][p1_gamma_idx] <- log(35)
    ## gamma rate
    upper[names(upper) == 'log_fsh_slx_pars'][p2_gamma_idx] <- log(2)
    
  }
  # array(exp(upper[names(upper) == 'log_fsh_slx_pars']), dim = dim(df$parms$log_fsh_slx_pars))
  # array(exp(lower[names(lower) == 'log_fsh_slx_pars']), dim = dim(df$parms$log_fsh_slx_pars))
  # # array(lower[names(lower) == 'log_fsh_slx_pars'], dim = dim = c(7,2,2))
  
  ## currently srv slx all logistic with a95, a50
  nsurvsel = dim(df$parms$log_srv_slx_pars)[1]
  ## lower for everything
  lower[names(lower) == 'log_srv_slx_pars'] <- log(0.0001)
  ## lower bound for p2 (a95)
  lower[names(lower) == 'log_srv_slx_pars'][c(c(1:nsurvsel,17:(16+nsurvsel))+nsurvsel)] <- log(70)
  ## upper bound for p1 (a50 or mean)
  upper[names(upper) == 'log_srv_slx_pars'][c(c(1:nsurvsel,17:(16+nsurvsel)))] <- log(70)
  ## upper bound for p2 (a95)
  upper[names(upper) == 'log_srv_slx_pars'][c(c(1:nsurvsel,17:(16+nsurvsel))+nsurvsel)]  <- log(70)
  lower[names(lower) == 'omega_0ij'] = 0
  upper[names(upper) == 'omega_0ij'] = 1
  
  ## sanity check
  
  ## last five flts p2 should be 10; p2 for first 4 fleets should be > p1
  
  ##  p2 should be > p1
  array(exp(upper[names(upper) == 'log_srv_slx_pars']), dim = dim(df$parms$log_srv_slx_pars))
  ## all zero and/or
  
  array(exp(lower[names(lower) == 'log_srv_slx_pars']), dim = dim(df$parms$log_srv_slx_pars))
  # upper[names(upper) == 'PSEL'] <- 9
  # upper[names(upper) == 'logh'] <- log(0.999)
  # upper[names(upper) == 'F0'] <- 2
  return(list("upper"=upper, "lower"=lower))
}
