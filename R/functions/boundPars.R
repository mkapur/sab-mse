boundPars <- function(obj, r0_lower = 10, boundSlx = c('fsh','srv')){
  
  ## bounds on repro pars ----
  
  lower <- obj$par-Inf
  upper <- obj$par+Inf
  
  lower[names(lower) == 'log_srv_slx_pars'] <-   lower[names(lower) == 'log_fsh_slx_pars'] <- 0
  upper[names(upper) == 'log_srv_slx_pars'] <-   upper[names(upper) == 'log_fsh_slx_pars'] <- log(100)
  
  lower[names(lower) == 'logh_k'] <- log(0.2) ## duh
  upper[names(upper) == 'logh_k'] <- log(0.99)
  
  lower[names(lower) == 'epsilon_tau'] <- 0
  # upper[names(upper) == 'epsilon_tau'] <- 1
  
  lower[names(lower) == 'logR_0k'] <- r0_lower
  
  lower[names(lower) == 'mort_k'] <- 0 ## duh
  upper[names(upper) == 'mort_k'] <- 1
  lower[names(lower) == 'logSDR'] <- log(0.0001)
  lower[names(lower) == 'omega_0ij'] <- log(0.0001)
  upper[names(upper) == 'omega_0ij'] <- log(0.999)

  lower[names(lower) == 'b'] <- 0
  upper[names(upper) == 'b'] <- 1
  
  lower[names(lower) == 'logq_f'] <- 0
  # upper[names(upper) == 'logq_f'] <- log(1)
  ## bounds on fsh slx ----
  
  ## brute force par locations; if mappy mirrors first ak fleets,
  ## assume that fleets 1 and 2 are now AK FIX mirror and AK TWL mirror
  # array(1:length(obj$par[names(obj$par) == "log_fsh_slx_pars"]), dim = c(4,2,2))
  
  
  if('fsh' %in% boundSlx){
    
    ## build fish bounds as if everything is there, then remove fixed fleets
    ## will have same dims as mappy and df$parms
    fsh_slx_map_lower <- fsh_slx_map_upper <- array(rep(NA,length(mappy$log_fsh_slx_pars)),
                         dim = dim(df$parms$log_fsh_slx_pars),
                         dimnames = dimnames(df$parms$log_fsh_slx_pars))
    
    ## if no fleets to fix, go with the normal slx bounds
    p1_logistic_idx <- c(1:2,15:16) #c(1:4,19:22)
    p2_logistic_idx <- p1_logistic_idx+df$nfleets_fish
    p1_norm_idx <- c(3,4,6,7,17,18,20,21)  #c(5,6,8,9,23,24,26,27)
    p2_norm_idx <- p1_norm_idx+df$nfleets_fish
    p1_gamma_idx <- c(5,19)#c(7,25)
    p2_gamma_idx <- p1_gamma_idx+df$nfleets_fish
    # #* fsh slx lower bounds ----
    
    # ## logistic p1 (a50)
    fsh_slx_map_lower[p1_logistic_idx] <-  log(35)
    ## logistic p2 (a95)
    fsh_slx_map_lower[p2_logistic_idx] <-  log(60)
    ## normal p1 (mean)
    fsh_slx_map_lower[p1_norm_idx] <-  log(15)
    ## normal p2 (sd)
    fsh_slx_map_lower[p2_norm_idx] <-  log(1)
    ## gamma shape (k*theta equals mean)
    fsh_slx_map_lower[p1_gamma_idx] <-log(15)
    ## gamma rate
    fsh_slx_map_lower[p2_gamma_idx] <-  log(2)
    # #* fsh slx upper bounds ----
    # ## logistic p1 (a50)
    fsh_slx_map_upper[p1_logistic_idx] <-  log(60)
    ## logistic p2 (a95)
    fsh_slx_map_upper[p2_logistic_idx] <-log(70)
    ## normal p1 (mean)
    fsh_slx_map_upper[p1_norm_idx] <-log(65)
    ## normal p2 (sd)
    fsh_slx_map_upper[p2_norm_idx] <- log(15)
    ## gamma shape (k*theta equals mean)
    fsh_slx_map_upper[p1_gamma_idx] <- log(35)
    ## gamma rate
    fsh_slx_map_upper[p2_gamma_idx] <- log(2)
    
    
    ##  check if slx was fixed at all
    if(length(grep("log_fsh_slx_pars", names(mappy))) != 0){
      ## identify positions of fixed fleet pars
      Nas <- which(is.na(mappy[[grep("log_fsh_slx_pars", names(mappy))]]))
      ## drop these from master bounds array
      fsh_slx_map_upper <- fsh_slx_map_upper[-Nas]
      fsh_slx_map_lower <- fsh_slx_map_lower[-Nas]
      
      ## only for sanity check
      # nfixedfleets <- length(Nas)/4
      # seeddim <- df$nfleets_fish-nfixedfleets
      # length(fsh_slx_map_upper) ==  length(mappy$log_fsh_slx_pars)-length(Nas)
      # array(fsh_slx_map_upper, dim = c(seeddim,2,1,2))
      # array(1:length(obj$par[names(obj$par) == "log_fsh_slx_pars"]), dim = c(6,2,2))
      ## identify which fleets were NA
      # rownames(mappy$log_fsh_slx_pars)
      # which(is.na(mappy[[grep("log_fsh_slx_pars", mappy)]]))
    }
    
    upper[names(upper) == "log_fsh_slx_pars"] <- fsh_slx_map_upper
    lower[names(lower) == "log_fsh_slx_pars"] <- fsh_slx_map_lower
    
  }
  # array(exp(upper[names(upper) == 'log_fsh_slx_pars']), dim = dim(df$parms$log_fsh_slx_pars))
  # array(exp(lower[names(lower) == 'log_fsh_slx_pars']), dim = dim(df$parms$log_fsh_slx_pars))
  # # array(lower[names(lower) == 'log_fsh_slx_pars'], dim = dim = c(7,2,2))
  
  # currently srv slx all logistic with a95, a50
  # nsurvsel = dim(df$parms$log_srv_slx_pars)[1]
  # array(1:length(  lower[names(lower) == 'log_srv_slx_pars']),
  # dim= c(dim(df$parms$log_srv_slx_pars)[1],2,2))
  if("srv" %in% boundSlx){

    
    ## first check if slx was fixed at all
    if(length(grep("log_srv_slx_pars", names(mappy)))  > 0){

      
      
      map_srvslx <- array(as.numeric(mappy$log_srv_slx_pars), 
                          dim = c(df$nfleets_surv+df$nfleets_acomp-4,2,max(df$srv_blks_size),2),
                          dimnames = dimnames(df$parms$log_srv_slx_pars))

      lwr.temp <- upr.temp <- map_srvslx ## now there are only values where we need to fill them
      ## identify which fleets were fixed by just looking at first block
      # mappy_srvslx <- mappy[[grep("log_srv_slx_pars", names(mappy))]]
      # Nas <- which(is.na(mappy_srvslx)[1:16])
      #   
      #   # which(is.na(mappy[[grep("log_srv_slx_pars", names(mappy))]][1:32]))
      # nfixedfleets <- length(Nas)/2
      # keptflts <- dimnames(df$parms$log_srv_slx_pars)[[1]][-Nas] ## return non-fixed fltnames
      # 
      # lwr.temp <- upr.temp <- array(1:length(obj$par[names(obj$par) == "log_srv_slx_pars"]), 
      #       dim = c(8-nfixedfleets,2,max(df$srv_blks_size),2),
      #       dimnames = list(keptflts,c('p1','p2'),
      #                       c('block',1:max(df$srv_blks_size)),c('Fem','Mal')))
      #
      lwr.temp['AK_VAST_W',"p1",1:df$srv_blks_size[,'AK_VAST_W'],] <- log(30); 
      upr.temp['AK_VAST_W',"p1",1:df$srv_blks_size[,'AK_VAST_W'],] <- 3.68887945; 
      lwr.temp['AK_VAST_W',"p2",1:df$srv_blks_size[,'AK_VAST_W'],] <- log(40); 
      upr.temp['AK_VAST_W',"p2",1:df$srv_blks_size[,'AK_VAST_W'],] <- 4.24849524
      
      lwr.temp['AK_VAST_E',"p1",1:df$srv_blks_size[,'AK_VAST_E'],] <- 3.401197;
      upr.temp['AK_VAST_E',"p1",1:df$srv_blks_size[,'AK_VAST_E'],] <- 3.68887945
      lwr.temp['AK_VAST_E',"p2",1:df$srv_blks_size[,'AK_VAST_E'],] <- 4.007333; 
      upr.temp['AK_VAST_E',"p2",1:df$srv_blks_size[,'AK_VAST_E'],] <- 4.24849524
      
      lwr.temp['BC_EARLY',"p1",1:df$srv_blks_size[,'BC_EARLY'],] <- 3.401197; 
      upr.temp['BC_EARLY',"p1",1:df$srv_blks_size[,'BC_EARLY'],] <- 3.68887945
      lwr.temp['BC_EARLY',"p2",1:df$srv_blks_size[,'BC_EARLY'],] <- 4.007333; 
      upr.temp['BC_EARLY',"p2",1:df$srv_blks_size[,'BC_EARLY'],] <- 4.24849524
      
      lwr.temp['BC_VAST',"p1",1:df$srv_blks_size[,'BC_VAST'],] <- 3.401197; 
      upr.temp['BC_VAST',"p1",1:df$srv_blks_size[,'BC_VAST'],] <-3.68887945
      lwr.temp['BC_VAST',"p2",1:df$srv_blks_size[,'BC_VAST'],] <- 4.007333; 
      upr.temp['BC_VAST',"p2",1:df$srv_blks_size[,'BC_VAST'],] <- 4.24849524
      
      lwr.temp['WC_VAST',"p1",1:df$srv_blks_size[,'WC_VAST'],] <- log(30);
      upr.temp['WC_VAST',"p1",1:df$srv_blks_size[,'WC_VAST'],] <- log(75)
      lwr.temp['WC_VAST',"p2",1:df$srv_blks_size[,'WC_VAST'],'Mal'] <- log(50);  
      lwr.temp['WC_VAST',"p2",1:df$srv_blks_size[,'WC_VAST'],'Fem'] <- log(45);
      upr.temp['WC_VAST',"p2",1:df$srv_blks_size[,'WC_VAST'],] <- log(75)
      
      lower[names(lower) == 'log_srv_slx_pars'] <- lwr.temp[!is.na(lwr.temp)]
      upper[names(upper) == 'log_srv_slx_pars'] <- upr.temp[!is.na(upr.temp)]
      ## general bounds on p1 and p2 for all fleets
      # seeddim <-  df$nfleets_surv+df$nfleets_acomp-4-nfixedfleets
      # lower[names(lower) == 'log_srv_slx_pars'][c(1:seeddim,(2*seeddim+1):(3*seeddim))] <- log(30) ## p1 
      # lower[names(lower) == 'log_srv_slx_pars'][c((seeddim+1):(2*seeddim),(3*seeddim+1):(4*seeddim))] <- log(30) ## p2
      # 
      # upper[names(upper) == 'log_srv_slx_pars'][c(1:seeddim,(2*seeddim+1):(3*seeddim))] <- log(55) ## p1 3.91202301
      # upper[names(upper) == 'log_srv_slx_pars'][c((seeddim+1):(2*seeddim),(3*seeddim+1):(4*seeddim))] <- log(80) ## p2 3.91202301
      
      ## just rational bounds
      # lower[names(lower) == 'log_srv_slx_pars'] <- 0
      # upper[names(upper) == 'log_srv_slx_pars'] <- log(80)

      
    } else if(length(grep("log_srv_slx_pars", names(mappy)))  == 0){
      
      seeddim = df$nfleets_surv+df$nfleets_acomp-4
      lower[names(lower) == 'log_srv_slx_pars'][c(1:seeddim,(2*seeddim+1):(3*seeddim))] <- log(30) ## p1
      lower[names(lower) == 'log_srv_slx_pars'][c((seeddim+1):(2*seeddim),(3*seeddim+1):(4*seeddim))] <- log(70) ## p2

      upper[names(upper) == 'log_srv_slx_pars'][c(1:seeddim,(2*seeddim+1):(3*seeddim))] <- log(30) ## p1
      upper[names(upper) == 'log_srv_slx_pars'][c((seeddim+1):(2*seeddim),(3*seeddim+1):(4*seeddim))] <- log(70) ## p2
      
      
      
      ## custom bounds
      # upper[names(upper) == 'log_srv_slx_pars'][c(3,19)] <- log(75) ## BCEARLY P1
      # upper[names(upper) == 'log_srv_slx_pars'][c(11,27)] <- log(75) ## BCEARLY P2
      # lower[names(lower) == 'log_srv_slx_pars'][c(11,27)] <- log(10) ## BCEARLY P2
      
      # upper[names(upper) == 'log_srv_slx_pars'][c(2,18)] <- log(80) ## AKVASTE P1
      # upper[names(upper) == 'log_srv_slx_pars'][c(3,19)] <- log(80) ## BCEARLY P1
      
      # upper[names(upper) == 'log_srv_slx_pars'][10] <- log(85) ## AKVASTE P2 fem only
      # upper[names(upper) == 'log_srv_slx_pars'][20] <- log(80) ## BCVAST P1 MALE only
      # upper[names(upper) == 'log_srv_slx_pars'][c(11,27)] <- log(85) ## BCEARLY P2
      
      # upper[names(upper) == 'log_srv_slx_pars'][c(11,27)] <- log(85) ## BCEARLY P2
      
      # lower[names(lower) == 'log_srv_slx_pars'] <- 0
      # upper[names(upper) == 'log_srv_slx_pars']<- log(80)

    }

  } ## end srv in boundslx

  return(list("upper"=upper, "lower"=lower))
}

