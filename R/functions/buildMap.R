buildMap <- function(toFix = c("omega_0ij","epsilon_tau", "mort_k", "logh_k"),
                     fixFlt = c("BC_LL","BC_TRAP","BC_TWL")){
  
  mappy <- list(); idx <- 1
  
  # for(i in 1:length(df$parms)){
  for(i in seq_along(toFix)){
    parloc <- which(names(df$parms) == toFix[i])
    if(toFix[i] != 'log_fsh_slx_pars'  & toFix[i] != 'log_srv_slx_pars'){
      mappy[[idx]] <- factor(rep(NA, length(df$parms[[parloc]])))
      names(mappy)[idx] <- toFix[i]
      idx = idx+1
      ## end non slx pars
    } else{
      if(toFix[i]  == 'log_fsh_slx_pars'){
        
          fsh_slx_map <- array(1:length(df$parms$log_fsh_slx_pars),
                               dim = dim(df$parms$log_fsh_slx_pars),
                               dimnames = dimnames(df$parms$log_fsh_slx_pars))
          if('all_fsh' %in% fixFlt) fsh_slx_map[1:length(df$parms$log_fsh_slx_pars)]  <- factor(NA)
          for(flt in fixFlt){
            fsh_slx_map[row.names(fsh_slx_map) == flt,1:2,,1:2] <- factor(NA)
          }
          ## also fix fleets which have no time block population
          for(i in seq_along(dimnames(df$parms$log_fsh_slx_pars)[[1]])){
            if(df$fsh_blks_size[i] != max(df$fsh_blks_size)){
              fsh_slx_map[row.names(fsh_slx_map) == dimnames(df$parms$log_fsh_slx_pars)[[1]][i],
                          1:2,(df$fsh_blks_size[i]+1):max(df$fsh_blks_size),1:2] <- factor(NA)
            } 
            
          }
          
          mappy[[idx]] <- factor(fsh_slx_map)
          names(mappy)[idx] <- 'log_fsh_slx_pars'
          
          idx = idx+1
        } else if( toFix[i]  == 'log_srv_slx_pars' ){
          # stop("fx not ready to automate fixing survey slx")
          srv_slx_map <- array(1:length(df$parms$log_srv_slx_pars),
                               dim = dim(df$parms$log_srv_slx_pars),
                               dimnames = dimnames(df$parms$log_srv_slx_pars))
          for(flt in fixFlt){
            srv_slx_map[row.names(srv_slx_map) == flt,1:2,,1:2] <- factor(NA) ## fix all designated fleets
          }
          ## also fix fleets which have no time block population
          for(i in seq_along(dimnames(df$parms$log_srv_slx_pars)[[1]])){
            if(df$srv_blks_size[i] != max(df$srv_blks_size)){
              srv_slx_map[row.names(srv_slx_map) == dimnames(df$parms$log_srv_slx_pars)[[1]][i],
                          1:2,(df$srv_blks_size[i]+1):max(df$srv_blks_size),1:2] <- factor(NA)
            } 
            
          }
          
          if('all_srv' %in% fixFlt) srv_slx_map[1:length(df$parms$log_srv_slx_pars)]  <- factor(NA)
          
          
          mappy[[idx]] <- factor(srv_slx_map)
          names(mappy)[idx] <-'log_srv_slx_pars'
          
          idx = idx+1
          
        }
      } ## end  slx pars
    } ## end i in to fix

  return(mappy)
}

# df$parms$logR_0k = rep(25,4)


# omega_0ij_map[1,] <- df$parms$omega_0ij[1,] ## estimate to/from C1 only
# 
# ## mirror selex in AK E/W 
# ## if you want to MIRROR selex, fill a value in the specific location which is identical for each fleet

# dimnames(fsh_slx_map)[[1]] <- df$fltnames_fish
# 
# fsh_slx_map[c(1,2)] <- 1 ## mirror p1 for females, W
# fsh_slx_map[c(19,20)] <- 2 ## mirror p1 for males, W
# fsh_slx_map[c(3,4)] <- 3 ## mirror p1 for females, E
# fsh_slx_map[c(21,22)] <- 4 ## mirror p1 for males, E
# fsh_slx_map[c(10,11)] <- 5## mirror p2 for females, W
# fsh_slx_map[c(28,29)] <- 6## mirror p2 for males, W
# fsh_slx_map[c(12,13)] <- 7## mirror p2 for females, E
# fsh_slx_map[c(30,31)] <- 8## mirror p2 for males, E
## fix BC selex (for use with -1 slx)