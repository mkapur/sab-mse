load_data_OM <- function(nspace = 6, 
                         nstocks = 4,
                         myear = 2019,
                         move = TRUE, 
                         LBins = 81,
                         logSDR = 1.4, 
                         bfuture = 0.5,
                         yr_future  = 0,
                         b = 0.5
){
  
  #' @nspace = Spatial areas
  #' @nstocks = demographic stocks
  #' @myear = Final year that contains empirical data
  #' @move = do you want to enable movement among spaces?
  #' @logSDR Standard deviation of recruitment
  #' @bfuture recruitment bias adjustment in the future - scalar
  #' @yr_future Create dummy data for future years
  
  ## ERROR TRAPS 
  if(nspace  <4 & nstocks >3) stop("mismatch between nspace and nstocks")
  if(nspace == 1 & move) stop("can't have movement with fully pooled model")

  years <- 1960:(myear+yr_future)
  nyear <- length(years)
  tEnd <- length(years)
  age <- 0:70 
  nage <- length(age)
  
  # Maturity ----
  load(here("input","input_data","OM_maturity_ak.rdata")) ## ak is age, stock
  # movement ----
  
  ## placeholder for X_ija -this will need to get converted from length
  if(move == FALSE){
    X_ijas <- array(0, dim = c(nspace,nspace,nage,2))
    for(s in 1:2){
      for(a in 1:dim(X_ijas)[[3]]){
        diag(X_ijas[,,a,s]) <- 1
      }## end age
    } ## end sex
  }## end move == FALSE
  if(move == TRUE) {
    load(here("input","input_data","X_ijas.rdata"))
  } ## end move == TRUE
  
  omega_ais <- array(0, dim = c(nage,nspace,2))## eigenvector for stable spatial distribution at age
  for(s in 1:2){
    for(a in 1:nage){
      # omega_ais[a,,s] <- round(eigen(X_ijas[,,a,s])$values,3)
      omega_ais[a,,s] <- round(eigen(X_ijas[,,a,s])$values/
                                 sum(eigen(X_ijas[,,a,s])$values),3)
      # omega_ais[a,,s][which(omega_ais[a,,s] < 0)] <- 0.05
      if(all(X_ijas[,,a,s] == 0)){ omega_ais[a,,s] <- 1/nspace} ## since we have
      ## no movement values for sub A5, just assume even dist
      ## I forced this to be 166 vs 167 because it was causing too many individuals
    }
  }
  # omega_0ij = matrix(0, nrow = nspace, ncol = nspace)
  omega_0ij = X_ijas[,,6,1]
  
  # Weight at length ----
  load(here("input","input_data","OM_wtatlen_kab.rdata")) ## a and be are pars of al^b
  
  # Growth ----
  load(here("input","input_data","OM_growthPars.rdata")) 
  
  ## Mortality ----
  load(here('input','input_data','M_k.rdata'))
  
  ## Fleet [names and nfleets] ----
  ## build fleets
  ## makes the master flag_fleets matrix
  ## and attendant indices for subsetting
  fltnames <- read.table(here("input","input_data","OM_fleetnames.txt"), header = TRUE) ## this is like flag_fleets
  fltnames_fish <- fltnames$NAME[fltnames$COMM]
  fltnames_surv <- fltnames$NAME[fltnames$SURV]
  fltnames_acomp <- fltnames$NAME[fltnames$ACOMP]
  fltnames_lcomp <- fltnames$NAME[fltnames$LCOMP]
  
  nfleets_fish <- length(fltnames$NAME[fltnames$COMM])
  nfleets_surv <- length(fltnames$NAME[fltnames$SURV])
  nfleets_acomp <- length(fltnames$NAME[fltnames$ACOMP])
  nfleets_lcomp <- length(fltnames$NAME[fltnames$LCOMP])
  
  
  
  # Catch ----
  catch <- round(read.csv(here("input","input_data","OM_catch.csv")),1)
  catch[is.na(catch)] <- -1.0
  catch_yf_error = array(0.1, dim = dim(catch))
  
  ## Discard ----
  load(here("input","input_data","OM_discard.csv")) ## loads as omdis
  
  ## Selex ----
  load(here('input','input_data',"OM_fish_selex_yafs.rdata"))
  load(here('input','input_data',"OM_surv_selex_yafs.rdata"))
  
  ## time blocks by fleet
  srv_blks <- matrix(0, nrow = tEnd, ncol = nfleets_surv+(nfleets_acomp-4))
  fsh_blks <- matrix(0, nrow = tEnd, ncol = nfleets_fish)
  # srv_blks <- matrix(tEnd-1, ncol = 1, nrow = nfleets_surv+(nfleets_acomp-4))
  # fsh_blks <-  matrix(tEnd-1, ncol = 1, nrow = nfleets_fish)
  
  selType_fish <- as.numeric(fltnames$SELTYPE[fltnames$COMM])-1
  ## note that the first two acomp fleets are already inside seltype fish
  ## only first ONE if AK fix not aggregated
  selType_surv <- as.numeric(c(fltnames$SELTYPE[fltnames$SURV],fltnames$SELTYPE[fltnames$ACOMP][c(2,4,5)]))-1
  selShape_fish <- c(0,2,2,2,3,2,2) ## 0 and 1 logistic, 2 dome normal, 3 dome gamma
  selShape_surv <- c(rep(0,nfleets_surv+(nfleets_acomp-4))) ## 0 and 1 logistic, 2 dome normal, 3 dome gamma
  if(length(selType_surv) != length(selShape_surv)) stop("seltype surv length doesn't match selshape surv")
  # Survey ----
  # survey <- read.csv(here("input","input_data",'OM_indices.csv'))
  # survey[is.na(survey)] <- -1.0## flag for numeric TMB checks
  # survey[,"BC_EARLY"] <-  survey[,"BC_EARLY"] +0.0111 ## flag for numeric TMB checks
  # survey <- round(survey,  1)
  # survey_err <- read.csv(here("input","input_data",'OM_indices_sigma.csv'))
  
  # survey <- read.csv(here("input","input_data",'OM_indices_BaseQ=GOA_Late.csv'))
  # survey[is.na(survey)] <- -1.0## flag for numeric TMB checks
  # survey[,"BC_EARLY"] <-  survey[,"BC_EARLY"] +0.0111 ## flag for numeric TMB checks
  # survey <- round(survey,  1)
  # survey_err <- read.csv(here("input","input_data",'OM_indices_sigma_BaseQ=GOA_Late.csv'))
  
  survey <- read.csv(here("input","input_data",'OM_indices_BaseQ=WCGBTS.csv'))
  survey[is.na(survey)] <- -1.0## flag for numeric TMB checks
  # survey[,"BC_EARLY"] <-  survey[,"BC_EARLY"] +0.0111 ## flag for numeric TMB checks
  # survey <- round(survey,  1)
  survey_err <- read.csv(here("input","input_data",'OM_indices_sigma_BaseQ=WCGBTS.csv'))
  
  
  ## Comps ----
  ## Len comps [these are arrays by fleet]
  load(here("input","input_data",'OM_lencomps_female.rdata'))
  load(here("input","input_data",'OM_lencomps_male.rdata'))
  
  ## Age comps. Note that AK is not sex-specific
  ## need to make these -1 for NA years
  # load(here("input","input_data",'OM_agecomps_female.rdata'))
  # load(here("input","input_data",'OM_agecomps_male.rdata'))
  load(here("input","input_data",'OM_agecomps_yafs.rdata'))
  ## Aging Error ----
  ## M X age
  load(here("input","input_data",'ageerr_ExpAge.rdata'))
  load(here("input","input_data",'ageerr_SD.rdata'))
  
  ## Phi objects ----
  ## setup phi (spatial matching matrix) depending on spatial setup
  ## MAKE A NSPACE == 3 AND 1 OPTION FOR COMBINING
  spmat <- data.frame(subarea = c('A1',"A3","B3","B2","C2","C1"),
                      stock = c("R4","R3","R3","R2","R2","R1"),
                      mgmt = c("AK","AK", rep("BC",2), rep("CC",2)))
  if(nspace == 6){ ## OM
    
    ## phi_survy
    phi_if_surv <- matrix(0, nrow = nfleets_surv, ncol = nspace)
    rownames(phi_if_surv) <- names(survey)
    colnames(phi_if_surv) <- rev(spmat$subarea)
    # phi_if_surv[1,6] <-  phi_if_surv[2,5] <-  
    #   phi_if_surv[3:4,3:4]<-  phi_if_surv[5,1:2] <- 1
    for(i in 1:nrow(phi_if_surv)){
      reg = substr(rownames(phi_if_surv)[i],1,2)
      if(reg == 'AK') {
        phi_if_surv[i,5:6] <- 1
      } else  if(reg == 'BC'){
        phi_if_surv[i,3:4] <- 1
      }else{
        phi_if_surv[i,1:2] <- 1
      }
    }
    
    phi_if_acomp <- matrix(0, nrow = nfleets_acomp, ncol = nspace)
    rownames(phi_if_acomp) <- fltnames_acomp
    colnames(phi_if_acomp) <- rev(spmat$subarea)
    for(i in 1:nrow(phi_if_acomp)){
      reg = substr(rownames(phi_if_acomp)[i],1,2)
      if(reg == 'AK') {
        phi_if_acomp[i,5:6] <- 1
      } else  if(reg == 'BC'){
        phi_if_acomp[i,3:4] <- 1
      }else{
        phi_if_acomp[i,1:2] <- 1
      }
    }
    
    # phi_if_acomp[1,6] <-  phi_if_acomp[2:3,5] <-  
      # phi_if_acomp[4:6,3:4]<-  phi_if_acomp[7:8,1:2] <- 1
    
    ## indicates the position of acomp fleet
    phi_ff_acomp <- matrix(-1, nrow = nfleets_acomp, ncol = 5) 
    rownames(phi_ff_acomp) <- fltnames_acomp
    colnames(phi_ff_acomp) <- c('fsh_slx_pos','srv_slx_pos',"nsamp_pos","commacomp_pos","survacomp_pos")
    ## MANUAL UPDATE WITH COLNAMES [seltypes are in same order as fltnames]
    phi_ff_acomp[which(rownames(phi_ff_acomp) %in% paste(fltnames_fish)),1] <- 
      which(grepl(paste(rownames(phi_ff_acomp), collapse = "|"), paste(fltnames_fish)))-1
    # phi_ff_acomp[,1] <- c(0,1,-1,5,-1,-1,7,8) ## position in fishery selex (-1 means not applicable)
    phi_ff_acomp[,2] <- c(-1,5,-1,6,7,-1,-1) ## Pos in survey
    phi_ff_acomp[,3] <- c(5:11) ## ordering for nsamp
    phi_ff_acomp[which(rownames(phi_ff_acomp) %in% paste(fltnames_fish)),4] <- 0:3
    # phi_ff_acomp[,4] <- c(0,1,-1,2,-1,-1,3,4) ## ordering for comm comps
    phi_ff_acomp[,5] <- c(-1,0,-1,1,2,-1,-1) ## ordering for surv comps (only 3)

    ## phi_fish
    phi_if_fish <- matrix(0, nrow = nfleets_fish, ncol = nspace) ## placeholder for fishing fleets
    rownames(phi_if_fish) <- names(catch)[2:ncol(catch)]
    colnames(phi_if_fish) <-  rev(spmat$subarea)
    ## works whether you have 7 or 9 fleets
    for(i in 1:nrow(phi_if_fish)){
      reg = substr(rownames(phi_if_fish)[i],1,2)
      if(reg == 'AK') {
        phi_if_fish[i,5:6] <- 1
      } else  if(reg == 'BC'){
        phi_if_fish[i,3:4] <- 1
      }else{
        phi_if_fish[i,1:2] <- 1
      }
    }
    # phi_if_fish[c(1,3),6] <- phi_if_fish[c(2,4),5] <-   phi_if_fish[5:7,3:4] <-  
    #   phi_if_fish[c(8,9),1:2] <-  1
    
    ## phi_im
    phi_im <- matrix(0, ncol = 3, nrow = nspace)
    colnames(phi_im) <- rev(unique(spmat$mgmt))
    rownames(phi_im) <- rev(spmat$subarea)
    phi_im[1:2,1] <- phi_im[3:4,2] <- phi_im[5:6,3] <- 1
    
    ## phi_ik
    phi_ki <-  matrix(0, ncol = nspace, nrow = nstocks) ## nesting of subareas within stocks, for recruitment purposes
    rownames(phi_ki) <- rev(unique(spmat$stock))
    colnames(phi_ki) <- rev(spmat$subarea)
    
    phi_ki[1,1] <-  phi_ki[2,2:3] <-  phi_ki[3,4:5]<-  phi_ki[4,6]  <- 1
    phi_ik2 <- matrix(apply(phi_ki,2, function(x)which(x == 1))-1) ## a vector for par subsetting, the columns are subareas
    
    ## phi_ij [eq 6]
    phi_ij <-  matrix(1, ncol = nspace, nrow = nspace) ## 0 indicates  subareas comprise THE SAME stock
    rownames(phi_ij) = colnames(phi_ij) = rev(spmat$subarea)
    diag(phi_ij) <- phi_ij[2:3,2:3] <- phi_ij[4:5,4:5] <- 0
    
    
    ## phi_fm
    phi_fm <- matrix(0, nrow = nfleets_fish, ncol = 3)
    rownames(phi_fm) = names(catch)[2:ncol(catch)]
    colnames(phi_fm) = rev(unique(spmat$mgmt))
    for(i in 1:nrow(phi_fm)){
      reg = substr(rownames(phi_fm)[i],1,2)
      if(reg == 'AK') {
        phi_fm[i,3] <- 1
      } else  if(reg == 'BC'){
        phi_fm[i,2] <- 1
      }else{
        phi_fm[i,1] <- 1
      }
    }
    ## same as above but for comps (mix of fisheries & surveys)
    phi_fm_acomp <- matrix(0, nrow = nfleets_acomp, ncol = 3)
    rownames(phi_fm_acomp) = fltnames_acomp
    colnames(phi_fm_acomp) = rev(unique(spmat$mgmt))
    for(i in 1:nrow(phi_fm_acomp)){
      reg = substr(rownames(phi_fm_acomp)[i],1,2)
      if(reg == 'AK') {
        phi_fm_acomp[i,3] <- 1
      } else  if(reg == 'BC'){
        phi_fm_acomp[i,2] <- 1
      }else{
        phi_fm_acomp[i,1] <- 1
      }
    }
    # phi_fm_acomp[1:3,3] <- phi_fm_acomp[4:6,2]  <- phi_fm_acomp[7:8,1]  <- 1
    phi_fm_acomp2 <- matrix(apply(phi_fm_acomp,1, function(x)which(x == 1))-1) ## a vector for par subsetting, the columns are survey fleets
    
    acomp_flt_type <- matrix(1, ncol = nfleets_acomp) ## 0 is commercial, 1 is survey
    colnames(acomp_flt_type) <- fltnames_acomp
    ## match colnames to fltnames which are commercial, assign those as 1
    acomp_flt_type[which(fltnames$COMM[which(grepl(paste(colnames(acomp_flt_type), collapse = "|"),fltnames$NAME))])] <- 0
    # acomp_flt_type[fltnames$NAME == rownames(acomp_flt_type) ] <- 1

    phi_lcomp_fm <- matrix(0, nrow = nfleets_lcomp, ncol = 3)
    rownames(phi_lcomp_fm) = fltnames_lcomp
    colnames(phi_lcomp_fm) = rev(unique(spmat$mgmt))
    # phi_lcomp_fm[1:6,3] <- phi_lcomp_fm[7:9,2]  <- phi_lcomp_fm[10,1]  <- 1
    for(i in 1:nrow(phi_lcomp_fm)){
      reg = substr(rownames(phi_lcomp_fm)[i],1,2)
      if(reg == 'AK') {
        phi_lcomp_fm[i,3] <- 1
      } else  if(reg == 'BC'){
        phi_lcomp_fm[i,2] <- 1
      }else{
        phi_lcomp_fm[i,1] <- 1
      }
    }
    ## tau_ki
    tau_ki <-  matrix(0, ncol = nspace, nrow = nstocks) ## nesting of subareas within stocks, for recruitment purposes
    rownames(tau_ki) <- rev(unique(spmat$stock))
    colnames(tau_ki) <- rev(spmat$subarea)
    tau_ki[1,1] <-   tau_ki[4,6]  <- 1 ## 100% of recruitment in stock
    tau_ki[2,2:3] <-  tau_ki[3,5:4] <-  c(0.75,0.25) ## A2 and C1 are larger
  } else {
    # phi_if_surv <- matrix(rbinom(nfleets_surv*nspace,1,0.5), byrow = TRUE, nrow = nfleets_surv, ncol = nspace) ## placeholder for alternative spatial stratifications
    # phi_if_fish <- matrix(c(0,1,1,1,1,0), nrow = nfleets_fish, ncol = nspace)  ## placeholder for fishing fleets
    # phi_ki <- matrix(c(1,0,0,1), byrow = TRUE, nrow = nstocks, ncol = nspace) ## placeholder for alternative spatial stratifications
    # phi_ik2 <- apply(phi_ki,2, function(x)which(x == 1))-1 ## a vector for par subsetting, the columns are subareas
    # 
    # # phi_fm <- matrix(0, nrow = nspace, ncol = 2)
    # # phi_fm[1,1] <- phi_fm[2,2] <- 1
    # phi_fm <- matrix(0, nrow = nfleets_fish, ncol = 3)
    # phi_fm[1:2,1] <- phi_fm[3,3]  <- 1
    # ## autogenerate stock-distinction matrix
    # phi_ij <-  matrix(NA, byrow = TRUE, ncol = nspace, nrow = nspace)
    # for(i in 1:nspace){
    #   for(j in 1:nspace){
    #     phi_ij[i,j] = ifelse(phi_ik2[i] == phi_ik2[j],0,1)
    #   }
    # }
    # tau_ki <- matrix(c(0.25,0.75,0.9,0.1), nrow = nstocks, byrow = TRUE, ncol = nspace) ## placeholder for alternative spatial stratifications
  }
  
  # b <- matrix(NA, nyear)
  # Yr <- 1946:max(years)
  # # Parameters 
  # yb_1 <- 1965 #_last_early_yr_nobias_adj_in_MPD
  # yb_2 <- 1971 #_first_yr_fullbias_adj_in_MPD
  # yb_3 <- 2016 #_last_yr_fullbias_adj_in_MPD
  # yb_4 <- max(years) #_first_recent_yr_nobias_adj_in_MPD
  # b_max <- 0.87 #_max_bias_adj_in_MPD
  # 
  # b[1] <- 0
  # for(j in 2:length(Yr)){
  #   
  #   if (Yr[j] <= yb_1){
  #     b[j] = 0}
  #   
  #   if(Yr[j] > yb_1 & Yr[j]< yb_2){
  #     b[j] = b_max*((Yr[j]-yb_1)/(yb_2-yb_1));
  #   }
  #   
  #   if(Yr[j] >= yb_2 & Yr[j] <= yb_3){
  #     b[j] = b_max}
  #   
  #   if(Yr[j] > yb_3 & Yr[j] < yb_4){
  #     b[j] = b_max*(1-(yb_3-Yr[j])/(yb_4-yb_3))
  #   }
  #   
  #   if(Yr[j] >= yb_4){
  #     b[j] = 0
  #   }
  #   # if (b[j]<b[j-1]){
  #   #   stop('why')
  #   # }
  # }  
  ## log_fsh_slx_pars ----
  # fsh_blks_size is a 1 x nfleets_surv ivector which indicates the number of timeblocks applicable
  # to each fleet.
  fsh_blks_size <- matrix(1, nrow = 1, ncol = nfleets_fish)
  colnames(fsh_blks_size) <- c( as.character(fltnames_fish))
  fsh_blks_size[,'WC_FIX'] <- 4
  fsh_blks_size[,'WC_TWL'] <- 4
  fsh_blks_size[,'AK_FIX'] <- 2
  # fsh_blks is an h x nfleets_fish imatrix with the MAX year of a given timeblock.
  # it will be a ragged array bc some fleets have fewer blocks.
  fsh_blks <- matrix(2019, nrow = max(fsh_blks_size), 
                     ncol = nfleets_fish)
  colnames(fsh_blks) <- c( as.character(fltnames_fish))
  fsh_blks[1:fsh_blks_size[,'WC_FIX'],'WC_FIX' ] <- c(1997,2003,2010,2019)
  fsh_blks[1:fsh_blks_size[,'WC_TWL'],'WC_TWL' ] <- c(1982,2003,2010,2019)
  fsh_blks[1:fsh_blks_size[,'AK_FIX'],'AK_FIX' ] <- c(1995,2019)
  
  ## fill all blocks with start pars
  log_fsh_slx_pars = array(NA, dim = c(nfleets_fish,2, max(fsh_blks_size),2), 
                           dimnames = list(c(paste(fltnames_fish)),
                                           c("p1","p2"),
                                           c(paste0('block',1: max(fsh_blks_size))),
                                           c('Fem','Mal')))
  
  ## ak fix is logistic
  log_fsh_slx_pars["AK_FIX","p1",1:fsh_blks_size[,"AK_FIX"], c('Fem','Mal')] <- log(50)
  log_fsh_slx_pars["AK_FIX","p2",1:fsh_blks_size[,"AK_FIX"], c('Fem','Mal')] <-  log(67)
  ## ak twl is dome normal
  log_fsh_slx_pars["AK_TWL","p1",1:fsh_blks_size[,"AK_TWL"], c('Fem','Mal')] <- log(45)
  log_fsh_slx_pars["AK_TWL","p2",1:fsh_blks_size[,"AK_TWL"], c('Fem','Mal')] <-  log(10)
  ## have custom start pars for bc fleets, which are dome normal and gamma
  log_fsh_slx_pars["BC_LL","p1",1:fsh_blks_size[,"BC_LL"], c('Fem','Mal')] <- log( 63.6959)
  log_fsh_slx_pars["BC_LL","p2",1:fsh_blks_size[,"BC_LL"], c('Fem','Mal')] <-  log(3.09715)
  log_fsh_slx_pars["BC_TRAP","p1",1:fsh_blks_size[,"BC_TRAP"], c('Fem','Mal')] <- log( 62.8329)
  log_fsh_slx_pars["BC_TRAP","p2",1:fsh_blks_size[,"BC_TRAP"], c('Fem','Mal')] <-  log(7.04483)
  log_fsh_slx_pars["BC_TWL","p1",1:fsh_blks_size[,"BC_TWL"], c('Fem','Mal')] <- log( 33.8898)
  log_fsh_slx_pars["BC_TWL","p2",1:fsh_blks_size[,"BC_TWL"], c('Fem','Mal')] <-  log(1.41494)
  ## Just copy dome normal values for wc
  log_fsh_slx_pars["WC_FIX","p1",1:fsh_blks_size[,"WC_FIX"], c('Fem','Mal')] <- log( 63.6959)
  log_fsh_slx_pars["WC_FIX","p2",1:fsh_blks_size[,"WC_FIX"], c('Fem','Mal')] <-  log(3.09715)
  log_fsh_slx_pars["WC_TWL","p1",1:fsh_blks_size[,"WC_TWL"], c('Fem','Mal')] <- log( 62.8329)
  log_fsh_slx_pars["WC_TWL","p2",1:fsh_blks_size[,"WC_TWL"], c('Fem','Mal')] <-  log(7.04483)

  ## log_srv_slx_pars ----
  
  # srv_blks_size is a 1 x nfleets_surv ivector which indicates the number of timeblocks applicable
  # to each fleet.
  srv_blks_size <- matrix(1, nrow = 1, ncol = nfleets_surv+nfleets_acomp-4)
  colnames(srv_blks_size) <- c( as.character(fltnames_surv), as.character(fltnames_acomp[c(2,4,5)]))
  srv_blks_size[,'WC_VAST'] <- 3
  # srv_blks is an h x nfleets_surv imatrix with the MAX year of a given timeblock.
  # it will be a ragged array bc some fleets have fewer blocks.
  srv_blks <- matrix(2019, nrow = max(srv_blks_size), 
                     ncol = nfleets_surv+nfleets_acomp-4)
  colnames(srv_blks) <- c( as.character(fltnames_surv), as.character(fltnames_acomp[c(2,4,5)]))
  srv_blks[1:srv_blks_size[,'WC_VAST'],'WC_VAST' ] <- c(1995,2010,2019)


  ## all of these are currently logistic with l/a50, and a delta
  log_srv_slx_pars =  array(0, dim = c( nfleets_surv+(nfleets_acomp-4),2,1,2),   
                            dimnames = list(c(paste(fltnames_surv),paste(fltnames_acomp[c(2,4,5)])),
                                           c("p1","p2"),
                                           c(paste0('block',1)),
                                           c('Fem','Mal')))
  ## males and females
  log_srv_slx_pars['AK_VAST_W','p1',,] <- c(35.36234,42.69320)
  log_srv_slx_pars['AK_VAST_W','p2',,] <- c(45.31228,44.42085)
  log_srv_slx_pars['AK_VAST_E','p1',,] <- 50
  log_srv_slx_pars['AK_VAST_E','p2',,] <- 65
  log_srv_slx_pars['BC_EARLY','p1',,] <- 75.62026
  log_srv_slx_pars['BC_EARLY','p2',,] <- 69.27150
  log_srv_slx_pars['BC_VAST','p1',,] <-  c(38.31500,50)
  log_srv_slx_pars['BC_VAST','p2',,] <- c(43.05036,69.27150)
  log_srv_slx_pars['WC_VAST','p1',,1:2] <- c(43.04513,39.59406)
  log_srv_slx_pars['WC_VAST','p2',,1:2] <- c(55.41398,41.78538)
  
  ## comps
  log_srv_slx_pars['AK_GOA_SURV','p1',,] <- 50
  log_srv_slx_pars['AK_GOA_SURV','p2',,] <- 65
  
  log_srv_slx_pars['BC_StRS',1,1,1] <- 64.2127 - 12.9197
  log_srv_slx_pars['BC_StRS',1,1,2] <-  64.2127 - 12.9197
  log_srv_slx_pars['BC_StRS',2,1,1] <- 70
  log_srv_slx_pars['BC_StRS',2,1,2] <- 70
  log_srv_slx_pars['BC_SS',1,1,1] <- 54.1045-4.55724
  log_srv_slx_pars['BC_SS',1,1,2] <- 54.1045-4.55724
  log_srv_slx_pars['BC_SS',2,1,1] <- 65
  log_srv_slx_pars['BC_SS',2,1,2] <- 65

  
  log_srv_slx_pars = log(log_srv_slx_pars)
  
  mort_k <- c(0.2,0.15,0.05,0.1)

  ## Parms List ----
  ## things that will get estimated later on, everthing else is FIXED
  ## note that these go from AK to WC
  parms <- list(
    logh_k = log(c(0.7,0.88,0.7,0.7)),

    logR_0k = rep(log(8*10e6),4), #c(log(8*10e6),log(8*10e6),10,10), ## sum wc = 12
    omega_0ij = omega_0ij,
    logq_f = rep(log(0.5), 5),
    b = rep(0,nyear),  
    logpi_acomp = rep(log(50),nfleets_acomp),
    logSDR = 1.4,
    ## structure is fleet x alpha, beta x time block (1 for now)x sex 
    log_fsh_slx_pars = log_fsh_slx_pars,
    log_srv_slx_pars = log_srv_slx_pars,
    epsilon_tau = rep(log(5),nspace),
    mort_k = mort_k ## mortality
  )
  
  ## initial matrix ----
  # load(here('input','input_data',"Neqn.rdata"))
  
  # alpha_g1 <- c(62.8329, 63.6959, 33.8898, 54.1045, 64.2127) ## trap ll twl std strs
  # beta_g1 <- c(7.04483, 3.09715, 1.41494, 4.55724, 12.9197)
  
  # parms <- list( # Just start all the simluations with the same initial conditions 
  #      logRinit = parms.scalar$logRinit+log(rmul),
  #      logh = parms.scalar$logh,
  #      logMinit = parms.scalar$logMinit,
  #      logSDsurv = parms.scalar$logSDsurv,
  #      #logSDR = log(1.4),
  #      logphi_catch = parms.scalar$logphi_catch,
  #      #logphi_survey = log(11.33),
  #      # logSDF = log(0.1),
  #      # Selectivity parameters 
  #      psel_fish = parms.sel$value[parms.sel$source == 'fish'],
  #      psel_surv = parms.sel$value[parms.sel$source == 'survey'],
  #      initN = initN,
  #      Rin = Rdev,
  #      PSEL = PSEL
  #    )
  load(here('input','input_data','unfished_ALK.rdata'))
  load(here('input','input_data','mla_yais.rdata')) ## from prelim runs, for ssb0
  
  ## Return df ----
  df <-list(    
    #* MODEL STRUCTURE ----
    nage = nage,
    age = age,
    nyear = nyear,
    tEnd = tEnd, # The extra year is to initialize 
    years = years,
    nstocks = nstocks,
    nspace = nspace,
    LBins = LBins,
    move = move,
    nmgmt_reg = ncol(phi_fm),
    
    #* FLEETS STRUCTURE ----
    nfleets_surv = nfleets_surv,
    nfleets_fish = nfleets_fish,
    nfleets_acomp = nfleets_acomp,
    nfleets_lcomp = nfleets_lcomp,
    selShape_fish = selShape_fish, ## 0 and 1 logistic, 2 dome normal, 3 dome gamma
    selShape_surv = selShape_surv,
    selType_fish = selType_fish, ## 0 for age, 1 for length-based
    selType_surv = selType_surv,
    
    fltnames_surv = fltnames_surv,
    fltnames_fish = fltnames_fish,
    fltnames_acomp = fltnames_acomp,
    fltnames_lcomp = fltnames_lcomp,
    
    phi_if_surv = phi_if_surv,
    phi_if_fish = phi_if_fish,
    phi_if_acomp = phi_if_acomp,
    phi_ff_acomp=phi_ff_acomp,
    phi_ki = phi_ki,
    phi_im = phi_im,
    phi_ik2 = t(phi_ik2),
    phi_ij = phi_ij,
    phi_fm = phi_fm,
    tau_ki = tau_ki,
    phi_fm_acomp = phi_fm_acomp,
    phi_fm_acomp2 = t(phi_fm_acomp2),
    phi_lcomp_fm =phi_lcomp_fm,
    acomp_flt_type = acomp_flt_type,
    
    #* DEMOG ----
    X_ijas = X_ijas,
    omega_ais = omega_ais,
    Linf_yk = growthPars$Linf_yk,
    kappa_yk = growthPars$kappa_yk,
    sigmaG_yk = growthPars$sigmaG_yk,
    L1_yk = growthPars$L1_yk,
    wtatlen_kab = wtatlen_kab,
    mat_ak = mat_ak, ## maturity
    # mort_k = mort_k, ## mortality
    unfished_ALK_F = unfished_ALK_F,
    mla_yais = mla_yais-1,
    
    #* DATA ----
    acomp_yafs_obs = OM_agecomps_yafs,
    surv_yf_obs = as.matrix(survey), # Make sure the survey has the same length as the catch time series
    surv_yf_err = as.matrix(survey_err), 
    age_error = as.matrix(ageerr_ExpAge[,2:ncol(ageerr_ExpAge)]),
    age_error_SD = as.matrix(ageerr_SD[,2:ncol(ageerr_SD)]),
    catch_yf_obs = as.matrix(catch),
    catch_yf_error = as.matrix(catch_yf_error),
    discard = omdis,
    
    #* SELEX ----
    fish_selex_yafs = OM_fish_selex_yafs,
    surv_selex_yafs = OM_surv_selex_yafs,
    fsh_blks = fsh_blks, 
    srv_blks = srv_blks,
    fsh_blks_size = fsh_blks_size,
    srv_blks_size = srv_blks_size,
    # Neqn = Neqn, ## solve(I-Mat2) ## load externally once M is setup
    # fsh_blks = t(fsh_blks), ## currently not ready to be fleet-specific
    # srv_blks = t(srv_blks),
    #* ADDL PARS ----
    parms = parms,
    bfuture = bfuture
  )
  return(df) ## this has all the data in a format ready for estimation
}
