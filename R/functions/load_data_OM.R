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
  omega_0ij = matrix(0, nrow = nspace, ncol = nspace)
  diag(omega_0ij) <- 1
  
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
  
  selShape_fish <- c(rep(0,4),2,2,3,2,2) ## 0 and 1 logistic, 2 dome normal, 3 dome gamma
  selShape_surv <- c(rep(0,10)) ## 0 and 1 logistic, 2 dome normal, 3 dome gamma
  
  selType_fish <- as.numeric(fltnames$SELTYPE[fltnames$COMM])-1
  ## note that the first two acomp fleets are already inside seltype fish
  selType_surv <- as.numeric(c(fltnames$SELTYPE[fltnames$SURV],fltnames$SELTYPE[fltnames$ACOMP][3:8]))-1
  
  
  nfleets_fish <- length(fltnames$NAME[fltnames$COMM])
  nfleets_surv <- length(fltnames$NAME[fltnames$SURV])
  nfleets_acomp <- length(fltnames$NAME[fltnames$ACOMP])
  nfleets_lcomp <- length(fltnames$NAME[fltnames$LCOMP])

  # Catch ----
  catch <- read.csv(here("input","input_data","OM_catch.csv"))
  catch[is.na(catch)] <- -1
  ## Discard ----
  load(here("input","input_data","OM_discard.csv")) ## loads as omdis
  
  ## Selex ----
  load(here('input','input_data',"OM_fish_selex_yafs.rdata"))
  load(here('input','input_data',"OM_surv_selex_yafs.rdata"))
  
  ## time blocks by fleet
  srv_blks <- matrix(0, nrow = tEnd, ncol = nfleets_surv+(nfleets_acomp-2))
  fsh_blks <- matrix(0, nrow = tEnd, ncol = nfleets_fish)
  
  selShape_fish <- c(rep(0,4),2,2,3,2,2) ## 0 and 1 logistic, 2 dome normal, 3 dome gamma
  selShape_surv <- c(rep(0,nfleets_surv+(nfleets_acomp-2))) ## 0 and 1 logistic, 2 dome normal, 3 dome gamma
  
  selType_fish <- as.numeric(fltnames$SELTYPE[fltnames$COMM])-1
  ## note that the first two acomp fleets are already inside seltype fish
  selType_surv <- as.numeric(c(fltnames$SELTYPE[fltnames$SURV],fltnames$SELTYPE[fltnames$ACOMP][3:8]))-1
  
  if(length(selType_surv) != length(selShape_surv)) stop("seltype surv length doesn't match selshape surv")
  # Survey ----
  survey <- read.csv(here("input","input_data",'OM_indices.csv'))
  survey[is.na(survey)] <- -1 ## flag for numeric TMB checks
  survey_err <- read.csv(here("input","input_data",'OM_indices_sigma.csv'))


  ## Comps ----
  ## Len comps [these are arrays by fleet]
  load(here("input","input_data",'OM_lencomps_female.rdata'))
  load(here("input","input_data",'OM_lencomps_male.rdata'))
  
  ## Age comps. Note that AK is not sex-specific, and therefore duplicatd
  load(here("input","input_data",'OM_agecomps_female.rdata'))
  load(here("input","input_data",'OM_agecomps_male.rdata'))

  ## Aging Error ----
  ## M X age
  load(here("input","input_data",'ageerr_ExpAge.rdata'))
  load(here("input","input_data",'ageerr_SD.rdata'))

  ## Phi objects ----
  ## setup phi (spatial matching matrix) depending on spatial setup
  ## MAKE A NSPACE == 3 AND 1 OPTION FOR COMBINING
  spmat <- data.frame(subarea = c('A1',"A2","B2","B1","C2","C1"),
                      stock = c("R4","R3","R3","R2","R2","R1"),
                      mgmt = c("AI","AK", rep("BC",2), rep("CC",2)))
  if(nspace == 6){ ## OM
    
    ## phi_survy
    phi_if_surv <- matrix(0, nrow = nfleets_surv, ncol = nspace)
    rownames(phi_if_surv) <- names(survey)
    colnames(phi_if_surv) <- spmat$subarea

    phi_if_surv[1,1] <-  phi_if_surv[2,2] <-  
      phi_if_surv[3:4,3:4]<-  phi_if_surv[5,5:6] <- 1
    
    phi_if_acomp <- matrix(0, nrow = nfleets_acomp, ncol = nspace)
    rownames(phi_if_acomp) <- fltnames_acomp
    colnames(phi_if_acomp) <- spmat$subarea
    phi_if_acomp[1,1] <-  phi_if_acomp[2:3,2] <-  
      phi_if_acomp[4:6,3:4]<-  phi_if_acomp[7:8,5:6] <- 1
    ## phi_fish
    phi_if_fish <- matrix(0, nrow = nfleets_fish, ncol = nspace) ## placeholder for fishing fleets
    rownames(phi_if_fish) <- names(catch)[2:ncol(catch)]
    colnames(phi_if_fish) <-  spmat$subarea
    
    phi_if_fish[c(1,3),1] <- phi_if_fish[c(2,4),2] <-   phi_if_fish[5:7,3:4] <-  
      phi_if_fish[c(8,9),5:6] <-  1
    
    ## phi_ik
    phi_ki <-  matrix(0, ncol = nspace, nrow = nstocks) ## nesting of subareas within stocks, for recruitment purposes
    rownames(phi_ki) <- unique(spmat$stock)
    colnames(phi_ki) <- spmat$subarea
    
    phi_ki[1,1] <-  phi_ki[2,2:3] <-  phi_ki[3,4:5]<-  phi_ki[4,6]  <- 1
    phi_ik2 <- matrix(apply(phi_ki,2, function(x)which(x == 1))-1) ## a vector for par subsetting, the columns are subareas
    
    ## phi_ij [eq 6]
    phi_ij <-  matrix(1, ncol = nspace, nrow = nspace) ## 0 indicates  subareas comprise THE SAME stock
    rownames(phi_ij) = colnames(phi_ij) = spmat$subarea
    diag(phi_ij) <- phi_ij[4,5] <- phi_ij[5,4] <- 0
    
      
    ## phi_fm
    phi_fm <- matrix(0, nrow = nfleets_fish, ncol = 3)
    rownames(phi_fm) = names(catch)[2:ncol(catch)]
    colnames(phi_fm) = c('AK','BC','WC')
    phi_fm[1:4,1] <- phi_fm[5:7,2]  <- phi_fm[8:9,3]  <- 1
    
    ## same as above but for comps (mix of fisheries & surveys)
    phi_fm_acomp <- matrix(0, nrow = nfleets_acomp, ncol = 3)
    rownames(phi_fm_acomp) = fltnames_acomp
    colnames(phi_fm_acomp) = c('AK','BC','WC')
    phi_fm_acomp[1:3,1] <- phi_fm_acomp[4:6,2]  <- phi_fm_acomp[7:8,3]  <- 1
    phi_fm_acomp2 <- matrix(apply(phi_fm_acomp,1, function(x)which(x == 1))-1) ## a vector for par subsetting, the columns are survey fleets
    
    acomp_flt_type <- matrix(0, ncol = nfleets_acomp) ## 0 is commercial, 1 is survey
    acomp_flt_type[3:8] <- 1
    colnames(acomp_flt_type) <- fltnames_acomp

    
    phi_lcomp_fm <- matrix(0, nrow = nfleets_lcomp, ncol = 3)
    rownames(phi_lcomp_fm) = fltnames_lcomp
    colnames(phi_lcomp_fm) = c('AK','BC','WC')
    phi_lcomp_fm[1:6,1] <- phi_lcomp_fm[7:9,2]  <- phi_lcomp_fm[10,3]  <- 1
    

    ## tau_ki
    tau_ki <-  matrix(0, ncol = nspace, nrow = nstocks) ## nesting of subareas within stocks, for recruitment purposes
    rownames(tau_ki) <- unique(spmat$stock)
    colnames(tau_ki) <- spmat$subarea
    tau_ki[1,1] <-   tau_ki[4,6]  <- 1 ## 100% of recruitment in stock
    tau_ki[2,2:3] <-  tau_ki[3,5:4] <-  c(0.75,0.25) ## A2 and C1 are larger
  } else {
    phi_if_surv <- matrix(rbinom(nfleets_surv*nspace,1,0.5), byrow = TRUE, nrow = nfleets_surv, ncol = nspace) ## placeholder for alternative spatial stratifications
    phi_if_fish <- matrix(c(0,1,1,1,1,0), nrow = nfleets_fish, ncol = nspace)  ## placeholder for fishing fleets
    phi_ki <- matrix(c(1,0,0,1), byrow = TRUE, nrow = nstocks, ncol = nspace) ## placeholder for alternative spatial stratifications
    phi_ik2 <- apply(phi_ki,2, function(x)which(x == 1))-1 ## a vector for par subsetting, the columns are subareas
   
    # phi_fm <- matrix(0, nrow = nspace, ncol = 2)
    # phi_fm[1,1] <- phi_fm[2,2] <- 1
    phi_fm <- matrix(0, nrow = nfleets_fish, ncol = 3)
    phi_fm[1:2,1] <- phi_fm[3,3]  <- 1
     ## autogenerate stock-distinction matrix
    phi_ij <-  matrix(NA, byrow = TRUE, ncol = nspace, nrow = nspace)
    for(i in 1:nspace){
      for(j in 1:nspace){
        phi_ij[i,j] = ifelse(phi_ik2[i] == phi_ik2[j],0,1)
      }
    }
    tau_ki <- matrix(c(0.25,0.75,0.9,0.1), nrow = nstocks, byrow = TRUE, ncol = nspace) ## placeholder for alternative spatial stratifications
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


  ## Parms List ----
  ## things that will get estimated later on, everthing else is FIXED
  parms <- list(
    logh_k = c(0.7,0.7,0.88,0.7),
    logR_0k = c(log(8*10e6),log(8*10e6),10,4), ## sum wc = 12
    omega_0ij = omega_0ij,
    logq_f = rep(log(0.5), 5),
    b = rep(0,nyear),
    logSDR = 1.4,
    ## structure is fleet x alpha, beta x time block (1 for now)x sex 
    log_fsh_slx_pars = array(0.2, dim = c(nfleets_fish,2,1,2)),
    log_srv_slx_pars =  array(0.4, dim = c( nfleets_surv+(nfleets_acomp-2),2,1,2))
  )


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
    phi_ki = phi_ki,
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
    mat_age = rep(0.2,nage), ## mortality
    unfished_ALK_F = unfished_ALK_F,
    mla_yais=mla_yais,
    
    #* DATA ----
    surv_yf_obs = as.matrix(round(survey)), # Make sure the survey has the same length as the catch time series
    surv_yf_err = as.matrix(survey_err), 
    age_error = as.matrix(ageerr_ExpAge[,2:ncol(ageerr_ExpAge)]),
    age_error_SD = as.matrix(ageerr_SD[,2:ncol(ageerr_SD)]),
    catch_yf_obs = as.matrix(catch),
    discard = omdis,
    
    #* SELEX ----
    fish_selex_yafs = OM_fish_selex_yafs,
    surv_selex_yafs = OM_surv_selex_yafs,
    fsh_blks = fsh_blks, ## currently not ready to be fleet-specific
    srv_blks = srv_blks,
    #* ADDL PARS ----
    parms = parms,
    bfuture = bfuture
  )
  return(df) ## this has all the data in a format ready for estimation
}
