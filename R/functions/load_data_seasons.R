## Load the hake data
# year and age input 
load_data_seasons <- function(nspace = 6, 
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
  
  years <- 1960:(myear+yr_future)
  nyear <- length(years)
  tEnd <- length(years)
  age <- 0:70 
  nage <- length(age)
 
  # Maturity ----
  load(here("input","input_data","OM_maturity_ak.rdata")) ## ak is age, stock

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
  
  selType_fish <- fltnames$SELTYPE[fltnames$COMM]
  selType_surv <- fltnames$SELTYPE[fltnames$SURV]
  
  
  nfleets_fish <- length(fltnames$NAME[fltnames$COMM])
  nfleets_surv <- length(fltnames$NAME[fltnames$SURV])
  nfleets_acomp <- length(fltnames$NAME[fltnames$ACOMP])
  nfleets_lcomp <- length(fltnames$NAME[fltnames$LCOMP])

  # Catch ----
  catch <- read.csv(here("input","input_data","OM_catch.csv"))
  
  ## Discard ----
  load(here("input","input_data","OM_discard.csv")) ## loads as omdis
  
  ## Selex: still need SURVEY ----
  load(here('input','input_data',"OM_fish_selex_yafs.rdata"))
  load(here('input','input_data',"OM_survey_selex_yafs.rdata"))
  
  # Survey ----
  # survey <- read.csv(here("input","data",'acoustic survey.csv'))
  survey <- read.csv(here("input","input_data",'OM_indices.csv'))
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
    
    ## phi_fish
    phi_if_fish <- matrix(0, nrow = nfleets_fish, ncol = nspace) ## placeholder for fishing fleets
    rownames(phi_if_fish) <- names(catch)[2:ncol(catch)]
    colnames(phi_if_fish) <-  spmat$subarea
    
    phi_if_fish[c(1,3),1] <- phi_if_fish[c(2,4),2] <-   phi_if_fish[5:7,3:4] <-  
      phi_if_fish[c(8,9),5:6] <-  1
    
    ## phi_ik
    phi_ik <-  matrix(0, ncol = nspace, nrow = nstocks) ## nesting of subareas within stocks, for recruitment purposes
    rownames(phi_ik) <- unique(spmat$stock)
    colnames(phi_ik) <- spmat$subarea
    
    phi_ik[1,1] <-  phi_ik[2,2:3] <-  phi_ik[3,4:5]<-  phi_ik[4,6]  <- 1
    phi_ik2 <- apply(phi_ik,2, function(x)which(x == 1))-1 ## a vector for par subsetting, the columns are subareas
    
    ## phi_ij [eq 6]
    phi_ij <-  matrix(1, ncol = nspace, nrow = nspace) ## 0 indicates  subareas comprise THE SAME stock
    rownames(phi_ij) = colnames(phi_ij) = spmat$subarea
    diag(phi_ij) <- phi_ij[4,5] <- phi_ij[5,4] <- 0
    
      
    ## phi_fm
    phi_fm <- matrix(0, nrow = nfleets_fish, ncol = 3)
    rownames(phi_fm) = names(catch)[2:ncol(catch)]
    colnames(phi_fm) = c('AK','BC','WC')
    phi_fm[1:2,1] <- phi_fm[3:5,2]  <- phi_fm[6:9,3]  <- 1
    
    ## tau_ki
    tau_ki <-  matrix(0, ncol = nspace, nrow = nstocks) ## nesting of subareas within stocks, for recruitment purposes
    rownames(tau_ki) <- unique(spmat$stock)
    colnames(tau_ki) <- spmat$subarea
    tau_ki[1,1] <-   tau_ki[4,6]  <- 1 ## 100% of recruitment in stock
    tau_ki[2,2:3] <-  tau_ki[3,4:5] <-  0.5 ## split 50/50 for now
  } else {
    phi_if_surv <- matrix(rbinom(nfleets_surv*nspace,1,0.5), byrow = TRUE, nrow = nfleets_surv, ncol = nspace) ## placeholder for alternative spatial stratifications
    phi_if_fish <- matrix(c(0,1,1,1,1,0), nrow = nfleets_fish, ncol = nspace)  ## placeholder for fishing fleets
    phi_ik <- matrix(c(1,0,0,1), byrow = TRUE, nrow = nstocks, ncol = nspace) ## placeholder for alternative spatial stratifications
    phi_ik2 <- apply(phi_ik,2, function(x)which(x == 1))-1 ## a vector for par subsetting, the columns are subareas
   
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
    logh_k = rep(log(0.2),4),
    logRinit = log(1e5)
  )

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
    
    #* FLEETS STRUCTURE ----
    nfleets_surv = nfleets_surv,
    nfleets_fish = nfleets_fish,
    nfleets_acomp = nfleets_acomp,
    nfleets_lcomp = nfleets_lcomp,
    selType_fish = selType_fish,
    selType_surv = selType_surv,
    
    fltnames_surv = fltnames_surv,
    fltnames_fish = fltnames_fish,
    fltnames_acomp = fltnames_acomp,
    fltnames_lcomp = fltnames_lcomp,

    phi_if_surv = phi_if_surv,
    phi_if_fish = phi_if_fish,
    phi_ik = phi_ik,
    phi_ik2 = phi_ik2,
    phi_ij = phi_ij,
    phi_fm = phi_fm,
    tau_ki = tau_ki,
    
    #* DEMOG ----
    X_ijas = X_ijas,
    omega_ais = omega_ais,
    Linf_yk = growthPars$Linf_yk,
    kappa_yk = growthPars$kappa_yk,
    sigmaG_yk = growthPars$sigmaG_yk,
    L1_yk = growthPars$L1_yk,
    wtatlen_kab = wtatlen_kab,
    mat_ak = mat_ak,
    #* DATA ----
    survey = survey, # Make sure the survey has the same length as the catch time series
    survey_err = survey_err, #ac.data$ss.error, # Make sure the survey has the same length as the catch time series
    survey = survey, #ac.data$ss.survey,
    age_error = ageerr_ExpAge,
    age_error_sd = ageerr_SD,
    catch = catch,
    discard = omdis,
    fish_selex_yafs = OM_fish_selex_yafs,
    surv_selex_yafs = OM_surv_selex_yafs,
    
    #* ADDL PARS ----
    parms = parms,
    b = b,
    bfuture = bfuture,
    logSDR = logSDR
  )
  return(df) ## this has all the data in a format ready for estimation
}
