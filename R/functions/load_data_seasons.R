## Load the hake data
# year and age input 
load_data_seasons <- function(nspace = 6, 
                              nstocks = 4,
                              myear = 2018,
                              move = TRUE, 
                              LBins = 81,

                              # movemaxinit = 0.35, 
                              # movefiftyinit = 6,
                              # nsurvey = 2, 
                              logSDR = 1.4, 
                              bfuture = 0.5,
                              # moveout = 0.85, 
                              # movesouth = 0.05,
                              # moveinit = NA, 
                              # moveslope = 0.9,
                              # selectivity_change = 0,
                              yr_future  = 0
                              # sel_hist = 1
                              ){
  
  #' @nspace = Spatial areas
  #' @nstocks = demographic stocks
  #' @myear = Final year that contains empirical data
  #' @move = do you want to enable movement among spaces?
  #' @logSDR Standard deviation of recruitment
  #' @bfuture recruitment bias adjustment in the future - scalar
  #' @yr_future Create dummy data for future years
  
  
  # if(is.na(moveinit)){
  #   if(nspace == 2){
  #   moveinit <-  c(0.25,0.75)
  #   } else if(nspace > 2){
  #     moveinit <- c(runif(nspace, 0.01, 0.45)) ## placeholder movement rates (random)
  #     
  #   }
  # }

  years <- 1960:(myear+yr_future)
  nyear <- length(years)
  tEnd <- length(years)*nseason
  age <- 0:70 # 0:95
  
  ## Age stuff
  nage <- length(age)
  # msel <- rep(1,nage)

  ## for later use
  # recruitmat <- matrix(0, nspace) # 4 is seasonality 
  # recruitmat[1] <- 1 # 10 percent change of spawning north
  # recruitmat[2] <- 1 # 90 percent of spawning south
  
  
  # Maturity ----
  load(here("input","input_data","OM_maturity_ak.rdata")) ## ak is age, stock
  
  
  ## Movement ----
  # movefifty <- movefiftyinit
  # 
  # movemax <- rep(movemaxinit,nseason)
  # 
  # movemat <- array(0, dim = c(nspace, nage, nseason, nyear)) # Chances of moving in to the other grid cell 
  # 
  # 
  # if(nspace == 1){
  #   move = FALSE
  # }else{
  #   move = TRUE
  # }
  # 
  # if(move == TRUE){
  #   for(j in 1:nspace){
  #     for(i in 1:nseason){
  #       movemat[j,,i,] <- movemax[i]/(1+exp(-moveslope*(age-movefifty))) ## saturation for movement, may fix instead
  #       
  #     }
  #   }
  #   
  #   movemat[,1:2,,] <- 0 # Recruits and 1 year olds don't move
  #   
  #   if(nseason == 4){ # For the standard model
  #     
  #     movemat[1,3:nage,2:3,] <- movesouth # Don't move south during the year
  #     movemat[1,3:nage,1,] <- moveout*0.5 # continuing south movement at spawning time
  #     
  #     movemat[1,3:nage,nseason,] <- moveout
  #     movemat[2,3:nage,nseason,] <- movesouth
  #   }
  #   # movemat[1,11:nage,nseason] <- 0
  #   # movemat[2,11:nage,nseason] <- 0
  #   
  #   
  #   
  #   # move.init <- array(0.5, dim = c(nspace, nage))
  #   # 
  #   # move.init[1,] <- 0.3
  #   # move.init[2,] <- 0.7
  #   move.init <- moveinit
    
  # }else{
  #   move.init <- 1
  # }
  
  ## placeholder for X_ija -this will need to get converted from length
  if(move == FALSE) X_ija <- array(rep(0, nspace*nspace*nage), c(nspace,nspace,nage))
  if(move == TRUE) {
    X_ija <- array(runif( nspace*nspace*nage, 0.01,0.05),  c(nspace,nspace,nage))
    
    ## for placeholder; if the rows are summing to greater than one set to zero
    for(n in 1:nage){
      for(i in 1:nspace){
        for(j in 1:nspace){
          if(sum(X_ija[i,1:j,n]) > 1)  X_ija[i,j:nspace,n] <- 0
        } ## end j space
      } ## end i space
    } ## end n ages
  } ## end move == TRUE

  omega_ai <- matrix(rep(0, nage*nspace), c(nage,nspace))## eigenvector for stable spatial distribution at age
  for(a in 1:nage){
    omega_ai[a,] <- eigen(X_ija[,,a])$values 
    omega_ai[a,][which(omega_ai[a,] < 0)] <- 0.05
  }

  
  # Weight at length ----
  load(here("input","input_data","OM_wtatlen_kab.rdata")) ## a and be are pars of al^b
  
  # Growth ----
  load(here("input","input_data","OM_growthPars.rdata")) 
  
  ## Mortality ----
  load(here('input','input_data','M_k.rdata'))
  # wage_ss <- read.csv(here("input","data",'wage_ss.csv'))
  # wage_ss <- wage_ss[wage_ss$Yr %in% years,]
  # wage_unfished <- read.csv(here("input","data",'unfished_waa.csv'))
  # wage_ssb <- wage_ss[wage_ss$Fleet == -2,paste('X',age, sep = '')]
  # wage_ssb[[1]] <- unname(wage_ssb[[1]])
  # wage_catch <- wage_ss[wage_ss$Fleet == 1 ,paste('X',age, sep = '')]
  # wage_survey <- wage_ss[wage_ss$Fleet == 2,paste('X',age, sep = '')]
  # wage_mid <- wage_ss[wage_ss$Fleet == -1,paste('X',age, sep = '')]
  
  # if(yr_future>0){ // Namiong issues here, fix in later update
  #   tmp_ssb <- matrix(rep(wage_ssb[1,], each = yr_future), nrow = yr_future)
  #   wage_ssb <- cbind(wage_ssb,tmp_ssb)
  #    
  #   
  #   wage_catch <- wage_ss[wage_ss$Fleet == 1 ,paste('X',age, sep = '')]
  #   wage_survey <- wage_ss[wage_ss$Fleet == 2,paste('X',age, sep = '')]
  #   wage_mid <- wage_ss[wage_ss$Fleet == -1,paste('X',age, sep = '')]
  #   
  #   
  #   
  #   
  # }
  # 
  
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
  
  ## nfleets should be input
  ## auto-read nfleets_surv etc from this
  # flag_fleets <- array(NA, dim = c(4,8,nyear) ) 
  ## row 1 = catch, row 2 = survey biomass, row 3 = acomp, row 4 = lcomp
  
  # ## populate catch fleets [needs to be yearly; perhaps automate upon readin of CSV
  # flag_fleets[,1:4,] <- c(1,0,0,0)
  # ## populate survey fleets
  # flag_fleets[,5:7,] <- c(0,1,0,0)
  # ## who has acomps?
  # flag_fleets[,8,] <- c(1,0,0,0)
  # ## who has lcomps?
  # flag_fleets[,1,] <- c(1,0,0,0)
  
  # nfleets_fish <- sum(flag_fleets[1,] == 1) -1 ## zero-indexed
  # nfleets_surv <- sum(flag_fleets[2,] == 1) -1 ## zero-indexed
  # nfleets_acomp <- sum(flag_fleets[3,] == 1) -1 ## zero-indexed
  # nfleets_lcomp <- sum(flag_fleets[4,] == 1) -1 ## zero-indexed
  
  
  ## make phi objects (which column corresponds to which fleet)
  # phi_fleet_catch <- matrix(NA, ncol = ncol(flag_fleets));   idx = 0 
  # for(i in 1:ncol(flag_fleets)){
  #   if(flag_fleets[1,i,] == 1){
  #     phi_fleet_catch[i] <- idx; idx = idx+1 ## up counter
  #   }
  # }
  #  phi_fleet_surv <- matrix(NA, ncol = ncol(flag_fleets));   idx = 0 
  # for(i in 1:ncol(flag_fleets)){
  #   if(flag_fleets[2,i,] == 1){
  #     phi_fleet_surv[i] <- idx; idx = idx+1 ## up counter
  #   }
  # }
  # 
  # phi_fleet_acomp <- matrix(NA, ncol = ncol(flag_fleets));   idx = 0 
  # for(i in 3:ncol(flag_fleets)){
  #   if(flag_fleets[3,i,] == 1){
  #     phi_fleet_acomp[i] <- idx; idx = idx+1 ## up counter
  #   }
  # }
  # phi_fleet_lcomp <- matrix(NA, ncol = ncol(flag_fleets));   idx = 0 
  # for(i in 4:ncol(flag_fleets)){
  #   if(flag_fleets[4,i,] == 1){
  #     phi_fleet_lcomp[i] <- idx; idx = idx+1 ## up counter
  #   }
  # }
  
 

  
  # Catch ----
  # catch <- read.csv(here("input","data",'hake_totcatch.csv'))
  # catch2 <- cbind(catch$year, catch$Fishery, catch$Fishery, catch$Fishery) ## multifleet placeholder
  catch <- read.csv(here("input","input_data","OM_catch.csv"))

  # Survey ----
  # survey <- read.csv(here("input","data",'acoustic survey.csv'))
  survey <- read.csv(here("input","input_data",'OM_indices.csv'))
  survey_err <- read.csv(here("input","input_data",'OM_indices_sigma.csv'))
  
  # survey <- read.csv(here("input","data",'survey.csv'))
  # survey2 <- as.matrix(read.csv("input/cleaned/clean_survey.csv"))## this needs to be built into load_data_seasons
  # survey_x2  <- rep(2, length(years)) ## we have obs from 1970+
  # survey_x2[1:10] <- -2 ## a -2 if no survey, 2 if survey occured
  # nfleets_surv <- ncol(survey) 

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

  
  # if(nseason == 4){
  # surveyseason <- 3
  # 
  # }else{
  #   surveyseason <- floor(nseason/2)
  # }
  # 
  # 
  # if(nseason == 1){
  #   surveyseason = 1
  # }

  ## Phi objects ----
  ## setup phi (spatial matching matrix) depending on spatial setup
  spmat <- data.frame(subarea = c('A1',"A2","B2","B1","C2","C1"),
                      stock = c("R4","R3","R3","R2","R2","R1"),
                      mgmt = c("AI","AK", rep("BC",2), rep("CC",2)))
  if(nspace == 6){ ## OM
    
    ## phi_survy
    phi_if_surv <- matrix(0, nrow = nfleets_surv, ncol = nspace)
    
    rownames(phi_if_surv) <- names(survey)
    colnames(phi_if_surv) <- rev(spmat$subarea)
    
    phi_if_surv[1,1:2] <-  phi_if_surv[2,3:4] <-  
      phi_if_surv[3,3:4]<-  phi_if_surv[4,5] <-  phi_if_surv[5,6] <- 1
    
    ## phi_fish
    phi_if_fish <- matrix(0, nrow = nfleets_fish, ncol = nspace) ## placeholder for fishing fleets
    rownames(phi_if_fish) <- names(catch)[2:ncol(catch)]
    colnames(phi_if_fish) <- rev(spmat$subarea)
    
    phi_if_fish[1:2,1:2] <-  phi_if_fish[3:5,3:4] <-  
      phi_if_fish[c(6,8),5] <-        phi_if_fish[c(7,9),6] <- 1
    
  
    ## phi_ik
    phi_ik <-  matrix(0, ncol = nspace, nrow = nstocks) ## nesting of subareas within stocks, for recruitment purposes
    rownames(phi_ik) <- unique(rev(spmat$stock))
    colnames(phi_ik) <- rev(spmat$subarea)
    
    phi_ik[1,1] <-  phi_ik[2,2:3] <-  phi_ik[3,4:5]<-  phi_ik[4,6]  <- 1
    phi_ik2 <- apply(phi_ik,2, function(x)which(x == 1))-1 ## a vector for par subsetting, the columns are subareas
    
    ## phi_ij
    phi_ij <-  matrix(0, ncol = nspace, nrow = nspace) ## 1 indicates if subareas comprise DISTINCT stocks
    rownames(phi_ij) = colnames(phi_ij) = rev(spmat$subarea)
    phi_ij[1,2:nspace] <- phi_ij[2,c(1,3:nspace)] <- phi_ij[3,c(1:2,4:nspace)] <- phi_ij[4,c(1:3,nspace)]<- phi_ij[5,c(1:3,nspace)] <- phi_ij[6,c(1:4)] <- 1
      
    ## phi_fm
    phi_fm <- matrix(0, nrow = nfleets_fish, ncol = 3)
    rownames(phi_fm) = names(catch)[2:ncol(catch)]
    phi_fm[1:2,1] <- phi_fm[3:5,2]  <- phi_fm[6:9,3]  <- 1
    
    ## tau_ik
    tau_ik <-  matrix(0, ncol = nspace, nrow = nstocks) ## nesting of subareas within stocks, for recruitment purposes
    rownames(tau_ik) <- unique(rev(spmat$stock))
    colnames(tau_ik) <- rev(spmat$subarea)
    tau_ik[1,1] <-   tau_ik[4,6]  <- 1 ## 100% of recruitment in stock
    tau_ik[2,2:3] <-  tau_ik[3,4:5] <-  0.5 ## split 50/50 for now
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
    tau_ik <- matrix(c(0.25,0.75,0.9,0.1), nrow = nstocks, byrow = TRUE, ncol = nspace) ## placeholder for alternative spatial stratifications
  }

  
  # age_survey <- as.data.frame(matrix(-1, nyear,dim(age_survey.df)[2]))
  # names(age_survey) <- names(age_survey.df)
  # age_survey$year <- years
  # age_catch <- as.data.frame(matrix(-1, nyear,dim(age_catch.df)[2]))
  # names(age_catch) <- names(age_catch.df)
  # age_catch$year <- years
  # 
  # for (i in 1:dim(age_survey.df)[1]){
  #   idx <- which(age_survey$year == age_survey.df$year[i])
  #   age_survey[idx,] <-age_survey.df[i,]
  #   
  # }
  # 
  # for (i in 1:dim(age_catch.df)[1]){
  #   idx <- which(age_catch$year == age_catch.df$year[i])
  #   age_catch[idx,] <-age_catch.df[i,]
  #   
  # }
  
  # Load parameters from the assessment 
  # initN <- rev(read.csv(here("input","data",'Ninit_MLE.csv')))[,1]
  # Rdev <- read.csv(here("input","data",'Rdev_MLE.csv'))[,1]
  # PSEL <- as.matrix(read.csv(here("input","data",'p_MLE.csv'))) ## time varying selex pars for fihsery(?)
  #Fin <- assessment$F0
  
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
  # 
  #b <- matrix(1, tEnd)
  # b <- as.matrix(read.csv(here("input","data",'b_input.csv'))) ## this is rec penalty
 
  
  # if(move == TRUE){
  #    mul <- 1.015
  # }else{
  #  mul <- 1
  #  }

  # load parameters specifically for hake 
  # parms.scalar <- read.csv(here("input","data","parms_scalar.csv"))
  # parms.sel <- read.csv(here("input","data",'selectivity.csv'))
  # initN <- as.matrix(read.table(here("input","data",'initN.csv'))) ## initial n at age?
  
  # Rdev <- as.matrix(read.csv(here("input","data",'Rdev.csv')))
  
  
  # if(sel_hist == 1){
  # PSEL <- as.matrix(read.csv(here("input","data",'PSEL.csv')))
  # }else{
  # PSEL <- matrix(0, 5, 28)
  # }
  
  
  # if(nseason == 4 & nspace == 2){
  # Fnseason <- matrix(NA, 2,4)
  # 
  # #Fnseason[1,] <- c(0.0,0.4,0.50,0.1) # Must add to one 
  # Fnseason[1,] <- c(0.001,0.188,0.603,0.208)
  # #Fnseason[2,] <- c(0.0,0.4,0.50,0.1) # Must add to onec
  # Fnseason[2,] <- c(0.000,0.317,0.382,0.302)/sum(c(0.000,0.317,0.382,0.302)) # Divide by sum to sum to 1 
  #   
  # }else{
  #   Fnseason <- matrix(NA, nspace, nseason)
  #   Fnseason[1:nspace,] <- 1/nseason # Equally distributed catch
  #   
  #   
  # }


  # rmul <-1
  # 
  # if(nspace == 2){
  #   rmul <- 1.1
  # }
  
  ## Parms List ----
  ## things that will get estimated later on, everthing else is FIXED
  parms <- list(
    logh_k = rep(0.2, 4),
    logRinit = 1e7
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
  #    
  #    psel<- matrix(NA,nspace, 5) 
  #    
  #    for(i in 1:nspace){
  #      #psel[i,] <- c(2.8476, 0.973,0.3861,0.1775,0.5048) # USA selectivity 
  #      psel[i,] <- parms$psel_fish
  #      
  #    }
  #   
  #    if(nspace == 2){
  #    psel[1,] <- c(1,1,1,1,1)
  #    
  #    }
       
       
# Flag if there's a selectivity change in that year     
     # selYear <- 1991
     # 
     # flag_sel <- rep(0,nyear)
     # flag_sel[which(years == selYear):which(years == myear)] <- 1
     
     
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
    
    #* FLEETS STRUCTURE ----
    nfleets_surv = nfleets_surv,
    nfleets_fish = nfleets_fish,
    phi_if_surv = phi_if_surv,
    phi_if_fish = phi_if_fish,
    phi_ik = phi_ik,
    phi_ik2 = phi_ik2,
    phi_ij = phi_ij,
    phi_fm = phi_fm,
    tau_ik = tau_ik,
    
    #* DEMOG ----
    X_ija = X_ija,
    omega_ai = omega_ai,
    Linf_yk = growthPars$Linf_yk,
    kappa_yk = growthPars$kappa_yk,
    sigmaG_yk = growthPars$sigmaG_yk,
    
    #* DATA ----
    survey = survey, # Make sure the survey has the same length as the catch time series
    survey_err = survey_err, #ac.data$ss.error, # Make sure the survey has the same length as the catch time series
    survey = survey, #ac.data$ss.survey,
    age_error = ageerr_ExpAge,
    age_error_sd = ageerr_SD,
    L1_yk = L1_yk,
    catch = catch,
    
    #* ADDL PARS ----
    parms = parms,
    b = b,
    bfuture = bfuture,
    logSDR = logSDR
    # survey2 = survey2,
    # flag_surv_bio = survey_x2, #ac.data$survey_x, # Is there a survey in that year?
                  # wage_ssb = t(wage_ssb),
                  # wage_catch = t(wage_catch),
                  # wage_survey = t(wage_survey),
                  # wage_mid = t(wage_mid),
                  # selidx = which(years == selYear),
                  #  Input parameters
                  # year_sel = length(1991:max(years)), # Years to model time varying sel
                  # Msel = msel,
                  # Matsel= as.numeric(mat),
          
                  # logQ = log(1.14135),   # Analytical solution
                  # Selectivity 
                  # Smin = 1,
                  # Smin_survey = 2,
                  # Smax = 6,
                  # Smax_survey = 6,
                  # flag_sel = flag_sel,
                  # surveyseason = surveyseason,
                  # nsurvey = nsurvey, # Frequency of survey years (e.g., 2 is every second year)
                  # survey

                  # flag_surv_acomp =ac.data$sflag,
                  # acomp
                  # age_survey = age_survey.tmp,
                  # age_survey2 = age_survey.tmp2,
                  # age_maxage = 15, # Max age for age comps 
                  # Catch
                  #                Catchobs = catch$Fishery, # Convert to kg
                  # ss_catch = ac.data$ss.catch,
                  # flag_catch =ac.data$cflag,
                  # age_catch = age_catch.tmp,
                  # variance parameters
                  # logSDcatch = log(0.01),
           
                  # logphi_survey = log(11.46),


                  #logh = log(0.8),
                  # Space parameters 
                  # smul = 0.5, # Annual survey timing 
                  # sigma_psel = 1.4,
                  # sum_zero = 0,

                  #TAC = TAC,
                  # movemat = movemat,
                  # move = move,
                  # recruitmat = recruitmat,
                  # move.init = move.init,
                  # movefifty = movefifty,
                  # movemax = movemax,
                  # movesouth = movesouth,
                  # moveout = moveout,
                  # moveslope = moveslope,
                  # F0 = Fin,
                  # psel = psel,
                  # parms = parms,
                  # Fnseason = Fnseason,
                  # selectivity_change = selectivity_change,
                  # Catch2 = catch2,


                  # Parameters from the estimation model 
              
  )
  
 
  
  # Catch.country <- read.csv(here("input","data",'catch_per_country.csv'))
  # df$Catch.country <- as.matrix(Catch.country[,2:3])[,c(2,1)]
  
  # df$Catch <- rowSums(df$Catch.country)
  
  # if(nyear > length(df$Catch)){
  # 
  #   
  #   df$Catch <- c(df$Catch,rep(mean(df$Catch), nyear-length(df$Catch)))
  # 
  #   
  # }
  
  # if(nyear >nrow(df$Catch.country)){
  #   df$Catch.country <- rbind(df$Catch.country,t(replicate(nyear-nrow(Catch.country),colMeans(df$Catch.country))))
  # }
  
  if(yr_future > 0){
    
    idx.future <- length(1966:myear)+seq(2,yr_future, by = df$nsurvey) # Years where survey occurs 
    
    df$survey_x <- c(df$survey_x,rep(-2, yr_future))
    df$survey_x[idx.future] <- 2
    
    df$survey_err <- c(df$survey_err,rep(1, yr_future))
    df$survey_err[idx.future] <- mean(df$survey_err[df$survey_err != 1])
    
    df$ss_survey <- c(df$ss_survey, rep(0,  yr_future))
    df$ss_survey[idx.future] <- mean(df$ss_survey[df$ss_survey != -1])
    df$flag_survey <- c(df$flag_survey, rep(-1,yr_future))
    df$flag_survey[idx.future] <- 1
    df$flag_catch[years > 2018] <- 1
    
    Rdevs <- rnorm(n = yr_future,mean = 0, sd = exp(df$logSDR))
    #Rdevs <- rep(0, yr_future)
    df$parms$Rin <- c(df$parms$Rin,Rdevs)
    
    # Bias adjustment 
    df$b <- c(df$b,rep(df$bfuture, yr_future))
  }
  
  return(df) ## this has all the data in a format ready for estimation
  
}
