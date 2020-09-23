## SHIRE Operating model
## Spatial, transboundary data generator for sablefish
## Returns values ready to be estimated in runsabassessment.cpp
## teplate: run.agedbased.true.catch

runOM_datagen <- function(df, seed = 731){
  
  ## Load data & define structure ----
  nspace <- df$nspace
  nstocks <- df$nstocks
  nage <- df$nage
  age <- df$age
  
  nfleets_surv <- df$nfleets_surv
  nfleets_fish <- df$nfleets_fish
  nfleets_acomp <- df$nfleets_acomp
  nfleets_lcomp <- df$nfleets_lcomp
  
  fltnames_surv = df$fltnames_surv
  fltnames_fish = df$fltnames_fish
  fltnames_acomp = df$fltnames_acomp
  fltnames_lcomp = df$fltnames_lcomp
  selType_fish = df$selType_fish
  selType_surv = df$selType_surv
  
  nmgmt_reg <- ncol(df$phi_fm)
  
  nyear <- df$tEnd
  year <- df$years
  tEnd <- length(df$years)
  
  ## these are TMB-ready; need to add one for indexing to work
  phi_if_surv <- df$phi_if_surv 
  phi_if_fish <- df$phi_if_fish 
  phi_fm <- df$phi_fm
  phi_ik <- df$phi_ik 
  inames <- colnames(phi_ik)
  knames <- rownames(phi_ik)
  mnames <- colnames(phi_fm)
  phi_ik2 <- df$phi_ik2 + 1 ## zero-indexed, add one
  tau_ki <- df$tau_ki 
  
  ## Biology
  M_k <- df$M_k
  # M <- rep(0.2,nage) #M0*Msel # Naural mortality at age
  mat_age <- rep(0.15, nage) ## mortality
  wtatlen_kab <- df$wtatlen_kab
  mat_ak <- df$mat_ak ## maturity by age and stock
  load(here('input','input_data','unfished_ALK.rdata')) ## from prelim runs, for ssb0
  ## Obs
  # catch_yf_obs <- df$Catch2
  catch_yf_obs <- df$catch
  surv_yf_obs <- df$survey
  
  # // movement //
  omega_ais <- df$omega_ais
  # omega_ais[,] <- 0.5
  
  X_ijas <- df$X_ijas
  
  # // growth //
  Linf_yk <- df$Linf_yk
  L1_yk <- df$L1_yk
  kappa_yk <- df$kappa_yk
  sigmaG_yk <- df$sigmaG_yk
  phi_ij <- df$phi_ij ##  matrix of whether i,j are from distinct stocks (0 otherwise)
  LBins <- df$LBins
  
  
  # M, selectivity 
  ## note these combine age and length sel, flag type via selType
  fish_selex_yafs <- df$fish_selex_yafs
  surv_selex_yafs <- df$surv_selex_yafs
  

  SDR <- exp(df$logSDR)
  b <- rep(1, tEnd)
  q = 0.5 ## placeholder
  # True values 
  # M0 <- 0.2 #exp(df$parms$logMinit) # no difference between males and females
  h_k <- rep(0.7,4) # exp(df$parms$logh_k)
  R_0k <- rep(exp(df$parms$logRinit), nspace) ## change this to better value
  
  ## Unfished Naa and SB0 ----
  ## note that omega makes this non-smooth
  N_0ais <- array(0, dim = c(nage, nspace, 2), dimnames = list(c(age),c(inames),c('Fem','Mal')))
  for(s in 1:2){ ## 1 is female, 2 is male
    for(k in 1:nstocks){
      for(i in 1:nspace){
        for(a in 1:(nage-1)){
          N_0ais[a,i,s] = N_0ais[a,i,s]+
            0.5* 
            omega_ais[a,i,s]*
            R_0k[k]*
            tau_ki[k,i]* ## downscale to subarea. C1 and A1 higher because they're own stock
            exp(-(mat_age[a]*age[a]))
          # cat(paste(c(s,k,i,a), sep = " "),     N_0ais[a,i,s] ,"\n")
        } ## end age < plus
        # // note the A+ group will be in slot A-1
        N_0ais[nage,i,s] =   omega_ais[a,i,s]*
          N_0ais[nage-1,i,s]*
          exp(-(mat_age[nage-1]*age[nage-1]))/
          (1-exp(-mat_age[nage]*age[nage]))
      } ## // end subareas
    }  ## // end stocks
  } ## end sexes
  # N_0ais[1:15,,1] %>% data.frame() %>% mutate(rowSums(.))
  ## females only
  SSB_0i <- rep(0, nspace);  SSB_0k <- rep(0, nstocks);
  for(i in 1:nspace){
    for(a in 1:(nage)){
      SSB_0i[i] = SSB_0i[i] +  N_0ais[a,i,1]*wtatlen_kab[phi_ik2[i],1]*
        unfished_ALK_F[a,i]^wtatlen_kab[phi_ik2[i],2]*mat_ak[a,phi_ik2[i]] ;
      for(k in 1:nstocks){
        SSB_0k[k] = SSB_0k[k] + phi_ik[k,i]*SSB_0i[i];
      } #// end stocks
    } #// end ages
  } # // end space
  
  ## Ninit ----
  # Ninit_ais <- matrix(NA, nrow = nage, ncol = nspace)
  Ninit_ais <- array(0, dim = c(nage, nspace, 2), dimnames = list(c(age),c(inames),c('Fem','Mal')))
  tildeR_initk <-  rep(1, nstocks)
  tildeR_yk <- matrix(1, nrow = tEnd, ncol = nstocks)
  
  for(s in 1:2){
    for(k in 1:nstocks){   
      for(i in 1:nspace){
        for(a in 1:(nage-1)){
          Ninit_ais[a,i,s] = Ninit_ais[a,i,s]+
            0.5* ## sex ratio
            omega_ais[a,i,s] * ## omega is zero for non-movement years leading to weird shapes.
            tau_ki[k,i] * 
            R_0k[k]*
            exp(-(mat_age[a]*age[a])) #* 
            # exp(-0.5*SDR*SDR+tildeR_initk[k])
        } #// end ages
        Ninit_ais[nage,i,s] = omega_ais[a,i,s] * 
                                 exp(-mat_age[nage]*age[nage-1])/
          (1-exp(-(mat_age[nage]*age[nage])))#* 
                                                               # exp(-0.5*SDR*SDR+tildeR_initk[k]))
      } #// end space
    } #// end stocks
  } # end sex
  
  Length_yais_beg <- Length_yais_mid <- N_yais_beg <- N_yais_mid <-N_yais_end <- array(NA, dim = c(tEnd, nage, nspace, 2),
                                                                          dimnames = list(c(year),
                                                                                          c(age),
                                                                                          c(inames),
                                                                                          c('Fem','Mal')))
  LengthAge_alyis_beg <- LengthAge_alyis_mid <- array(NA,dim = c(nage,LBins,tEnd,nspace,2),
                                                      dimnames = list(c(age),c(1:LBins), c(year),
                                                                      c(inames),c('Fem','Mal')))
  
  SSB_yi <- matrix(0, nrow = tEnd, ncol = nspace); 
  SSB_yk <- matrix(0, nrow = tEnd, ncol = nstocks)
  R_yi <- matrix(0, nrow = tEnd, ncol = nspace)
  R_yk <- matrix(0, nrow = tEnd, ncol = nstocks)
  
  colnames(SSB_yi) <- colnames(R_yi) <- inames
  colnames(SSB_yk) <- colnames(R_yk) <- knames
  rownames(SSB_yi) <- rownames(SSB_yk) <- rownames(R_yi) <- rownames(R_yk) <- year
  
  
  niter <- 50 ## F iterations
  F1_yf <- F2_yf <- array(0, dim = c(tEnd, nfleets_fish, niter+1),
                          dimnames = list(c(year),
                                          c(fltnames_fish),
                                          c(1:(niter+1)))) ## storage for intermediate guesses
  Freal_yf <-  matrix(0, nrow = tEnd, ncol = nfleets_fish) ## storage for final guess
  Zreal_ya <- matrix(0, nrow = tEnd, ncol = nage) 
  Zreal_yai <- array(0, dim = c(tEnd, nage, nspace))
  F_area_yfi <-  array(0, dim = c(tEnd, nfleets_fish, nspace), dimnames = list(c(year),
                                                                               c(fltnames_fish),
                                                                               c(inames)))
  F_ym <- matrix(0, nrow = tEnd, ncol = nmgmt_reg)
  # Catch_yf_est <- array(0, dim = c(tEnd,nage,nfleets_fish))
  # CatchN <- matrix(0, nrow = tEnd, ncol = nfleets_fish)
  
  catch_yaf_pred <-  array(0, dim = c(tEnd, nage, nfleets_fish),
                           dimnames = list(c(year),
                                           c(age),
                                           c(paste(fltnames_fish))))
  catch_yaif_pred <-  array(0, dim = c(tEnd, nage, nspace, nfleets_fish),
                            dimnames = list(c(year),
                                            c(age),
                                            c(inames),
                                            c(paste(fltnames_fish))))
  catch_yf_pred <- N_avail_yf <-  matrix(0, nrow= tEnd, ncol = nfleets_fish,
                                         dimnames = list(c(year), paste(fltnames_fish)))
  N_weight_yfi <- catch_yfi_pred <-  array(0, dim = c(tEnd, nfleets_fish, nspace),
                                           dimnames = list(c(year), paste(fltnames_fish), c(inames)))
  Nsamp_acomp_yf <-  survey_yf_pred <- matrix(0, nrow= tEnd, ncol = nfleets_surv,
                                              dimnames = list(c(year), paste(fltnames_surv)))
  
  ## start year loop ----
  # for(y in 1:25){
  for(y in 1:(tEnd-1)){
    cat(y,"\n")
    ## Year 0 ----
    if(y == 1){
      # for(k in 1:nstocks){   
      for(i in 1:nspace){
        for(s in 1:2){ ## sexes
          Length_yais_beg[1,1,i,s] <- L1_yk[y,phi_ik2[i],s]
          N_yais_beg[1,1,i,s] <- Ninit_ais[1,i,s]
          N_yais_mid[1,1,i,s] <- N_yais_beg[1,1,i,s]*exp(-mat_age[1]/3)
          for(a in 2:(nage-1)){ ## fill A0 in position 1 later
            for(j in 1:nspace){           
              pLeave = 0.0;  NCome = 0.0; # // reset for new age
              if(i != j){
                pLeave = pLeave + X_ijas[i,j,a,s]; #// will do 1-this for proportion which stay
                NCome = NCome + X_ijas[j,i,a,s]*Ninit_ais[a,j,s]; #// actual numbers incoming
              }
            } #// end subareas j
            # // this is the synthesis syntax; 10 is placeholder for LMIN
            # // likely need a lower L1 at age stock-specific and linear before that age
            
            Length_yais_beg[y,a,i,s] = Linf_yk[1,phi_ik2[i],s]+(L1_yk[y,phi_ik2[i],s]-Linf_yk[1,phi_ik2[i],s])*
              exp(-kappa_yk[1,phi_ik2[i],s]*a)
            Length_yais_mid[y,a,i,s] = Linf_yk[1,phi_ik2[i],s]+(L1_yk[y,phi_ik2[i],s]-Linf_yk[1,phi_ik2[i],s])*
              exp(-0.5*kappa_yk[1,phi_ik2[i],s]*a)
            N_yais_beg[y,a,i,s] = ((1-pLeave)*Ninit_ais[a,i,s] + NCome)*exp(-mat_age[a]/3)
          } #// end ages
          ## // plus group includes those already at A AND age into A
          for(j in 1:nspace){   
            pLeave = 0.0;  NCome = 0.0; # // reset for new age
            if(i != j){
              pLeave = pLeave + X_ijas[i,j,nage,s]
              NCome = NCome + X_ijas[j,i,nage,s]*(Ninit_ais[nage,j,s] + Ninit_ais[nage-1,j,s])  #// if M becomes spatial use M_aj here
            }
          } #// end subareas j
          N_yais_beg[y,nage,i,s] =  ((1-pLeave)*(Ninit_ais[nage,i,s] + Ninit_ais[nage-1,i,s]) +  NCome)*exp(-mat_age[nage]/3)
          
          Length_yais_beg[y,nage,i,s] = Linf_yk[1,phi_ik2[i],s]+(L1_yk[y,phi_ik2[i],s]-Linf_yk[1,phi_ik2[i],s])*
            exp(-kappa_yk[1,phi_ik2[i],s]*nage-1)
          Length_yais_mid[y,nage,i,s]  = Linf_yk[1,phi_ik2[i],s]+(L1_yk[y,phi_ik2[i],s]-Linf_yk[1,phi_ik2[i],s])*
            exp(-0.5*kappa_yk[1,phi_ik2[i],s]*nage-1)
          
        } #// end sexes
      } #// end  subareas i
    } ## end y == 1
    # if(any(is.na(N_yais_beg[y,,,]))) stop('NA ON year', y,"\n")
    ## SSB_y ----
    
    for(i in 1:nspace){
      SSB_yi[y,i] <- 0
      for(a in 1:(nage)){
        SSB_yi[y,i] <- SSB_yi[y,i] +  N_yais_beg[y,a,i,1]*wtatlen_kab[phi_ik2[i],1]*
          Length_yais_beg[y,a,i,1]^wtatlen_kab[phi_ik2[i],2]*mat_ak[a,phi_ik2[i]]
        if(is.na(SSB_yi[y,i])) stop("NA ON SSB_YI year ",y," space ",i," age ",a)
      } #// end ages
    } #// end space
    for(k in 1:nstocks){
      SSB_yk[y,k]<- 0
      for(i in 1:nspace){
        SSB_yk[y,k] <- SSB_yk[y,k] + phi_ik[k,i]*SSB_yi[y,i] 
      } # // end stocks
    } #// end space
    # cat(sum(SSB_yk[y,]),"\n")
    # cat(sum(SSB_yi[y,]),"\n")
    if(is.na(sum(SSB_yk[y,]))) stop("NA ON SSB_YK",y) 
    
    ## Ryi, Ryk Recruits ----
    # next year based on present SSB
    omega_0ij <- rep(1, nspace)
    for(i in 1:nspace){
      for(k in 1:nstocks){
        # // SSB_yk already has summation
        R_yk[y,k] = (4*h_k[k]*R_0k[k]*SSB_yk[y,k])/
          (SSB_0k[k]*(1-h_k[k])+ 
             SSB_yk[y,k]*(5*h_k[k]-1))#*exp(-0.5*b[y]*SDR*SDR+tildeR_yk[y,k])
        # if(R_yk[y,k] == 0) stop(paste("RYK IS ZER ON,",y,k,"\n"))
      } # // end stocks
      R_yi[y,i] = R_yk[y,phi_ik2[i]]*tau_ki[phi_ik2[i],i]*omega_0ij[i] #// downscale to subarea including age-0 movement
      N_yais_beg[y+1,1,i,1:2] = 0.5*R_yi[y,i] #// fill age-0 recruits
    } ### end space
    # cat(sum(R_yk[y,]),"\n")
    # cat(sum(R_yi[y,]),"\n")
    # cat(sum(N_yais_beg[y+1,1,,]),"\n")
    
    #N- and Nominal Length ----
    # at-age for the middle of this year and beginning of next 
    for(s in 1:2){
      for(i in 1:nspace){
        N_yais_mid[y,1,i,s] <- N_yais_beg[y,1,i,s]*exp(-mat_age[1]/3)
        # Length_yais_beg[y,1,i,s] <- L1_yk[y,phi_ik2[i],s]
        ## linear growth below A4 as in synthesis
        len.step <- ifelse( L1_yk[y,phi_ik2[i],s] < 4,  L1_yk[y,phi_ik2[i],s], 4)  ## lmin is size at age 0 (?)
        len.slope <- ( L1_yk[y,phi_ik2[i],s] - len.step) / 4 ## linear growth slope
        # Length at the start of the year (cm)
        Length_yais_beg[y,1:4,i,s] <- len.step[1]+len.slope*(seq(1, (4 + 1), 1) - 1)[1:4]  
        Length_yais_beg[y,5,i,s] <- L1_yk[y,phi_ik2[i],s] ## L1 Corresponds to age 4 per analysis
        Length_yais_mid[y,1:5,i,s] <- Length_yais_beg[y,1:5,i,s] + (Linf_yk[y,phi_ik2[i],s]-Length_yais_beg[y,1:5,i,s]*
                                                                      (1-exp(-0.5*kappa_yk[y,phi_ik2[i],s])))
        for(a in 2:(nage-1)){ ## note that TMB starts at pos 1 which is age 1 which is pos 2 here
          pLeave = 0.0;  NCome = 0.0
          for(j in 1:nspace){           
            if(i != j){
              pLeave = pLeave + X_ijas[i,j,a,s]; ### will do 1-this for proportion which stay
              NCome = NCome + X_ijas[j,i,a,s]*N_yais_beg[y,a,j,s]; ### actual numbers incoming
            }
          } ### end subareas j         
          N_yais_mid[y,a,i,s] = N_yais_beg[y,a,i,s]*exp(-mat_age[a]/3) 
          N_yais_beg[y+1,a,i,s] = ((1-pLeave)*N_yais_beg[y,a-1,i,s] + NCome)*exp(-mat_age[a]/3) ## this exponent needs to be Ztuned eventually
        } ## end ages for N

        for(a in 6:nage-1){
          ## as in document: next year A1 == this year A0 plus growth
          Length_yais_beg[y+1,a,i,s] = Length_yais_beg[y,a-1,i,s] + (Linf_yk[y,phi_ik2[i],s]-Length_yais_beg[y,a-1,i,s])*
            (1-exp(-kappa_yk[y,phi_ik2[i],s]))
          Length_yais_mid[y,a,i,s] = Length_yais_beg[y,a,i,s] + (Linf_yk[y,phi_ik2[i],s]-Length_yais_beg[y,a,i,s]*
                                                                   (1-exp(-0.5*kappa_yk[y,phi_ik2[i],s])))
        } ## end ages for L
        # ## plus groups
        pLeave = 0.0;  NCome = 0.0
        for(j in 1:nspace){           
          if(i != j){
            pLeave <- pLeave + X_ijas[i,j,nage-1,s] 
            NCome <- NCome + X_ijas[j,i,nage-1,s]*(N_yais_beg[y,nage,j,s] + N_yais_beg[y,nage-1,j,s]) 
          } ## end i != j
        } ## end subareas j
        N_yais_mid[y,nage,i,s] = N_yais_beg[y,nage,i,s]*exp(-mat_age[nage]/3)
        N_yais_beg[y+1,nage,i,s] =   ((1-pLeave)*( N_yais_beg[y,nage,i,s]+ N_yais_beg[y,nage-1,i,s]) + NCome)*exp(-mat_age[nage]/3);
        ## plus group weighted average (we already have the numbers at age)
        Length_yais_beg[y+1,nage,i,s] = ( N_yais_beg[y+1,nage-1,i,s]*
                                            (Length_yais_beg[y,nage-1,i,s]+
                                               (Linf_yk[y,phi_ik2[i],s]-Length_yais_beg[y,nage-1,i,s]*(1-exp(-kappa_yk[y,phi_ik2[i],s])))) +
                                            N_yais_beg[y+1,nage-1,i,s]*
                                            (Length_yais_beg[y,nage,i,s]+
                                               (Linf_yk[y,phi_ik2[i],s]-Length_yais_beg[y,nage,i,s])*(1-exp(-kappa_yk[y,phi_ik2[i],s]))))/
          (N_yais_beg[y+1,nage-1,i,s] + N_yais_beg[y+1,nage,i,s])
        
        Length_yais_mid[y+1,nage,i,s] = (N_yais_mid[y,nage-1,i,s]*
                                           (Length_yais_beg[y,nage-1,i,s]+(Linf_yk[y,phi_ik2[i],s]-Length_yais_beg[y,nage-1,i,s]*(1-exp(-0.5*kappa_yk[y,phi_ik2[i],s])))) +
                                           N_yais_mid[y,nage,i,s]*
                                           (Length_yais_beg[y,nage,i,s]+(Linf_yk[y,phi_ik2[i],s]-Length_yais_beg[y,nage,i,s])*(1-exp(-0.5*kappa_yk[y,phi_ik2[i],s]))))/
          (N_yais_mid[y,nage-1,i,s] + N_yais_mid[y,nage,i,s])
      } ## end subareas i
      # cat(y+1,sum(Length_yais_mid[y+1,,,s]),"\n")
    } ## end sexes
    
    ## reweight length-at-age based on movement from other stocks ----
    for(s in 1:2){
      for(i in 1:nspace){  
        for(a in 1:(nage)){ ## note that TMB starts at pos 1 which is age 1 which is pos 2 here
          LCome = 0.0; NCome = 0.0
          for(j in 1:nspace){           
            if(i != j){
              LCome = LCome + phi_ij[i,j]*N_yais_beg[y,a,j,s]*Length_yais_beg[y,a,j,s] ## for numerator
              NCome = NCome + phi_ij[i,j]*N_yais_beg[y,a,j,s] ## for denom
            }
          } ## end subareas j
          Length_yais_beg[y+1,a,i,s] = (N_yais_beg[y,a,i,s]*Length_yais_beg[y,a,i,s] + LCome)/
            (N_yais_beg[y,a,i,s]+NCome)
          # if(is.na(Length_yais_beg[y+1,1,i,s])) stop("NA Length_yais_beg on year ",y, " age1 ", "space ",i)
        } ## end ages
      } ## end subareas i
    } ## sexes
    
    ## prob of length-at-age ----
    for(s in 1:2){
      for(i in 1:nspace){  
        for(a in 1:(nage)){
          LengthAge_alyis_beg[a,1,y,i,s] = pnorm(1,  Length_yais_beg[y,a,i,s], sigmaG_yk[y,phi_ik2[i],s]);
          LengthAge_alyis_mid[a,1,y,i,s] = pnorm(1,  Length_yais_mid[y,a,i,s], sigmaG_yk[y,phi_ik2[i],s]);
          for(l in 2:(LBins-1)){
            LengthAge_alyis_beg[a,l,y,i,s] = pnorm(l+1,  Length_yais_beg[y,a,i,s], sigmaG_yk[y,phi_ik2[i],s]) -
              pnorm(l,  Length_yais_beg[y,a,i,s], sigmaG_yk[y,phi_ik2[i],s])
            LengthAge_alyis_mid[a,l,y,i,s] = pnorm(l+1,  Length_yais_mid[y,a,i,s], sigmaG_yk[y,phi_ik2[i],s]) -
              pnorm(l,  Length_yais_mid[y,a,i,s], sigmaG_yk[y,phi_ik2[i],s])
          } ## end LBins
          LengthAge_alyis_beg[a,LBins,y,i,s] = 1-pnorm(LBins, Length_yais_beg[y,a,i,s], sigmaG_yk[y,phi_ik2[i],s]);
          LengthAge_alyis_mid[a,LBins,y,i,s] = 1-pnorm(LBins, Length_yais_mid[y,a,i,s], sigmaG_yk[y,phi_ik2[i],s]);
          if(is.na(   LengthAge_alyis_beg[a,LBins,y,i,s])) stop('NA ON ', a,l,y,i,s,"\n")
          
        } ## end ages
      } ## end nspace
    } ## end sex
    
    }  ## end years testing
    # for(y in 15:25){
  for(y in 1:(tEnd-1)){
    ## Hybrid F tuning  ----
    # v1 <- 0.99;   Fmax <- 3; ##corresponds to an Fmax of 3
    # v1 = 0.865; Fmax <- 2##corresponds to an Fmax of 2
    v1 = 0.7; Fmax = 1.5
    # v1 = 0.65; Fmax <- 1.15
    # v1 = 0.25
    v2 <- 30;
    
    ##storage for intermediate guesses, reset each year
    catch_afk_TEMP <- array(0, dim = c(nage, nfleets_fish, niter+1), 
                            dimnames = list(c(age),
                                            c(paste(fltnames_fish)),
                                            c(1:(niter+1)))) 
    
    Adj <- Z_a_TEMP <- Z_a_TEMP2 <- NULL
    for(fish_flt in 2){
    # for(fish_flt in 1:nfleets_fish){
      
      selMult = c(0.5,0.3,1,1,1/4,1/4,1/10,1,1)[fish_flt]
      # selMult <- ifelse(fish_flt %in% c(1,2),0.5,ifelse(fish_flt %in% 5:7,0.25,1)) ## trying to tweak selex to fit
      
    # for(fish_flt in 1:nfleets_fish){
      catch_yaf_pred[y,,fish_flt] <- catch_yf_pred[y,fish_flt] <- catch_yfi_pred[y,fish_flt,] <-
        catch_yaif_pred[y,,,fish_flt] <- 0
      if(is.na(catch_yf_obs[y, fish_flt+1])) next() ## skip if no catch
      ## putative biomass available
      denom = 0
      for(i in 1:nspace){
        # if(phi_imat_age[i, m] == 0) next() ## skip area if not in mgmt reg
        ## this needs to deal accurately with sex-selex length OR age
        if(selType_fish[fish_flt] == 'AGE'){
          denom <- denom + (phi_if_fish[fish_flt, i] *
                              sum(selMult*fish_selex_yafs[y,,fish_flt,]*N_yais_beg[y,,i,]*
                                    wtatlen_kab[phi_ik2[i],1]*
                                    Length_yais_beg[y,,i,]^wtatlen_kab[phi_ik2[i],2],
                              catch_yf_obs[y, fish_flt+1]))
        } else if(selType_fish[fish_flt] == 'LEN'){
          # per AEP 
          denom <- denom + (phi_if_fish[fish_flt, i] *
                              ## not sure about this; M and F n_age x most likely len at age into weight = biomass??
                              sum(
                                # selMult*fish_selex_yafs[y, , fish_flt, 1] * 
                                  N_yais_beg[y, , i, 1] *
                                  LengthAge_alyis_beg[, , y, i, 1] *
                                  wtatlen_kab[phi_ik2[i], 1] *
                                  apply(LengthAge_alyis_beg[, , y, i, 1],1,which.max)^ ## most likely lengths at age ^
                                              wtatlen_kab[phi_ik2[i], 2],
                                # selMult*fish_selex_yafs[y, , fish_flt, 2] * 
                                  N_yais_beg[y, , i, 2] *
                                  LengthAge_alyis_beg[, , y, i, 1] *
                                  wtatlen_kab[phi_ik2[i], 1] *
                                  apply(LengthAge_alyis_beg[, , y, i, 2],1,which.max)^
                                              wtatlen_kab[phi_ik2[i], 2],
                              catch_yf_obs[y, fish_flt + 1]))
          
          ## old way
          # denom <- denom + (phi_if_fish[fish_flt, i] *
          #                     ## not sure about this; M and F n_age x most likely len at age into weight = biomass??
          #                     sum(
          #                       selMult*fish_selex_yafs[y, , fish_flt, 1] * N_yais_beg[y, , i, 1] *
          #                         wtatlen_kab[phi_ik2[i], 1] *
          #                         apply(LengthAge_alyis_beg[, , y, i, 1],1,which.max)^ ## most likely lengths at age ^
          #                                     wtatlen_kab[phi_ik2[i], 2],
          #                       selMult*fish_selex_yafs[y, , fish_flt, 2] * N_yais_beg[y, , i, 2] *
          #                         wtatlen_kab[phi_ik2[i], 1] *
          #                         apply(LengthAge_alyis_beg[, , y, i, 2],1,which.max)^ 
          #                                     wtatlen_kab[phi_ik2[i], 2],
          #                     catch_yf_obs[y, fish_flt + 1]))
     
        } ## end len sel
        # cat(i,denom,"\n")
      }
      ## make an initial guess for Ff using obs catch - need to update selex whihc is 1.0 now
      ## make this guess by M, and sum over phi_im
      F1_yf[y, fish_flt, 1] <-    catch_yf_obs[y, fish_flt+1]/denom
      latest_guess <-    F1_yf[y, fish_flt, 1]
  
      ## k iterations ----
      for(k in 2:(niter+1)){
        catch_afk_TEMP[,fish_flt,k] <- 0
        
        ## modify the guess Eq 20
        term0 = 1/(1+exp(v2*( latest_guess - v1)));
        term1 = latest_guess*term0;
        term2 = v1*(1-term0);
        F1_yf[y,fish_flt,k] = -log(1-(term1+term2))

        ## note that the k in catch_afk is ITERS, not stocks
        for(i in 1:nspace){
          for(a in 1:nage){
            Z_a_TEMP[a] <- sum(fish_selex_yafs[y, a, fish_flt, ]*F1_yf[y,fish_flt,k]) + mat_age[a]
            if(selType_fish[fish_flt] == 'AGE'){
              
              catch_afk_TEMP[a,fish_flt,k] <-    catch_afk_TEMP[a,fish_flt,k] +
                (F1_yf[y,fish_flt,k]/(Z_a_TEMP[a]))*
                (1-exp(-Z_a_TEMP[a]))*
                phi_if_fish[fish_flt, i]*
                sum(selMult*fish_selex_yafs[y,a,fish_flt,]* N_yais_beg[y,a,i,]*
                      wtatlen_kab[phi_ik2[i],1]*
                      Length_yais_beg[y,a,i,]^wtatlen_kab[phi_ik2[i],2])
              
              
            } else if(selType_fish[fish_flt] == 'LEN'){
              LAA <- ifelse( which.max(LengthAge_alyis_beg[a, , y, i, 1]) > length(fish_selex_yafs[y,      , fish_flt, 1]),
                             length(fish_selex_yafs[y,    nage  , fish_flt, 1]),
                             which.max(LengthAge_alyis_beg[a, , y, i, 1]))
              
              # per AEP
              catch_afk_TEMP[a,fish_flt,k] <-    catch_afk_TEMP[a,fish_flt,k] +
                (F1_yf[y,fish_flt,k]/(Z_a_TEMP[a]))*
                (1-exp(-Z_a_TEMP[a]))*
                phi_if_fish[fish_flt, i]*
                sum(
                # selMult*fish_selex_yafs[y,LAA, fish_flt, 1] *
                       N_yais_beg[y, a, i, 1] *
                       LengthAge_alyis_beg[a, , y, i, 1] *
                       wtatlen_kab[phi_ik2[i], 1] *
                       LAA ^
                       wtatlen_kab[phi_ik2[i], 2],
                     # selMult*fish_selex_yafs[y, LAA, fish_flt, 2] *
                       N_yais_beg[y, a, i, 2] *
                       LengthAge_alyis_beg[a, , y, i, 1] *
                       wtatlen_kab[phi_ik2[i], 1] *
                       LAA ^
                       wtatlen_kab[phi_ik2[i], 2] )
              
              # catch_afk_TEMP[a,fish_flt,k] <-    catch_afk_TEMP[a,fish_flt,k] +
              #   (F1_yf[y,fish_flt,k]/(Z_a_TEMP[a]))*
              #   (1-exp(-Z_a_TEMP[a]))*
              #   phi_if_fish[fish_flt, i]*
              #   sum( selMult*fish_selex_yafs[y,LAA, fish_flt, 1] *
              #       N_yais_beg[y, a, i, 1] *
              #       wtatlen_kab[phi_ik2[i], 1] *
              #       which.max(LengthAge_alyis_beg[a, , y, i, 1]) ^
              #       wtatlen_kab[phi_ik2[i], 2],
              #       selMult*fish_selex_yafs[y, LAA, fish_flt, 2] *
              #       N_yais_beg[y, a, i, 2] *
              #       wtatlen_kab[phi_ik2[i], 1] *
              #       which.max(LengthAge_alyis_beg[a, , y, i, 2]) ^
              #       wtatlen_kab[phi_ik2[i], 2] )
            } ## end seltype Len
          } ## end ages
        } ## end space
        # cat(sum(catch_afk_TEMP[,fish_flt,k]),"\n")
        ## Calc Adj Eq 22
        Adj[k] <- catch_yf_obs[y, fish_flt+1]/sum(catch_afk_TEMP[,fish_flt,k])
        # cat(k," ",Adj[k],"\n")
        ## Get new Z given ADJ - need to add discard here
        for (a in 1:nage)   {
          Z_a_TEMP2[a] <- Adj[k] * 
            sum(fish_selex_yafs[y, a, fish_flt, ] * F1_yf[y, fish_flt, k]) +
            mat_age[a]
        }
        
        ## Second Guess for F (EQ 24)
        denom = 0
        for(i in 1:nspace){
          for(a in 1:nage){
            if(selType_fish[fish_flt] == 'AGE'){
              denom <- denom + phi_if_fish[fish_flt, i] *
                sum(selMult*fish_selex_yafs[y,a,fish_flt,]* N_yais_beg[y,a,i,]*
                      wtatlen_kab[phi_ik2[i],1]*
                      Length_yais_beg[y,a,i,]^wtatlen_kab[phi_ik2[i],2])*
                (1-exp(-Z_a_TEMP2[a])) * (F1_yf[y,fish_flt,k]/(Z_a_TEMP2[a]))
            } else if(selType_fish[fish_flt] == 'LEN'){
              LAA <- ifelse( which.max(LengthAge_alyis_beg[a, , y, i, 1]) > length(fish_selex_yafs[y,      , fish_flt, 1]),
                             length(fish_selex_yafs[y,    nage  , fish_flt, 1]),
                             which.max(LengthAge_alyis_beg[a, , y, i, 1]))
              
              denom <- denom + phi_if_fish[fish_flt, i] *
                sum( 
                  # selMult*fish_selex_yafs[y,LAA, fish_flt, 1] *
                       N_yais_beg[y, a, i, 1] *
                       LengthAge_alyis_beg[a, , y, i, 1] *
                       wtatlen_kab[phi_ik2[i], 1] *
                       LAA ^
                       wtatlen_kab[phi_ik2[i], 2],
                     # selMult*fish_selex_yafs[y, LAA, fish_flt, 2] *
                       N_yais_beg[y, a, i, 2] *
                       LengthAge_alyis_beg[a, , y, i, 1] *
                       wtatlen_kab[phi_ik2[i], 1] *
                       LAA ^
                       wtatlen_kab[phi_ik2[i], 2] )*
                (1-exp(-Z_a_TEMP2[a])) * (F1_yf[y,fish_flt,k]/(Z_a_TEMP2[a]))
            } ## end seltype Len
          } ## end age
          # cat(i,denom,"\n")
        } ## end space
        # F2_yf[y, fish_flt, k] <- F2_yf[y, fish_flt, k-1] + catch_yf_obs[y, fish_flt+1]/denom
        ## this is the ratio between catch and putative biomass given Z update
        F2_yf[y, fish_flt, k] <- catch_yf_obs[y, fish_flt+1]/denom
        
        ## Modify the guess again Eq 25
        # term0 = 1/(1+exp(v2*( F2_yf[y,fish_flt,k] - v1)));
        term0 = 1/(1+exp(v2*( F2_yf[y,fish_flt,k] - v1*Fmax)));
        term1 = F2_yf[y,fish_flt,k]*term0;
        term2 = v1*(1-term0);
        F2_yf[y,fish_flt,k] = -log(1-(term1+term2))
        # cat(F2_yf[y,fish_flt,k],"\n")
        latest_guess <- F2_yf[y,fish_flt,k]
        # cat(k, latest_guess,"\n")
      } ## end hybrid F iterations
      
      ## Define F, Z and predicted catches ----
      Freal_yf[y, fish_flt] <- latest_guess ## final as Freal_yf
      if(is.na(Freal_yf[y, fish_flt])) stop("NA F ON",y,"\t", fish_flt,"\n")
      
      ## annoying multi-loops for F in area
      # N_avail_yf[y,fish_flt] <- 0
      ## get total N exploitable by this fleet
      # for(i in 1:nspace){
      #   N_avail_yf[y,fish_flt] <- N_avail_yf[y,fish_flt] + sum( phi_if_fish[fish_flt, i]*
      #                                                             sum(N_yais_beg[y,,i,]))
      # }
      ## get ratio of N in area & reweight F
      ## will just return Freal and 1 for single-area fisheries
      # for(i in 1:nspace){
      #   N_weight_yfi[y,fish_flt, i] <- sum(phi_if_fish[fish_flt, i]* sum(N_yais_beg[y,,i,]))/ N_avail_yf[y,fish_flt]
      #   F_area_yfi[y,fish_flt,i] <- Freal_yf[y, fish_flt]*N_weight_yfi[y,fish_flt, i]
      # }
      
      ## add together for mgmt regions
      for(m in 1:nmgmt_reg){
        F_ym[y,m] <- F_ym[y,m]+phi_fm[fish_flt,m]*Freal_yf[y, fish_flt]
      }
   
      for(i in 1:nspace){
        for(a in 1:nage){
          if(selType_fish[fish_flt] == 'AGE'){
            Zreal_ya[y,a] <-   Freal_yf[y, fish_flt] + mat_age[a] ## should this include all fleets?
            
            catch_yaf_pred[y,a,fish_flt] <- catch_yaf_pred[y,a,fish_flt] +
              (Freal_yf[y, fish_flt]/(Zreal_ya[y,a]))*(1-exp(-Zreal_ya[y,a]))*
              phi_if_fish[fish_flt, i]*
              
              sum(selMult*fish_selex_yafs[y,a,fish_flt,]*N_yais_beg[y,a,i,]*
                    wtatlen_kab[phi_ik2[i],1]*
                    Length_yais_beg[y,a,i,]^wtatlen_kab[phi_ik2[i],2])
            
            Zreal_yai[y,a,i]  <-  F_area_yfi[y,fish_flt,i] + mat_age[a]

            catch_yaif_pred[y,a,i,fish_flt] <- (F_area_yfi[y,fish_flt,i]/
                                                  (  Zreal_yai[y,a,i] ))*(1-exp(-  Zreal_yai[y,a,i] ))*
              phi_if_fish[fish_flt, i]*

              sum(fish_selex_yafs[y,a,fish_flt,]* N_yais_beg[y,a,i,]*
                    wtatlen_kab[phi_ik2[i],1]*
                    Length_yais_beg[y,a,i,]^wtatlen_kab[phi_ik2[i],2])
            
          } else if(selType_fish[fish_flt] == 'LEN'){
            ## if the expected length at age is greater than we have selex for, just use the last value
            LAA <- ifelse( which.max(LengthAge_alyis_beg[a, , y, i, 1]) > length(fish_selex_yafs[y,      , fish_flt, 1]),
                           length(fish_selex_yafs[y,    nage  , fish_flt, 1]),
                           which.max(LengthAge_alyis_beg[a, , y, i, 1])
                           )
            
            Zreal_ya[y,a] <-   Freal_yf[y, fish_flt] + mat_age[a] ## should this include all flets?
            
            catch_yaf_pred[y,a,fish_flt] <- catch_yaf_pred[y,a,fish_flt] +
              (Freal_yf[y, fish_flt]/(Zreal_ya[y,a]))*(1-exp(-Zreal_ya[y,a]))*
              phi_if_fish[fish_flt, i]*
              sum( 
                # selMult*fish_selex_yafs[y,LAA, fish_flt, 1] *
                     N_yais_beg[y, a, i, 1] *
                     LengthAge_alyis_beg[a, , y, i, 1] *
                     wtatlen_kab[phi_ik2[i], 1] *
                     LAA ^
                     wtatlen_kab[phi_ik2[i], 2],
                   # selMult*fish_selex_yafs[y, LAA, fish_flt, 2] *
                     N_yais_beg[y, a, i, 2] *
                     LengthAge_alyis_beg[a, , y, i, 1] *
                     wtatlen_kab[phi_ik2[i], 1] *
                     LAA ^
                     wtatlen_kab[phi_ik2[i], 2] )
            
            Zreal_yai[y,a,i]  <-  F_area_yfi[y,fish_flt,i] + mat_age[a]/3

            catch_yaif_pred[y,a,i,fish_flt] <-       catch_yaif_pred[y,a,i,fish_flt] +(F_area_yfi[y,fish_flt,i]/
                                                  (  Zreal_yai[y,a,i] ))*(1-exp(-  Zreal_yai[y,a,i] ))*
              phi_if_fish[fish_flt, i]*
              sum(
                # fish_selex_yafs[y,LAA, fish_flt, 1] *
                  N_yais_beg[y, a, i, 1] *
                  wtatlen_kab[phi_ik2[i], 1] *
                  which.max(LengthAge_alyis_beg[a, , y, i, 1]) ^
                  wtatlen_kab[phi_ik2[i], 2],
                # fish_selex_yafs[y, LAA, fish_flt, 2] *
                  N_yais_beg[y, a, i, 2] *
                  wtatlen_kab[phi_ik2[i], 1] *
                  which.max(LengthAge_alyis_beg[a, , y, i, 2]) ^
                  wtatlen_kab[phi_ik2[i], 2] )
            
          } ## end seltype len
        } ## end ages for predicted catch
        catch_yfi_pred[y,fish_flt,i] <- sum(catch_yaif_pred[y,,i,fish_flt])
      } ## end nspace for predicted catch
      catch_yf_pred[y,fish_flt] <- sum(catch_yaf_pred[y,,fish_flt])
      # if(catch_yf_pred[y,fish_flt]  > 2*catch_yf_obs[y,fish_flt+1] &
      #    catch_yf_obs[y,fish_flt+1] != 0 ) stop('cpred huge on flt ',  fltnames_fish[fish_flt], 'year ',y)
      # if(catch_yf_pred[y,fish_flt]  < catch_yf_obs[y,fish_flt+1]/2 &
      #    catch_yf_obs[y,fish_flt+1] != 0 ) stop('cpred tiny on flt ',  fltnames_fish[fish_flt], 'year ',y)
      cat(  y," PRED ",catch_yf_pred[y,fish_flt]," OBS ",catch_yf_obs[y,fish_flt+1]," ",fish_flt,
            " RATIO ", catch_yf_pred[y,fish_flt]/catch_yf_obs[y,fish_flt+1],"\n")
    } ## end fishery fleets

    ## N_yais_end ----
    ## now update N_yais_end (terminal post-fishing biomass)
    #N- and Nominal Length ----
    # at-age for the middle of this year and beginning of next 
    for(s in 1:2){
      for(i in 1:nspace){
        for(a in 1:nage){
        ## Z real has F and first third of biomass, mid has second third
        ## assumes no more movement
        N_yais_end[y,a,i,s] <- N_yais_mid[y,a,i,s]*exp(-mat_age[a]/3+Zreal_yai[y,a,i])
        } ## end ages for L
      } ## end subareas i
      # cat(y+1,sum(Length_yais_mid[y+1,,,s]),"\n")
    } ## end sexes
    
    
    # } ## yr test
    # head(catch_yf_pred)
    
    ## survey biomass ----
    ## Estimate survey biomass at midyear 
    
  
  
  for(y in 1:(tEnd-1)){   
    for( sur_flt in 1:nfleets_surv){
      if(is.na(surv_yf_obs[y,sur_flt])) next()
      Nsamp_acomp_yf[y,sur_flt] <- survey_yf_pred[y,sur_flt] <- 0
      # if(is.na(surv_yf_obs[y,sur_flt])){  next(cat('skipping',y, fltnames_surv,"\n"))}
      # selMult <- c(1/4,1/22,1/25,1/5,0.9)[sur_flt]
      selMult <- rep(1,nfleets_surv)[sur_flt]
      
      for(i in 1:nspace){
        for(a in 1:nage){
          # if(selType_surv[sur_flt] == 'AGE'){
          survey_yf_pred[y,sur_flt] <-  survey_yf_pred[y,sur_flt] +
            q*
            phi_if_surv[sur_flt,i]*
            sum(
                selMult*surv_selex_yafs[y,a,sur_flt,]*
                    N_yais_mid[y,a,i,]*
                  wtatlen_kab[phi_ik2[i],1]*
                  Length_yais_mid[y,a,i,]^wtatlen_kab[phi_ik2[i],2])

          Nsamp_acomp_yf[y,sur_flt] <-
            Nsamp_acomp_yf[sur_flt]  +
           phi_if_surv[sur_flt,i]* sum(selMult*surv_selex_yafs[y,a,sur_flt,]*N_yais_mid[y,a,i,]); ## To use with age comps; may need to change phi to sum acomp surveys
        } ## end ages
      } ## end nspace
      
      cat(  y," PRED ",survey_yf_pred[y,sur_flt]," OBS ",
            surv_yf_obs[y,sur_flt]," ",sur_flt,' RATIO ', 
            survey_yf_pred[y,sur_flt]/surv_yf_obs[y,sur_flt],  "\n")
    } ## end surv fleets
        # } else if(selType_surv[sur_flt] == 'LEN'){
        #     LAA <- ifelse( which.max(LengthAge_alyis_beg[a, , y, i, 1]) > length(fish_selex_yafs[y,      , sur_flt, 1]),
        #                    length(fish_selex_yafs[y,    nage  , sur_flt, 1]),
        #                    which.max(LengthAge_alyis_beg[a, , y, i, 1]))
        #     
        #     survey_yf_pred[y,sur_flt] <-  survey_yf_pred[y,sur_flt] +
        #       q*
        #       phi_if_surv[sur_flt,i]*
        #       sum( selMult*surv_selex_yafs[y,LAA, sur_flt, 1] *
        #              N_yais_mid[y, a, i, 1] *
                     # wtatlen_kab[phi_ik2[i], 1] *
        #              which.max(LengthAge_alyis_mid[a, , y, i, 1]) ^
        #              wtatlen_kab[phi_ik2[i], 2],
        #            selMult*surv_selex_yafs[y, LAA, sur_flt, 2] *
        #              N_yais_mid[y, a, i, 2] *
        #              wtatlen_kab[phi_ik2[i], 1] *
        #              which.max(LengthAge_alyis_mid[a, , y, i, 2]) ^
        #              wtatlen_kab[phi_ik2[i], 2] )
        #     
        #     Nsamp_acomp_yf[y,sur_flt] <-  Nsamp_acomp_yf[sur_flt]  +
        #       phi_if_surv[sur_flt,i]*
        #     sum(selMult*surv_selex_yafs[y, LAA, sur_flt, 1]*N_yais_mid[y,a,i,1],
        #         selMult*surv_selex_yafs[y, LAA, sur_flt, 2]*N_yais_mid[y,a,i,2])
        # } ## end LEN selex
 
    
    ## survey age comps w error
    
    ## age comps in catches
    
    
    # Catch_yaf_est(y,a,fish_flt) = (Freal_yf(a)/(Z(a)))*(1-exp(-Z(a)))*
    #   phi_if_fish(fish_flt, i)* N_yai_beg(y,a,i)*wage_catch(a,y); ## do this by fleet with phi
    # CatchN_yaf(y,a,fish_flt) = (Freal_yf(a)/(Z(a)))*(1-exp(-Z(a)))* phi_if_fish(fish_flt, i)* N_yai_beg(y,a,i);## Calculate the catch in kg
    # Catch_yf_est(y,fish_flt)= Catch_yf_est[y,a,fish_flt] + Catch_yaf_est(y,a,fish_flt); ## sum over the current catch at age
    # CatchN(y,fish_flt) <- CatchN_yaf[y,a,] CatchN_yaf(y,a,fish_flt);
    
    
    
    
  } ## END YEARS
  # df.out   <- list(N.save = Nsave,
  #                  SSB = SSB,
  #                  N.save.age = N.save.age,
  #                  R.save = R.save,
  #                  V.save = V.save,
  #                  SSB.all = SSB.all,
  #                  Catch.save.age = Catch.save.age,
  #                  CatchN.save.age = CatchN.save.age,
  #                  Catch = Catch,
  #                  Catch.age = Catch.age,
  #                  Catch.quota = Catch.quota,
  #                  Catch.quota.N = Catch.quota.N,
  #                  Fout = Fout.save,
  #                  age_comps_OM = age_comps_OM,
  #                  age_catch = age_comps_catch,
  #                  SSB_0 = SSB_0,
  #                  N0 = N0,
  #                  SSB.weight = SSB.weight,
  #                  survey.true = survey.true,
  #                  Z = Z.save,
  #                  survey = as.numeric(survey),
  #                  age_comps_surv = age_comps_surv,
  #                  age_comps_country = age_comps_surv_space,
  #                  age_comps_catch_space = age_comps_catch_space,
  #                  Fseason = Fseason.save,
  #                  Fsel = Fsel.save,
  #                  Ninit = Ninit,
  #                  SSB0 = SSB_0)
  # 
  # return(df.out)
  
  
}  ## END FUNC

# Ninit <- rep(NA,nage)
# Ninit_dev <- (df$parms$initN)
# Ninit[1] <- R0
# Ninit[2:(nage-1)] <-R0 * exp(-mat_age[2:(nage-1)]*age[2:(nage-1)])*exp(-0.5*SDR^2*0+Ninit_dev[1:(nage-2)])
# Ninit[nage] <- R0*exp(-(mat_age[nage-1]*age[nage-1]))/(1-exp(-mat_age[nage]))*exp(-0.5*SDR^2*0+Ninit_dev[nage-1])# Plus group (ignore recruitment dev's in first year )
# 
# Create containers to save the data
# SSB_init <- NA
# 
# for(i in 1:nspace){
#   SSB_init[i] <- sum(df$Matsel*Ninit*move.init[i], na.rm =T)*0.5
# }
# 




## true survey biomass


## catches

#Ninit[1] <- sum((4*h*R_0*SSB_init/(SSB_0*(1-h)+ SSB_init*(5*h-1)))*exp(-0.5*1*SDR^2+df$parms$Rin[1]), na.rm = T)
# year_1 <- c(year,max(year)+1)
# 
# SSB <- matrix(NA,nyear, nspace, 
#               dimnames = list(year = df$years,
#                               space = 1:nspace))
# SSB.all <- array(NA, dim = c(nyear, nspace,nseason),
#                  dimnames = list(year = year, space = 1:nspace, season = 1:nseason))
# SSB.weight <- matrix(NA,nyear, nspace,
#                      dimnames = list(year = year, space = 1:nspace))
# Biomass.save <- matrix(NA,nyear, nspace, 
#                        dimnames= list(year = year, space = 1:nspace))
# Catch <- matrix(NA,nyear, dimnames = list(year = year))
# Catch.age <- matrix(NA,nage,nyear, dimnames = list(age = age, year = year))
# CatchN <- matrix(NA,nyear, dimnames = list(year = year))
# CatchN.age <- matrix(NA,nage,nyear, dimnames = list(age =age, year = year))
# 
# 
# R.save <- matrix(NA,nyear, nspace, dimnames = list(year = year, space = 1:nspace))
# Fsel.save <- array(NA,dim = c(nage,nyear,nspace), dimnames = list(age = age, year = year, space = 1:nspace))
# Fseason.save <- array(NA,dim = c(nage, nyear, nspace,nseason), dimnames = list(age = age, year = year, space = 1:nspace,
#                                                                                season = 1:nseason))
# Fout.save <- array(NA, dim = c(nyear,nseason,nspace), 
#                    dimnames = list(year = year, season = 1:nseason, space = 1:nspace))
# 
# N.save.age <- array(NA,dim = c(nage,nyear+1, nspace, nseason), 
#                     dimnames = list(age = age, year = year_1, space = 1:nspace, season = 1:nseason))
# N.save.age.mid <- array(NA,dim = c(nage,nyear+1, nspace, nseason), 
#                         dimnames = list(age = age, year = year_1, space = 1:nspace, season = 1:nseason))
# R.save <- matrix(NA, nyear, nspace)
# V.save <- array(NA,dim = c(nyear, nspace, nseason), dimnames = list(
#   year = year, space = 1:nspace, season = 1:nseason))
# 
# Catch.save.age <- array(NA,dim = c(nage,nyear, nspace, nseason), 
#                         dimnames = list(age = age, year = year, space = 1:nspace, season = 1:nseason))
# CatchN.save.age <- array(NA,dim = c(nage,nyear, nspace, nseason), 
#                          dimnames = list(age = age, year = year, space = 1:nspace, season =1:nseason))
# Catch.quota <- array(NA, dim = c(nyear, nspace, nseason), 
#                      dimnames = list(year = year, space = 1:nspace, season =1:nseason))
# Catch.quota.N <- array(0, dim = c(nyear, nspace, nseason), dimnames = list(year = year, space = 1:nspace,
#                                                                            season = 1:nseason))
# 
# survey <- array(NA,dim = c(nyear), dimnames = list(year = year))
# survey.true <- array(NA, dim = c(nspace, nyear), dimnames = list(space = 1:nspace, year = year))
# surv.tot <- matrix(NA, nyear,nspace, dimnames = list(year = year, space = 1:nspace))
# 
# age_comps_surv <- array(NA, dim = c(df$age_maxage,nyear), dimnames = list(age = 1:df$age_maxage,
#                                                                           year = year)) # 
# age_comps_surv_space <- array(NA, dim = c(df$age_maxage,nyear,nspace), dimnames = list(
#   age = 1:df$age_maxage, year = year))
# 
# N.survey <- matrix(NA,df$age_maxage ,nyear, dimnames = list(age = 1:df$age_maxage,
#                                                             year= year))
# 
# age_comps_catch <- array(NA, dim = c(df$age_maxage,nyear), dimnames = list(age = 1:df$age_maxage,
#                                                                            year = year))
# age_comps_catch_space <- array(NA, dim = c(df$age_maxage,nyear,nspace), dimnames = list(
#   age = 1:df$age_maxage, year = year, space = 1:nspace))
# 
# age_comps_OM <- array(NA, dim = c(nage,nyear, nspace,nseason), 
#                       dimnames = list(age = age, year= year, space = 1:nspace, season = 1:nseason))
# 
# Z.save <- array(NA, dim = c(df$nage, nyear,nspace,nseason), dimnames = list(age= age, year = year, space = 1:nspace,
#                                                                             season = 1:nseason))
# 
# Z.save[,1,1,1] <- M
# Catch.age[,1] <- 0 # Assumed no fishing before data started 
# Catch[1]<- 0
# 
# CatchN[1] <- 0
# CatchN.age[,1] <- 0
# 
# survey[1] <- 1 # Surveys start later
# 
# for (space in 1:nspace) {
#   survey.true[space, 1] <-
#     sum(N.save.age[, 1, space, df$surveyseason] * surv.sel * q * df$wage_survey[, 1])
# }
# 
# idx.save <- seq(1,tEnd, by = nseason)
# 
# # Distribute over space 
# Ninit <- rep(NA,nage)
# names(Ninit) <- age
# Ninit_dev <- (df$parms$initN)
# 
# Ninit[2:(nage-1)] <-R0 * exp(-Mage[2:(nage-1)])*exp(-0.5*SDR^2*0+Ninit_dev[1:(nage-2)])
# Ninit[nage] <- R0*exp(-(mat_age[nage]*age[nage]))/(1-exp(-mat_age[nage]))*exp(-0.5*SDR^2*0+Ninit_dev[nage-1])# Plus group (ignore recruitment dev's in first year )
# 
#p.save <-matrix(NA,tEnd)


# for (space in 1:nspace){
#   # if (season == 1){
#   N.save.age[,1,space,1] <- Ninit*move.init[space] # Just to initialize 
#   N.save.age.mid[,1,space,1] <- N.save.age[,1,space,1]*exp(-0.5*(M/nseason))
#   # }else{
#   #   N.save.age[,1,space,season] <- N.save.age[,1,space,season-1]*exp(-M/nseason)
#   #   N.save.age.mid[,1,space,season] <- N.save.age[,1,space,season]*exp(-0.5*(M/nseason))
#   # }
#   # }
# }



# 
# 
# Fspace <- c(0.155612,0.7388) # Contribution of Total catch (add to one)    #Z <- (Fyear+Myear)
# Fnseason <- df$Fnseason
# pope.mul <- nseason/1*0.5
# pope.mul <- 0.50
# 
# if(nseason == 1){
#   
#   Fnseason <- matrix(rep(1, df$nspace))
#   
# }
# 
# 
# for (yr in 1:nyear){ # Loop over years add one year for initial distribution
#   
#   #if(year[yr] < year[df$selYear] | year[yr] > 2017){
#   # }else{
#   #   psel <- df$parms$psel_fish+df$parms$PSEL[,yr-df$selYear+1]
#   # }
#   if(year[yr] < 2019){
#     w_catch <- df$wage_catch[,yr]
#     w_surv <- df$wage_survey[,yr]
#     w_mid <- df$wage_mid[,yr]
#     w_ssb <- df$wage_ssb[,yr]
#   }else{
#     w_catch <- df$wage_catch[,1]
#     w_surv <- df$wage_survey[,1]
#     w_mid <- df$wage_mid[,1]
#     w_ssb <- df$wage_ssb[,1]
#   }
#   
#   
#   Ry <- df$parms$Rin[yr]
#   
#   
#   # Fyear <- F0[yr]*Fsel
#   Myear <- M # Natural mortality 
#   
#   ## add these to load data seasons 
#   # Fnseason <- matrix(1, nseason)
#   # Fnseason <- Fnseason/sum(Fnseason)
#   # Fnseason <- c(0,0.5,0.5,0)
#   
#   
#   if(df$move == FALSE){
#     Fspace <- 1 # All catches in the south
#   }
#   
#   Mseason <- Myear/nseason # M is distributed throughout the year
#   
#   
#   # fix Ssb and recruitment in all areas 
#   for(space in 1:nspace){
#     SSB.weight[yr,space] <- sum(N.save.age[,yr,space,1]*as.numeric(w_ssb), na.rm = TRUE)*0.5
#     SSB[yr,space] <- SSB.weight[yr,space] #sum(N.save.age[,yr,space,1]*Mat.sel, na.rm = TRUE)
#     
#     SSB.all[1,space,1]<- sum(N.save.age[,1,space,1]*Mat.sel, na.rm = TRUE)*0.5
#     
#     # Recruitment only in season 1  
#     R <- (4*h*R_0[space]*SSB[yr,space]/
#             (SSB_0[space]*(1-h)+ SSB[yr,space]*(5*h-1)))*exp(-0.5*df$b[yr]*SDR^2+Ry)#*recruitmat[space]
#     
#     N.save.age[1,yr,space,1] <- R
#     R.save[yr,space] <- R
#   }
#   
#   
#   for (season in 1:nseason){
#     for (space in 1:nspace){
#       
#       # Get the selectivity of the season and area 
#       psel <- df$psel[space,] 
#       
#       
#       if(df$flag_sel[yr] == 1){
#         pseltmp <- psel+df$parms$PSEL[,yr-df$selidx+1]*df$sigma_psel
#       }else{
#         pseltmp <- psel
#       }
#       
#       
#       if(year[yr] >2018){
#         
#         # if(df$selectivity_change == 0){
#         #   if(space == 1){
#         #     pseltmp <- c(1,1,1,1,1)
#         #   }else{
#         #   pseltmp <- psel
#         #   }
#         # }
#         
#         if(df$selectivity_change ==1){
#           if(space == 1){
#             #pseltmp <- psel
#             pseltmp <- c(1,1,1,1,1)
#           }else{
#             pseltmp <- c(0.05,0.05,0,0,0)
#           }
#         }
#         
#         if(df$selectivity_change ==2){
#           pseltmp <- df$psel[2,]+df$parms$PSEL[,ncol(df$parms$PSEL)]*df$sigma_psel}
#         
#       }
#       
#       #p.save[yr] <- sum(pseltmp)
#       # 
#       Fsel <- getSelec(age,pseltmp,df$Smin,df$Smax) # Constant over space right now 
#       rm(pseltmp)
#       
#       Fsel.save[,yr,space] <- Fsel
#       
#       if(nspace > 1){
#         if(df$years[yr]<= 2018){
#           Catch_space <- df$Catch.country[yr,space]
#         }else{
#           Catch_space <- df$Catch[yr]*Fspace[space]  
#         }
#       }else{
#         Catch_space <- df$Catch[yr]
#       }
#       
#       
#       E.temp <- Catch_space*Fnseason[space, season]#*Fspace[space] # Catch distribution in the year
#       B.tmp <-  sum(N.save.age[,yr,space,season]*exp(-Mseason*pope.mul)*w_catch*Fsel) # Get biomass from previous year
#       N.tmp <- N.save.age[,yr,space,season]#
#       V.save[yr,space,season] <- B.tmp
#       Catch.quota[yr,space,season] <- E.temp
#       
#       if(E.temp/B.tmp >= .9){
#         if(df$years[yr] < 2018){
#           stop(paste('Catch exceeds available biomass in year:',year,' and season', season, 'area', space)) # Stop if in the past 
#         }
#         #print(paste('Catch exceeds available biomass in year:',year,' and season', season, 'area', space))
#         E.temp <- 0.75*B.tmp
#         Catch.quota.N[yr,space,season] <- 1
#         #if(df$years[yr] > 2026){
#         #stop('danger')
#         #  }
#         
#       }
#       
#       Fout <- getF(E.temp,B.tmp,Mseason = Mseason, Fsel = Fsel, N.tmp = N.tmp, w_catch = w_catch, 
#                    method = 'Hybrid')
#       
#       Fout <- Fout
#       #Fout <- df$parms$F0[yr]
#       
#       if(E.temp>0){
#         Fseason <- Fout*Fsel
#         Fnew <- Fout  
#         Z <- Fnew*Fsel+Mseason
#         Fseason <- Fnew*Fsel
#       }else{
#         Fseason <- 0
#       }
#       
#       Fout.save[yr,season,space] <- Fout # terminal fishing mortality 
#       
#       Fseason.save[,yr,space,season] <- Fseason
#       
#       Z <- Mseason+Fseason
#       Z.save[,yr,space,season]<- Z
#       
#       
#       ## MK retool this so it auto-detects how many other spaces to sum across for inmigration
#       if(((space-1) == 0)){
#         spaceidx <- 2
#       }
#       if(space == nspace){
#         spaceidx <- nspace-1
#       }
#       if(space > 1 & space < nspace){
#         spaceidx <- c(space-1,space+1)
#       }
#       
#       if(df$move == FALSE){
#         spaceidx <- 1
#       }
#       
#       if(season <nseason){ ## doesnt matter if nseas = 1
#         
#         N.save.age[,yr,space,season+1] <- N.save.age[,yr,space,season]*exp(-Z)-
#           N.save.age[, yr,space,season]*exp(-Z)*(movemat[space,,season,yr])+ # Remove the ones that leave
#           N.save.age[, yr,spaceidx,season]*exp(-Z)*(movemat[spaceidx,,season,yr])# add the ones come to the surrounding areas
#         
#         age_comps_Omat_age[,yr,space,season] <- N.save.age[, yr,space,season]/sum(N.save.age[, yr,space,season])
#         
#         SSB.all[yr,space,season]<- sum(N.save.age[,yr,space,season]*Mat.sel, na.rm = T)
#         Catch.save.age[, yr,space, season] <- (Fseason/(Z))*(1-exp(-(Z)))*N.save.age[,yr,space,season]*w_catch
#         CatchN.save.age[, yr,space, season] <- (Fseason/(Z))*(1-exp(-(Z)))*N.save.age[,yr,space,season]
#         
#         
#       }else{
#         ## MK this needs to be updated to sum across >1 other space otherwise non conformable arrays!
#         N.save.age[2:(nage-1),yr+1,space,1] <- N.save.age[1:(nage-2),yr,space,season]*exp(-Z[1:(nage-2)])-
#           N.save.age[1:(nage-2), yr,space,season]*exp(-Z[1:(nage-2)])*(movemat[space,1:(nage-2),season,yr])+ # Remove the ones that leave
#           N.save.age[1:(nage-2), yr,spaceidx,season]*exp(-Z[1:(nage-2)])*(movemat[spaceidx,1:(nage-2),season,yr])
#         # add the ones come to the surrounding areas
#         
#         # Plus group 
#         Nsurvive.plus <- (N.save.age[nage-1, yr,space, nseason]*exp(-Z[nage-1])+
#                             N.save.age[nage, yr,space, nseason]*exp(-Z[nage]))
#         
#         Nout.plus <- Nsurvive.plus*(movemat[space,nage, season,yr]) # Leaving
#         
#         
#         Nin.plus <- (N.save.age[nage-1, yr,spaceidx,nseason]*exp(-Z[nage-1])+
#                        N.save.age[nage, yr,spaceidx,nseason]*exp(-Z[nage]))*
#           (movemat[spaceidx,nage, season,yr]) # Incoming
#         
#         N.save.age[nage,yr+1,space,1] <- Nsurvive.plus- Nout.plus + Nin.plus
#         
#         
#         age_comps_Omat_age[,yr,space,season] <- N.save.age[, yr,space,season]/sum(N.save.age[, yr,space,season])
#         
#         SSB.all[yr,space,season]<- sum(N.save.age[,yr,space,season]*Mat.sel, na.rm = T)
#         Catch.save.age[, yr,space, season] <- (Fseason/(Z))*(1-exp(-(Z)))*N.save.age[,yr,space,season]*w_catch
#         CatchN.save.age[, yr,space, season] <- (Fseason/(Z))*(1-exp(-(Z)))*N.save.age[,yr,space,season]
#         
#         
#       }
#       
#       
#       if(is.na(SSB[yr,space])){
#         stop('SSB is NA')
#       }
#     } 
#     
#     if(Catch.quota[yr,space,season]>0){
#       if((sum(Catch.save.age[, yr,space, season])/Catch.quota[yr,space,season]) > 1.05){
#         stop('F estimation overshoots more than 10%')
#       }
#     }
#     
#   } # End of season loop
#   
#   
#   
#   #Catch.age[,idx]  <- (Fyear/(Fyear+Myear))*(1-exp(-(Fyear+Myear)))*rowSums(N.save.age[,idx,,1])*w_catch # Calculate the catch in kg 
#   
#   if(nseason>1){
#     Catch.age[,yr] <- apply(Catch.save.age[,yr,,],MARGIN = 1,FUN = sum)
#     Catch[yr] <- sum(Catch.save.age[,yr,,])  
#     
#     CatchN.age[,yr] <- apply(CatchN.save.age[,yr,,],MARGIN = 1,FUN = sum)
#     CatchN[yr] <- sum(CatchN.save.age[,yr,,])  
#   }else{
#     
#     if(nspace == 1){
#       Catch.age[,yr] <- Catch.save.age[,yr,,]
#       Catch[yr] <- sum(Catch.save.age[,yr,,])
#       
#       CatchN.age[,yr] <- CatchN.save.age[,yr,,]
#       CatchN[yr] <- sum(CatchN.save.age[,yr,,])
#     }else{
#       Catch.age[,yr] <- rowSums(Catch.save.age[,yr,,])
#       Catch[yr] <- sum(Catch.save.age[,yr,,])
#       
#       CatchN.age[,yr] <- rowSums(CatchN.save.age[,yr,,])
#       CatchN[yr] <- sum(CatchN.save.age[,yr,,])
#     }
#   }  
#   
#   if(nseason == 1){
#     Msurveymul <- 0.5
#   }else{
#     Msurveymul <- 0
#   }
#   
#   for (space in 1:nspace){
#     survey.true[space,yr] <- sum(N.save.age[,yr,space,df$surveyseason]*
#                                    exp(-Msurveymul*Z.save[,yr,space,df$surveyseason])*surv.sel*q*w_surv)
#     
#   }
#   
#   
#   #  }
#   # Save the survey 
#   # Survey is conducted in the start of the year
#   # }else{
#   #   Msurveymul <- 0.5
#   # }
#   
#   if(df$move == FALSE){
#     Nsurv <- N.save.age[,yr,,df$surveyseason]*
#       exp(-Msurveymul*Z.save[,yr,space,df$surveyseason])
#   }else{
#     Nsurv <- rowSums(N.save.age[,yr,,df$surveyseason]*
#                        exp(-Msurveymul*Z.save[,yr,space,df$surveyseason]))
#   }
#   
#   if (df$flag_surv_acomp[yr] == 1){
#     
#     
#     if(year[yr] > 2018){
#       err <- rnorm(n = 1,mean = 0, sd = surv.sd)
#       surv <- exp(log(sum(Nsurv*surv.sel*q*w_surv))+err) # If the xtra factor is not included the mean is > 1
#     }else{
#       surv <- sum(Nsurv*surv.sel*q*w_surv)
#     }
#     survey[yr] <- surv
#   }else{
#     survey[yr] <- 1
#   }
#   
#   Ntot.year <- Nsurv
#   
#   surv.tmp <- sum(Ntot.year*surv.sel*q)
#   
#   if(df$flag_surv_acomp[yr] == 1){
#     age_comps_surv[1,yr] <- 0 # No year 1 recorded
#     
#     age_comps_surv[1:(df$age_maxage-1),yr] <-  (Ntot.year[2:(df$age_maxage)]*surv.sel[2:(df$age_maxage)]*q)/surv.tmp
#     age_comps_surv[df$age_maxage,yr] <- sum(Ntot.year[(df$age_maxage+1):nage]*surv.sel[(df$age_maxage+1):nage]*q)/surv.tmp
#   }else{
#     age_comps_surv[,yr] <- NA
#   }
#   
#   for(space in 1:nspace){
#     Ntot.year <- N.save.age[,yr,space,df$surveyseason]
#     surv.tot[yr,space]  <- sum(Ntot.year*surv.sel*q*exp(-Msurveymul*Z.save[,yr,space,df$surveyseason]))
#     
#     age_comps_surv_space[1,yr,space] <- 0 # No year 1 recorded
#     
#     age_comps_surv_space[1:(df$age_maxage-1),yr,space] <-  
#       (Ntot.year[2:(df$age_maxage)]*surv.sel[2:(df$age_maxage)]*q)/surv.tot[yr,space]
#     age_comps_surv_space[df$age_maxage,yr,space] <- 
#       sum(Ntot.year[(df$age_maxage+1):nage]*surv.sel[(df$age_maxage+1):nage]*q)/surv.tot[yr,space]
#     
#     if(nseason>1){
#       Catch.tmp <- rowSums(CatchN.save.age[, yr,space,])
#     }else{
#       Catch.tmp <- CatchN.save.age[, yr,space,]
#     }
#     
#     Catch.tot <- sum(CatchN.save.age[,yr,space,])
#     
#     age_comps_catch_space[1:(df$age_maxage-1),yr,space] <- Catch.tmp[2:(df$age_maxage)]/Catch.tot
#     age_comps_catch_space[df$age_maxage,yr,space] <- sum(Catch.tmp[(df$age_maxage+1):nage])/Catch.tot
#     
#     
#     
#   }
#   
#   #
#   if(df$flag_catch[yr] == 1){
#     age_comps_catch[1:(df$age_maxage-1),yr] <-  CatchN.age[2:(df$age_maxage),yr]/CatchN[yr]
#     age_comps_catch[df$age_maxage,yr] <- sum(CatchN.age[(df$age_maxage+1):nage,yr])/CatchN[yr]
#     
#   }else{
#     age_comps_catch[,yr] <- NA
#   }
#   
#   
# }# End of year loop
# #}
# 
# if(df$move == FALSE){
#   Nsave <- N.save.age[,,,nspace]
#   SSB.save <- SSB
# }else{  
#   Nsave <- apply(N.save.age[,,,1],2,rowSums)
#   SSB.save <- rowSums(SSB)
# }
# 
# # Add names to output 
# year_1 <- c(df$years,max(df$years+1))
# 
# 


## note this has NO estimation/optimization, just the deterministic outputs of the OM given the structure at hand.
# df.out   <- list(N.save = Nsave, 
#                  SSB = SSB, 
#                  N.save.age = N.save.age,
#                  R.save = R.save,
#                  V.save = V.save,
#                  SSB.all = SSB.all,
#                  Catch.save.age = Catch.save.age,
#                  CatchN.save.age = CatchN.save.age,
#                  Catch = Catch, 
#                  Catch.age = Catch.age, 
#                  Catch.quota = Catch.quota,
#                  Catch.quota.N = Catch.quota.N,
#                  Fout = Fout.save,
#                  age_comps_OM = age_comps_OM,
#                  age_catch = age_comps_catch,
#                  SSB_0 = SSB_0, 
#                  N0 = N0,
#                  SSB.weight = SSB.weight,
#                  survey.true = survey.true,
#                  Z = Z.save,
#                  survey = as.numeric(survey),
#                  age_comps_surv = age_comps_surv,
#                  age_comps_country = age_comps_surv_space,
#                  age_comps_catch_space = age_comps_catch_space,
#                  Fseason = Fseason.save,
#                  Fsel = Fsel.save, 
#                  Ninit = Ninit,
#                  SSB0 = SSB_0)
# 
# return(df.out)



