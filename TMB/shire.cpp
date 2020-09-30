// M KAPUR mod N Jacobsen
// 2020 kapurm@uw.edu
#include <TMB.hpp>
#include <iostream>

// TO DO
// correct growth for early ages
// make selex fleet specific
// introduce selectivity estimation
// need error on tau
// everything sex specific
// calculate reference points
// M vs myear -- at age or what?
// estimate SDR (currently parameter)


template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA AND PARAMETERS BY CATEGORY //
  
  // structure //
  DATA_INTEGER(nspace); // number of subreas to track
  DATA_INTEGER(nstocks); // number of stocks (demography)
  DATA_INTEGER(nage); // Plus group
  DATA_VECTOR(age); // ages
  DATA_INTEGER(tEnd); // number of years modeled
  DATA_VECTOR(year); // number of years modeled
  int nyear = year.size();
  
  DATA_INTEGER(nfleets_surv); // number of survey fleets
  DATA_INTEGER(nfleets_fish); //number of fishery fleets
  DATA_INTEGER(nfleets_acomp); // number of age comp fleets
  DATA_INTEGER(nfleets_lcomp); //number of len comp fleets
  DATA_INTEGER(nmgmt_reg); // mgmt regions (3)

  // DATA_ARRAY(phi_if_surv); // turn on/off subareas for survey fleets
  DATA_ARRAY(phi_if_fish); // turn on/off subareas for fishery fleets
  DATA_ARRAY(phi_ki); // 0/1 nesting of subareas i into stocks k (rows)
  DATA_IVECTOR(phi_ik2); // vector stating which subarea (col) belongs to each stock k (value)
  DATA_ARRAY(tau_ki); // downscaling from stocks to sub-areas
  // 
  // biology // 
  DATA_VECTOR(mat_age); // natural mortality at age

  // movement //
  PARAMETER_VECTOR(omega_0ij); // estimated age-0 movment among areas (used upon recruit gen)
  DATA_ARRAY(omega_ais); // eigenvect of movement between subareas for ages > 0
  DATA_ARRAY(X_ijas); // prob trans between subareas at age
  
  // growth //
  DATA_ARRAY(unfished_ALK_F); // for use in SSB0 calcs
  DATA_ARRAY(wtatlen_kab); // aL^b values by stock
  DATA_ARRAY(Linf_yk); // sex, stock, year specific
  DATA_ARRAY(L1_yk); // length at age 4 by stock; linear before this
  DATA_ARRAY(kappa_yk);
  DATA_ARRAY(sigmaG_yk); // perhaps turn to parameter later
  DATA_ARRAY(phi_ij); // matrix of whether i,j are from distinct stocks (0 otherwise)
  DATA_INTEGER(LBins); // maximum length bin in CM
  DATA_ARRAY(mla_yais); // for early or late periods, most likely length at age(for selex)
  
  array<Type> Length_yais_beg(tEnd+1,nage,nspace,2); // placeholder for true lengths-at-age
  array<Type> Length_yais_mid(tEnd+1,nage,nspace,2); // placeholder for true lengths-at-age
  array<Type> LengthAge_alyis_beg(nage,LBins,tEnd+1,nspace,2); // placeholder for true age-length dist
  array<Type> LengthAge_alyis_mid(nage,LBins,tEnd+1,nspace,2); // placeholder for true age-length dist

  // repro //
  PARAMETER_VECTOR(logR_0k); // Recruitment at equil by stock
  PARAMETER_VECTOR(logh_k); // Steepness by stock
  // DATA_SCALAR(logSDR); // Can it be estimated as a fixed effect?
  DATA_ARRAY(mat_ak); // maturity at age for stock
  
  // repro storage
  vector<Type> R(tEnd);
  array<Type>  R_yk(tEnd,nstocks); // stock-level recruitment (bev-holt)
  array<Type>  R_yi(tEnd,nspace); // subarea-level recruitment (downscaled)
  
  // observations //
  PARAMETER_VECTOR(b); // bias adjustment factor
  // DATA_INTEGER(nage); // Last age included in age comps

  // DATA_ARRAY(age_catch); // Age comps in catch -- should be by fleet
  
  // Selectivity
  DATA_IVECTOR(selType_fish); // 0 == AGESEL, 1= LENSEL
  DATA_IVECTOR(selShape_fish); 
  DATA_IVECTOR(selType_surv); // 0 == AGESEL, 1 = LENSEL
  DATA_IVECTOR(selShape_surv); 
  PARAMETER_ARRAY(log_fsh_slx_pars);       // Fishery selectivity (selShape controls parameterization)
  PARAMETER_ARRAY(log_srv_slx_pars);       // Survey selectivity (selShape controls parameterization)
  
  // Switch for selectivity type: 0 = a50, a95 logistic; 1 = a50, slope logistic
  // DATA_INTEGER(slx_type)
    // Predicted selectivity
  array<Type> fsh_slx_yafs(nyear, nage, nfleets_fish,2);           // Fishery selectivity-at-age by sex (on natural scale)
  array<Type> srv_slx_yafs(nyear, nage, nfleets_surv,2);           // Survey selectivity-at-age by sex(on natural scale)
  // Time varying parameter blocks (indexed as h) - each vector contains the terminal years of
  // each time block. Used for both selectivity and catchability
  DATA_IVECTOR(fsh_blks);        // fishery  
  DATA_IVECTOR(srv_blks);       // survey
    
  // Survey Biomass
  // DATA_VECTOR(survey_err);
  // DATA_VECTOR(survey); // Acoustic survey - vector of obs x year 
  // DATA_ARRAY(survey_bio_f_obs); // year x fleets relbio
  // PARAMETER(logSDsurv); // Survey uncertainty
  
  // Survey Comps
  // DATA_VECTOR(survey); // Age comp sample size
  // DATA_ARRAY(age_survey); // Age compositions, age x year
  // DATA_ARRAY(age_error); // Age compositions, age x fleet (row 1 = true age, row 2= biased age, row 3 = sd)
  
  // DATA_ARRAY(survey_acomp_f_obs); // Observed survey compositions,  year x age x nfleets_acomp {right now just survnfleets}
  // array<Type> acomp_yaf_temp(tEnd, age_maxage, nfleets_acomp); // placeholder for aging error calcs
  // array<Type> survey_acomp_f_est(tEnd, age_maxage, nfleets_acomp); //when error multiplied by Nage
  // vector<Type> Nsamp_acomp_f(nfleets_acomp); // placeholder for number sampled by comp survey (pre dirichlet weighting)
  
  // Catches
  DATA_ARRAY(catch_yf_obs); // obs catch by year and fleet
  // DATA_SCALAR(logSDcatch); // Error on catch, should be by fleet
  // PARAMETER(logphi_catch);
  
  // Catch Comps
  // array<Type> catch_acomp_f_est(tEnd, age_maxage, nfleets_fish); // estimated catch comps; uses derived quants
  
  // Catch storage
  array<Type> catch_yaf_pred(tEnd, nage, nfleets_fish);  // estimated catches at age by fleet
  array<Type> catch_yf_pred(tEnd,nfleets_fish); 
  array<Type> catch_yfi_pred(tEnd,nfleets_fish,nspace); 
  array<Type> catch_yaif_pred(tEnd,nage,nspace,nfleets_fish);  
  
  array<Type> CatchN_yaf(tEnd,nage,nfleets_fish);

  // array<Type> CatchN(tEnd,nfleets_fish);
  
  // F tuning storage
  int niter = 50;
  array<Type> F1_yf(tEnd,nfleets_fish+1, niter+1); // intermediate f guess storage
  array<Type> F2_yf(tEnd,nfleets_fish+1, niter+1); // intermediate f guess storage 
  array<Type> Freal_yf(tEnd,nfleets_fish); // final tuned fleet and yr specific F
  array<Type> Zreal_ya(tEnd,nage); // temp tuned fleet Z by y and age
  array<Type> Zreal_yai(tEnd,nage,nspace); // temp tuned fleet Z by y and age and area
  array<Type> F_area_yfi(tEnd,nfleets_fish,nspace); // temp tuned fleet Z by y and age
  array<Type> F_ym(tEnd,nmgmt_reg); 

  // array<Type> Ftuned_yf
  // array<Type> Ftuned_ym(tEnd,nfleets_fish); // summation of Fs nested in mgmt regions
  // array<Type> Ztuned_yf(tEnd, nage, nspace); // final tuned subarea and yr specific Z
  // vector<Type> Fyear(tEnd);
  // vector<Type> Freal(nage);
  // vector<Type> Z(nage);

  // Priors
  // DATA_SCALAR(Bprior);
  // DATA_SCALAR(Aprior);
  
  
  // biology storage
  array<Type> Ninit_ais(nage,nspace,2); // initial numbers at age in subarea, just once
  array<Type> N_0ais(nage, nspace,2); // numbers in year 0 at age in subarea
  vector<Type> SSB_0k(nstocks); // virgin spawnbio by stock
  vector<Type> SSB_0i(nspace); // virgin spawnbio by subarea
  array<Type> N_yais_beg( tEnd+1, nage, nspace,2); N_yais_beg.setZero();
  array<Type> N_yais_mid( tEnd+1, nage, nspace,2); N_yais_mid.setZero();
  array<Type> N_yais_end( tEnd+1, nage, nspace,2); N_yais_end.setZero();
  array<Type> SSB_yk(tEnd,nstocks);
  array<Type> SSB_yi(tEnd,nspace);
  // array<Type> surv_yf_pred(tEnd,nfleets_surv); // this is actually predicted
  // array<Type> Zsave(nage,tEnd);
  
  PARAMETER(logSDR);
  // PARAMETER(logphi_survey);
  
  // Transform out of log space
  // Type SDsurv = exp(logSDsurv);
  // Type SDcatch = exp(logSDcatch);
  Type SDR = exp(logSDR);
  vector<Type> R_0k = exp(logR_0k);
  vector<Type> h_k = exp(logh_k);
  // 
  // Type Minit = exp(logMinit);
  // Type q = exp(logQ);
  // Type phi_survey = exp(logphi_survey);
  // Type phi_catch = exp(logphi_catch);
  
  //  Minor calculations
  // vector<Type> M = Minit*Msel; // Natural mortality
  // vector<Type> Myear = M*Msel; // Natural mortality (if we want to change it later)
  // vector<Type> Zzero = M; // total mortality without fishing
  // vector<Type> logF(tEnd);
  array<Type> tildeR_yk(tEnd,nstocks); // recdevs
  vector<Type> tildeR_initk(nstocks); // recdevs for init
 
 
 // from SEAK
 
 // Fishery selectivity
 // Number of parameters in the chosen selectivity type: 
 int npar_slx = log_fsh_slx_pars.dim(1); // dim = array dimensions; 1 = # columns in array = # params in slx_type
 // // Preliminary calcs to bring parameters out of log space
 array<Type> fsh_slx_pars(log_fsh_slx_pars.dim);
 fsh_slx_pars.setZero();
 for (int k = 0; k < 2; k++) {
   for (int h = 0; h < fsh_blks.size(); h++) {
     for (int n = 0; n < npar_slx; n++) {
       fsh_slx_pars(h,n,k) = exp(log_fsh_slx_pars(h,n,k));
     }
   }
 }
 // // Notes on the following syntax: the do while allows you to estimate parameters within a y block. It
 // // "does" the looping over year and age "while" within the y block, then
 // // iterates to the next block. Year is not in a for loop because it is
 // // iterated by the do statement.
 // 
 // // The switch for slx_shape allows you to change parameterization SHAPE. This could
 // // easily be expanded to accomodate any selectivity type (the fsh_slx_pars
 // // allows for a flexible number of parameters and y blocks)
 // 
 int i = 0;
 for(int y = 0; y < fsh_blks.size(); y++){
   do{
     for(int fish_flt =0;fish_flt<(nfleets_fish);fish_flt++){
       switch (selType_fish(fish_flt)) { // 0 is age, 1 is leng
       case 0: // enter age based sel
         for (int s = 0; s < 2; s++) {
           for (int a= 0; a < nage; a++) {
             // Selectivity switch (case 0 or 1 references the value of slx_type)
             switch (selShape_fish(fish_flt)) {
             case 0: // Logistic with a50 and a95, where fsh_slx_pars(y,0,s) = a50 and fsh_slx_pars(y,1,s) = a95
               fsh_slx_yafs(i,a,fish_flt,s) = Type(1.0) / ( Type(1.0) + exp(-log(Type(19)) * (a - fsh_slx_pars(y,0,s)) / (fsh_slx_pars(y,1,s) - fsh_slx_pars(y,0,s))) );
               break;
             case 1: // Logistic with a50 and slope, where fsh_slx_pars(y,0,s) = a50 and fsh_slx_pars(y,1,s) = slope.
               //  *This is the preferred logistic parameterization b/c it reduces parameter correlation*
               fsh_slx_yafs(i,a,fish_flt,s)  = Type(1.0) / ( Type(1.0) + exp( Type(-1.0) * fsh_slx_pars(y,1,s) * (a - fsh_slx_pars(y,0,s)) ) );
               break;
             case 2: // Dome Normal with alpha (mean) and beta (sd)
               fsh_slx_yafs(i,a,fish_flt,s)  = exp(-(0.5 * (a - fsh_slx_pars(y,2,s))/pow(fsh_slx_pars(y,3,s),2)));
               break; 
             case 3: // Dome Gamma with alpha (mean) and beta (sd)
               fsh_slx_yafs(i,a,fish_flt,s)  =  pow(a, (fsh_slx_pars(y,2,s) - 1)) * exp(-a/fsh_slx_pars(y,3,s));
               break;
             } // end switch selShape
           } // end ages
         } // end sex
         break;
       case 1: // enter length based sel
         for (int s = 0; s < 2; s++) {
           for (int l= 0; l < LBins; l++) {
             // switch (selShape_fish(fish_flt)) {
             fsh_slx_yafs(i,l,fish_flt,s) = 1.0; //placeholder
             break;
             // } // end switch selShape
           } // end length
         } // end sex
         break;
       } // end selType
     } // end fish_flt
     i++;
   } while (i <= fsh_blks(y));
 } // end y blocks
 // 
 // // std::cout << fsh_slx(1,1,1) << "\n Fishery selectivity \n";
 // 
 // // Survey selectivity - see notes on syntax in fishery selectivity section
 // 
 // // Preliminary calcs to bring parameters out of log space
 // array<Type> srv_slx_pars(log_srv_slx_pars.dim);
 // srv_slx_pars.setZero();
 // 
 // for (int k = 0; k < 2; k++) {
 //   for (int h = 0; h < srv_blks.size(); h++) {
 //     for (int n = 0; n < npar_slx; n++) { 
 //       srv_slx_pars(h,n,k) = exp(log_srv_slx_pars(h,n,k));
 //     }
 //   }
 // }
 // 
 // i = 0;     // re-set i to 0 (do not redeclare)
 // 
 // for(int h = 0; h < srv_blks.size(); h++){
 //   do{
 //     for (int k = 0; k < 2; k++) {
 //       for (int j = 0; j < nage; j++) {
 //         
 //         // Selectivity switch (case 0 or 1 references the value of slx_type)
 //         switch (slx_type) {
 //         
 //         case 0: // Logistic with a50 and a95, where srv_slx_pars(h,0,k) = a50 and srv_slx_pars(h,1,k) = a95
 //           srv_slx(i,j,k) = Type(1.0) / ( Type(1.0) + exp(-log(Type(19)) * (j - srv_slx_pars(h,0,k)) / (srv_slx_pars(h,1,k) - fsh_slx_pars(h,0,k))) );
 //           break;
 //           
 //         case 1: // Logistic with a50 and slope, where srv_slx_pars(h,0,k) = a50 and srv_slx_pars(h,1,k) = slope.
 //           //  *This is the preferred logistic parameterization b/c it reduces parameter correlation*
 //           srv_slx(i,j,k) = Type(1.0) / ( Type(1.0) + exp( Type(-1.0) * srv_slx_pars(h,1,k) * (j - srv_slx_pars(h,0,k)) ) );
 //           
 //           break;
 //         }
 //       }
 //     }
 //     i++;
 //   } while (i <= srv_blks(h));
 // }
 // selectivity
 
 
  //  END DATA & PARS, BEGIN MODEL //
  
  // recdevs placeholder
  for(int k=0;k<(nstocks);k++){
    for(int y=0;y<(tEnd-1);y++){
      tildeR_yk(y,k) =0;
    }
    tildeR_yk(tEnd-1,k) =0;
  }
  for(int k=0;k<(nstocks);k++){
      tildeR_initk(k) =0;
  }

  // Equilibrium Unfished numbers-at-age, subarea (outside of y loop) 
  // identical to Ninit except no recdevs
  N_0ais.setZero();   
  for(int s=0;s<2;s++){
    for(int k=0;k<(nstocks);k++){
      for(int i=0;i<(nspace);i++){ 
        for(int a=0;a<(nage-1);a++){
          N_0ais(a,i,s) = 0.5*omega_ais(a,i,s)*R_0k(k)*tau_ki(k,i)*exp(-(mat_age(a)*age(a))); // compound multiply duh
        }  // note the A+ group will be in slot A-1
        N_0ais(nage-1,i,s) = omega_ais(nage-1,i,s)* N_0ais(nage-2,i,s)*exp(-sum(mat_age)) 
          /(Type(1.0)-exp(-mat_age(nage-1)));
      } // end subareas
    } // end stocks
  } // end sex
  
  // // Equilibrium Unfished SSB, stock (outside of y loop)
  for(int i=0;i<(nspace);i++){
    for(int a=0;a<nage;a++){ // Loop over ages
      SSB_0i(i) += mat_age(a)*
        N_0ais(a,i,0)*
        wtatlen_kab(phi_ik2(i),1)*
        pow(unfished_ALK_F(a,i),wtatlen_kab(phi_ik2(i),2))*
        mat_ak(a,phi_ik2(i));
      for(int k=0;k<(nstocks);k++){
        SSB_0k(k) += phi_ki(k,i)*SSB_0i(i);
      } // end stocks
    } // end ages
  } // end space

  // The first year of the simulation is initialized with the following age distribution
  Ninit_ais.setZero(); 
  for(int y=0;y<(10*nage);y++){
    for(int s=0;s<2;s++){
      for(int i=0;i<(nspace);i++){
        for(int a=0;a<(nage-1);a++){
          Ninit_ais(a,i,s) += 0.5* omega_ais(a,i,s) * 
            tau_ki(phi_ik2(i),i) * 
            R_0k(phi_ik2(i))* exp(-(mat_age(a)*age(a))) *
            exp(-0.5*SDR*SDR+tildeR_initk(phi_ik2(i)));
        } // end ages
        Ninit_ais(nage-1,i) += (omega_ais(nage-1,i,s) * Ninit_ais(nage-2,i,s) *
          exp(-mat_age(nage-1)*age(nage-1)))/(Type(1.0)-exp(-sum(mat_age))* exp(-0.5*SDR*SDR+tildeR_initk(phi_ik2(i))));
      } // end space
    } // end sex
  } // end yinit

  for(int y=0;y<(tEnd);y++){ // Start y loop
    // model year zero, use last year of Ninit_ai, and equil movement (omega) and downscaling (tau)
    // note we are assuming unfished here as the exponent is M only
    // add length at age initial here so that y loop can be removed
    if (y == 0){
      for(int i=0;i<(nspace);i++){
        for(int s=0;s<2;s++){
          Length_yais_beg(0,0,i,s) <- L1_yk(y,phi_ik2(i),s);
          N_yais_beg(0,0,i,s) <- Ninit_ais(0,i,s);
          N_yais_mid(0,0,i,s) <- N_yais_beg(0,0,i,s)*exp(-mat_age(0)/2);
          for(int a=1;a<(nage-1);a++){ // we will fill recruits (a0) later
            Type pLeave = 0.0; Type NCome = 0.0; // reset for new age
            for(int j=0;j<(nspace);j++){ 
              if(i != j){
                pLeave += X_ijas(i,j,a,s); // will do 1-this for proportion which stay
                NCome += X_ijas(j,i,a,s)*Ninit_ais(a,j,s); // actual numbers incoming
              } // end i != j
            } // end subareas j
            Length_yais_beg(y,a,i,s) = Linf_yk(0,phi_ik2(i),s)+(L1_yk(0,phi_ik2(i),s)-Linf_yk(0,phi_ik2(i),s))*
              exp(-kappa_yk(0,phi_ik2(i),s)*a);
            Length_yais_mid(y,a,i,s) =  Linf_yk(0,phi_ik2(i),s)+(L1_yk(0,phi_ik2(i),s)-Linf_yk(0,phi_ik2(i),s))*
              exp(-0.5*kappa_yk(0,phi_ik2(i),s)*a);
            N_yais_beg(y,a,i,s) = ((1-pLeave)*Ninit_ais(a,i,s) + NCome)*exp(-mat_age[a]/2);
          } // end ages
          Type pLeave = 0.0; Type NCome = 0.0; // reset for plusgroup age
          for(int j=0;j<(nspace);j++){ 
            if(i != j){
              pLeave += X_ijas(i,j,nage-1,s);
              NCome += X_ijas(j,i,nage-1,s)*(Ninit_ais(0,nage-1,j,s) + Ninit_ais(0,nage-2,j,s)) ; // if M becomes spatial use M_aj here
            } // end i != j
          } // end subareas j
          N_yais_beg(y,nage-1,i) =  ((1-pLeave)*(Ninit_ais(0,nage-1,i)+Ninit_ais(0,nage-2,i)) + NCome)*exp(-mat_age(nage-1));
          Length_yais_beg(y,nage-1,i) = Linf_yk(0,phi_ik2(i),s)+(L1_yk(0,phi_ik2(i),s)-Linf_yk(0,phi_ik2(i),s))*
            exp(-kappa_yk(0,phi_ik2(i))*(nage-1));
          Length_yais_mid(y,nage-1,i,s) =  Linf_yk(0,phi_ik2(i),s)+(L1_yk(0,phi_ik2(i),s)-Linf_yk(0,phi_ik2(i),s))*
            exp(-0.5*kappa_yk(0,phi_ik2(i),s)*(nage-1));
        } // end sexes
      } // end subareas i
    } // end y == 0

    Type lenstep = 0.0;Type lenslope = 0.0;
    // N- and Nominal Length - at-age for the middle of this year and beginning of next
    for(int s=0;s<2;s++){
      for(int i=0;i<(nspace);i++){
        N_yais_mid(y,0,i,s) = N_yais_beg(y,0,i,s)*exp(-mat_age(0)/2);
        // linear growth below A4 as in synthesis
        if(L1_yk(y,phi_ik2(i),s) < 3){
          Type lenstep = L1_yk(y,phi_ik2(i),s);
          Type lenslope = (L1_yk(y,phi_ik2(i),s) - lenstep) / 3;
        } else if(L1_yk(y,phi_ik2(i),s) >= 3){
          Type lenstep = 3.0;
          Type lenslope = (L1_yk(y,phi_ik2(i),s) - lenstep) / 3;
        }
        for(int a=0;a<4;a++){
          Length_yais_beg(y,a,i,s) = lenstep+lenslope*a;
        } // end linear age
        Length_yais_beg(y,4,i,s) <-  L1_yk(y,phi_ik2(i),s);
        for(int a=0;a<4;a++){
          Length_yais_mid(y,a,i,s) <- Length_yais_beg(y,a,i,s) + (Linf_yk(y,phi_ik2(i),s)-Length_yais_beg(y,a,i,s)*
            (1-exp(-0.5*kappa_yk(y,phi_ik2(i),s))));
        } // end linear age
        for(int a=1;a<(nage-1);a++){
          Type pLeave = 0.0; Type NCome = 0.0;
          for(int j=0;j<(nspace);j++){
            if(i != j){
              pLeave += X_ijas(i,j,a,s);
              NCome += X_ijas(j,i,a,s)*N_yais_beg(y,a,j,s);
            } // end i != j
          } // end subareas j
          N_yais_mid(y,a,i,s) = ((1-pLeave)*N_yais_beg(y,a,i,s) + NCome)*exp(-mat_age(a)/2);
        } // end ages for N
        for(int a=5;a<(nage-1);a++){
          Length_yais_beg(y+1,a,i,s)  = Length_yais_beg(y,a-1,i,s) + (Linf_yk(y,phi_ik2(i),s)-Length_yais_beg(y,a-1,i,s))*
            (1-exp(-kappa_yk(y,phi_ik2(i),s)));
          Length_yais_mid(y,a,i,s)= Length_yais_beg(y,a,i,s) + (Linf_yk(y,phi_ik2(i),s)-Length_yais_beg(y,a,i,s))*
            (1-exp(-0.5*kappa_yk(y,phi_ik2(i),s)));
        } 
        // plus groups
        Type pLeave = 0.0; Type NCome = 0.0;
        for(int j=0;j<(nspace);j++){
          if(i != j){
            pLeave += X_ijas(i,j,nage-1,s);
            NCome += X_ijas(j,i,nage-1)*(N_yais_beg(y,nage-1,j,s) + N_yais_beg(y,nage-2,j,s));
          } // end i != j
        } // end subareas j
       N_yais_mid(y,nage,i,s) =((1-pLeave)*N_yais_beg(y,nage,i,s) + NCome)*exp(-mat_age(nage)/2);
        // plus group weighted average (we already have the numbers at age)
        Length_yais_beg(y,nage-1,i,s) = (N_yais_beg(y,nage-2,i,s)*
          (Length_yais_beg(y,nage-2,i,s)+(Linf_yk(y,phi_ik2(i),s)-Length_yais_beg(y,nage-2,i,s)*(1-exp(-kappa_yk(y,phi_ik2(i),s))))) +
          N_yais_beg(y,nage-1,i,s)*
          (Length_yais_beg(y,nage-1,i,s)+(Linf_yk(y,phi_ik2(i),s)-Length_yais_beg(y,nage-1,i))*(1-exp(-kappa_yk(y,phi_ik2(i),s)))))/
            (N_yais_beg(y,nage-2,i,s) + N_yais_beg(y,nage-1,i),s);
        
        Length_yais_mid(y,nage-1,i,s) = (N_yais_mid(y,nage-2,i,s)*
          (Length_yais_beg(y,nage-2,i,s)+(Linf_yk(y,phi_ik2(i),s)-Length_yais_beg(y,nage-2,i,s)*(1-exp(-0.5*kappa_yk(y,phi_ik2(i),s))))) +
          N_yais_mid(y,nage-1,i,s)*
          (Length_yais_beg(y,nage-1,i,s)+(Linf_yk(y,phi_ik2(i),s)-Length_yais_beg(y,nage-1,i,s))*(1-exp(-0.5*kappa_yk(y,phi_ik2(i),s)))))/
            (N_yais_mid(y,nage-2,i,s) + N_yais_mid(y,nage-1,i),s);
      } // end subareas i
    } // end sexes
    
    // prob of length-at-age
    for(int s=0;s<2;s++){
      for(int i=0;i<(nspace);i++){
        for(int a=1;a<(nage);a++){
          LengthAge_alyis_beg(a,0,y,i,s) = pnorm(Type(1.0),  Length_yais_beg(y,a,i,s), sigmaG_yk(y,phi_ik2(i),s));
          LengthAge_alyis_mid(a,0,y,i,s) = pnorm(Type(1.0),  Length_yais_mid(y,a,i,s), sigmaG_yk(y,phi_ik2(i),s));
          for(int l=1;l<(LBins-1);l++){
            LengthAge_alyis_beg(a,l,y,i,s) = pnorm(Type(l+1),  Length_yais_beg(y,a,i,s), sigmaG_yk(y,phi_ik2(i),s)) -
              pnorm(Type(l),  Length_yais_beg(y,a,i,s), sigmaG_yk(y,phi_ik2(i),s));
            LengthAge_alyis_mid(a,l,y,i,s) = pnorm(Type(l+1),  Length_yais_mid(y,a,i,s), sigmaG_yk(y,phi_ik2(i),s)) -
              pnorm(Type(l),  Length_yais_mid(y,a,i,s), sigmaG_yk(y,phi_ik2(i),s));
          } // end LBins
          LengthAge_alyis_beg(a,LBins-1,y,i,s) = 1-pnorm(Type(LBins-1), Length_yais_beg(y,a,i,s), sigmaG_yk(y,phi_ik2(i),s));
          LengthAge_alyis_mid(a,LBins-1,y,i,s) = 1-pnorm(Type(LBins-1), Length_yais_mid(y,a,i,s), sigmaG_yk(y,phi_ik2(i),s));
        } // end ages
      } // end nspace
    } // end sex
    
      // Catch at beginning of year
      // Hybrid F tuning inputs & temp storage
      Type v1 = 0.7; Type v2 = 30; Type Fmax = 1.5;
      int niter = 50;

      array<Type> catch_afk_TEMP(nage, nfleets_fish, niter+1);
      
      for(int fish_flt =0;fish_flt<(nfleets_fish);fish_flt++){
        catch_yaf_pred.setZero();
        catch_yf_pred.setZero();
        catch_yfi_pred.setZero();
        catch_yaif_pred.setZero();
        
        Type denom = 0;
        for(int s=0;s<2;s++){
          for(int i=0;i<(nspace);i++){
            switch(selType_fish(fish_flt)){
            case 0: // age sel
              for(int a=1;a<(nage);a++){
                denom += phi_if_fish(fish_flt,i)*
                  fsh_slx_yafs(y,a,fish_flt,s)*
                  N_yais_mid(y,a,i)*
                  wtatlen_kab(phi_ik2(i),1)*
                  pow(Length_yais_mid(y,a,i,s),wtatlen_kab(phi_ik2(i),2))+
                  catch_yf_obs(y,fish_flt+1);
              } // end age
              break;
            case 1: // length sel
              for(int a=1;a<(nage);a++){
                for(int l=1;l<(LBins);l++){
                  denom += phi_if_fish(fish_flt,i)*
                    fsh_slx_yafs(y,a,fish_flt,s)*
                    N_yais_mid(y,a,i)*
                    LengthAge_alyis_mid(a,l,y,i,s)*
                    wtatlen_kab(phi_ik2(i),1)*
                    pow(LengthAge_alyis_mid(a,l,y,i,s),wtatlen_kab(phi_ik2(i),2))+
                    catch_yf_obs(y,fish_flt+1);
                } // end length
              } // end age
              break;
            } // end selType_fish
          } // end space
        } // end sex
        F1_yf(y,fish_flt,1) = catch_yf_obs(y, fish_flt)/denom;
        Type latest_guess = F1_yf(y,fish_flt,1);
        
        // k iterations 
        for(int k=2;k<(niter+1);k++){    
          // modify the guess Eq 20
          Type term0 = 1/(1+exp(v2*( latest_guess - v1)));
          Type term1 = latest_guess*term0;
          Type term2 = v1*(1-term0);
          F1_yf(y,fish_flt,k) = -log(1-(term1+term2));
          vector<Type>Z_a_TEMP(nage);

          for(int i=0;i<(nspace);i++){
            switch(selType_fish(fish_flt)){
            case 0: // age sel
              for(int a=0;a<(nage);a++){
                // catch_afk_TEMP(a,fish_flt,k).setZero();
                for(int s=0;s<2;s++){
                  Z_a_TEMP[a] += fsh_slx_yafs(y, a, fish_flt, s)*F1_yf(y,fish_flt,k) + mat_age(a);
                } // end sex for z a temp
                for(int s=0;s<2;s++){
                  catch_afk_TEMP(a,fish_flt,k) +=
                    F1_yf(y,fish_flt,k)/Z_a_TEMP[a]*
                    (1-exp(-Z_a_TEMP[a]))*
                    phi_if_fish(fish_flt,i)*
                    fsh_slx_yafs(y,a,fish_flt,s)*
                    N_yais_mid(y,a,i)*
                    wtatlen_kab(phi_ik2(i),1)*
                    pow(Length_yais_mid(y,a,i,s),wtatlen_kab(phi_ik2(i),2));
                } // end sex
              } // end age
              break;
            case 1: // length sel
              for(int a=1;a<(nage);a++){
                for(int s=0;s<2;s++){
                  Z_a_TEMP[a] += fsh_slx_yafs(y, a, fish_flt, s)*F1_yf(y,fish_flt,k) + mat_age(a);
                } // end sex for z a temp
              } // end age
              for(int l=0;l<(LBins);l++){
                for(int a=0;a<(nage);a++){
                  for(int s=0;s<2;s++){
                    catch_afk_TEMP(a,fish_flt,k) +=
                      F1_yf(y,fish_flt,k)/Z_a_TEMP[a]*
                      (1-exp(-Z_a_TEMP[a]))*
                      phi_if_fish(fish_flt,i)*
                      fsh_slx_yafs(y,l,fish_flt,s)* 
                      N_yais_mid(y,a,i)*
                      LengthAge_alyis_mid(a,l,y,i,s)*
                      wtatlen_kab(phi_ik2(i),1)*
                      pow(mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),2));
                  } // end sex
              } // end age
              } // end lbins
              break;
            } // end selType_fish
          } // end space
          vector<Type>Adj(niter);
          for(int a=0;a<(nage);a++){
            Adj(k) += catch_yf_obs(y,fish_flt+1)/catch_afk_TEMP(a,fish_flt,k);
          }
          // Get new Z given ADJ - need to add discard here
          vector<Type>Z_a_TEMP2(nage);
          Z_a_TEMP2.setZero();
          for(int a=0;a<(nage);a++){
            for(int s=0;s<2;s++){
            Z_a_TEMP2(a) += Adj(k)  *
              fsh_slx_yafs(y, a, fish_flt, s) * F1_yf(y, fish_flt, k) +
              mat_age(a);
            } // end sex
          } // end age
      } // end k iters
      
      
  } // temp yend
  
      //   for(int a=0;a<nage;a++){
      
      //       // make an initial guess for F using obs catch - need to update selex
      //      
      //         (phi_if_fish(fish_flt, i) * N_yais_beg(y,a,i)*wage_catch(a,y) *  selectivity_save(a,y) + catch_yf_obs(y, fish_flt));
      // 
      //       for(int fiter=0; fiter<10;fiter++){
      //         // modify the guess (overwrite)
      // 
      //         term0 = 1/(1+exp(v2*( F1_yf(y,fish_flt) - v1)));
      //         term1 = F1_yf(y,fish_flt)*term0;
      //         term2 = v1*(1-term0);
      //         F1_yf(y,fish_flt) = -log(1-(term1+term2));
      // 
      //       } // end hybrid F iterations
      //       Catch_yaf_est(y,a,fish_flt) = (Freal(a)/(Z(a)))*(1-exp(-Z(a)))* phi_if_fish(fish_flt, i)* N_yais_beg(y,a,i)*wage_catch(a,y); // do this by fleet with phi
      //       CatchN_yaf(y,a,fish_flt) = (Freal(a)/(Z(a)))*(1-exp(-Z(a)))* phi_if_fish(fish_flt, i)* N_yais_beg(y,a,i);// Calculate the catch in kg
      //       Catch_yf_est(y,fish_flt) += Catch_yaf_est(y,a,fish_flt); // sum over the current catch at age
      //       CatchN(y,fish_flt) += CatchN_yaf(y,a,fish_flt);
      //     } // end fishery fleets
      //   } // end ages
      // } // end nspace
        } // temporary yend
    
    
    
  //   
  //   // calculate SSB using N at beginning of year
  // 
  //   for(int i=0;i<(nspace);i++){
  //     for(int a=0;a<nage;a++){ // Loop over ages
  //       SSB_yi(y,i) += N_yais_beg(y,a,i)*wage_ssb(a,y)*0.5; // for storage
  //       for(int k=0;k<(nstocks);k++){  
  //         SSB_yk(y,k) += phi_kik,i)*N_yais_beg(y,a,i)*wage_ssb(a,y)*0.5; // hat
  //       } // end stocks
  //     } // end ages
  //   } // end space
  //   
  //   
  //   // generate recruits (N age = 0) this year based on present SSB
  //   for(int i=0;i<(nspace);i++){
  //     for(int k=0;k<(nstocks);k++){  
  //       // SSB_yk already has summation
  //       R_yk(y,k) = (4*h_k(k)*R_0k(k)*SSB_yk(y,k)
  //                                      /(SSB_0k(k)*(1-h_k(k))+ 
  //         SSB_yk(y,k)*(5*h_k(k)-1)))*exp(-0.5*b(y)*SDR*SDR+tildeR_yk(y,k));
  //     } // end stocks
  //     R_yi(y,i) = R_yk(y,phi_ik2(i))*tau_ki(phi_ik2(i),i)*omega_0ij(i); // downscale to subarea including age-0 movement
  //     N_yais_beg(y,0,i) =  R_yi(y,i); // fill age-0 recruits
  //   } // end space
  //   
  
  //     
  //     // reweight length-at-age based on movement from other stocks
  //     for(int i=0;i<(nspace);i++){
  //       for(int a=1;a<(nage);a++){
  //         Type LCome = 0.0; Type NCome = 0.0;
  //         for(int j=0;j<(nspace);j++){
  //           // sum up other areas applicable to this subarea + age
  //           if(i != j){
  //             LCome += phi_ij(i,j)*N_yais_beg(y,a,j)*Length_yai_beg(y,a,j); // for numerator
  //             NCome += phi_ij(i,j)*N_yais_beg(y,a,j); // for denom
  //           }
  //         } // end subareas j
  //         Length_yai_beg(y+1,a,i) = (N_yais_beg(y,a,i)*Length_yai_beg(y,a,i) + LCome)/(N_yais_beg(y,a,i)+NCome);
  //       } // end ages
  //     } // end subareas i
  //     
  //       

  //   
  
  //   
  //   // Estimate survey biomass at midyear
  //   for(int i=0;i<(nspace);i++){
  //     for(int a=0;a<nage;a++){   
  //       for(int sur_flt =0;sur_flt<(nfleets_surv);sur_flt++){
  //         survey_bio_f_est(y,sur_flt) += surveyselc(a)*wage_survey(a,y)*phi_if_surv(sur_flt,i)*N_yais_mid(y,a,i)*q; // need to include phi matrix to conditionally sum biomass over i 
  //         Nsamp_acomp_f(sur_flt) += surveyselc(a)*phi_if_surv(sur_flt,i)*N_yais_mid(y,a,i); // To use with age comps; may need to change phi to sum acomp surveys
  //       } // end surv fleets
  //     } // end ages
  //   } // end nspace
  //   
  //   
  //   
  //   
  //   // estimate age comps in surveys
  //   // need to include error here
  //   for(int surv_flt_acomp =0;surv_flt_acomp<(nfleets_acomp);surv_flt_acomp++){
  //     if(flag_surv_acomp(y) == 1){ // flag if there is an age measurement this year
  //       for(int i=0;i<(nspace);i++){
  //         for(int a=0;a<(age_maxage-1);a++){ // Loop over other ages
  //           if(a == 1){
  //             // first determine aging error offset
  //             // note that the first row has the a-tilde, the second row has the SD by fleet
  //             acomp_yaf_temp(y,a,surv_flt_acomp) = pnorm(age(a),  age_error(0,a,surv_flt_acomp), age_error(1,a,surv_flt_acomp));
  // 
  //           } else if(a< age_maxage){
  //             acomp_yaf_temp(y,a,surv_flt_acomp) = pnorm(Type(a+1),  age_error(0,a,surv_flt_acomp), age_error(1,a,surv_flt_acomp)) -
  //               pnorm(age(a),  age_error(0,a,surv_flt_acomp), age_error(1,a,surv_flt_acomp));
  //             // survey_acomp_f_est(y,a,surv_flt_acomp) = (surveyselc(a+1)*phi_if_surv(surv_flt_acomp,i)*N_yais_mid(y,a+1,i))/Nsamp_acomp_f(surv_flt_acomp); // estimated comps based on nbeg, should be fleet accrued
  //           } // end else
  //         } // end ages
  //           acomp_yaf_temp(y,age_maxage-1,surv_flt_acomp) = Type(1.0) -
  //             pnorm(Type(age_maxage-1), age_error(0,age_maxage-1,surv_flt_acomp), age_error(1,age_maxage-1,surv_flt_acomp) );
  //           //   
  //           //   survey_acomp_f_est(y,age_maxage-1,surv_flt_acomp) += (surveyselc(a+1)*phi_if_surv(surv_flt_acomp,i)*N_yais_mid(y,a+1,i))/Nsamp_acomp_f(surv_flt_acomp); // placeholder note the indexing on ntot might be off
  //           // 
  //           // 
  //           
  // 
  //       } // end nspace
  //     }  // end flag
  //   } // end acomp survey fleets
  //   
  //   // estimate age comps in catches
  //   // need to include error here
  //   if(flag_catch(y) == 1){ // Flag if  there was a measurement that year
  //     for(int fish_flt =0;fish_flt<(nfleets_fish);fish_flt++){
  //       for(int a=0;a<(nage-1);a++){ // Loop over ages for catch comp
  //         if(a<age_maxage){
  //           catch_acomp_f_est(y,a,fish_flt) = (CatchN_yaf(y,a+1,fish_flt)/CatchN(y,fish_flt)); // Catch comp (1 bc the data starts at age = 1)
  //         }else{
  //           catch_acomp_f_est(y,age_maxage-1,fish_flt) += (CatchN_yaf(y,a+1,fish_flt)/CatchN(y,fish_flt));
  //         } // end else
  //       } // end ages
  //     } // end fish_flt
  //   } // end flag
  // } // END y LOOP
  // 
  // 
  // // LIKELIHOODS //
  // // using namespace density;
  // Type ans_survey=0.0;
  // ////Save the observation model estimates
  // 
  // // Likelihood: survey biomass
  // for(int surv_flt =0;surv_flt<(nfleets_surv);surv_flt++){
  //   for(int y=1;y<tEnd;y++){ // Survey Surveyobs
  //     if(flag_surv_bio(y) == 2){
  //       ans_survey += -dnorm(log(survey_bio_f_est(y,surv_flt)), 
  //                            log(survey_bio_f_obs(y,surv_flt)), 
  //                            SDsurv+survey_err(y), TRUE); // the err also needs to be by flt
  //     } // end survey flag
  //     
  //   } // end y
  // } // end surv_flt
  // 
  // 
  // // Likelihood: catches
  // Type ans_catch = 0.0;
  // for(int fish_flt =0;fish_flt<(nfleets_fish);fish_flt++){
  //   for(int y=0;y<tEnd;y++){ // Total Catches
  //     ans_catch += -dnorm(log(Catch_yf_est(y,fish_flt)+1e-6), log(catch_yf_obs(y,fish_flt)+1e-6), SDcatch, TRUE); 
  //     
  //   }
  // }
  // 
  // REPORT(ans_catch)
  //   
  // // Likelihood: age comps in survey
  // Type ans_survcomp = 0.0;
  // Type ans_catchcomp = 0.0;
  // 
  // 
  // vector<Type>sum1(tEnd);
  // vector<Type>sum2(tEnd);
  // 
  // sum1.setZero();
  // sum2.setZero();
  // 
  // for(int surv_flt_acomp =0;surv_flt_acomp<(nfleets_acomp);surv_flt_acomp++){
  //   for(int y=1;y<tEnd;y++){ // Loop over available years
  //     if(flag_surv_acomp(y) == 1){ // Flag if  there was a measurement that year
  //       for(int a=1;a<age_maxage;a++){ // Loop over other ages (first one is empty for survey)
  //         // NOTE THAT THE survey_acomp_f_obs ARE IN A X y X FLEET, which is not the typical ordering
  //         // need to replace SS_survey(y) with the number of samples from each year x fleet
  //         sum1(y) += lgamma(ss_survey(y)*survey_acomp_f_obs(a,y,surv_flt_acomp)+1);
  //         sum2(y) += lgamma(ss_survey(y)*survey_acomp_f_obs(a,y,surv_flt_acomp) + phi_survey*ss_survey(y)*survey_acomp_f_est(y,a,surv_flt_acomp)) -
  //           lgamma(phi_survey*ss_survey(y)*survey_acomp_f_est(y,a,surv_flt_acomp));
  //       } // end ages
  //       ans_survcomp += lgamma(ss_survey(y)+1)-sum1(y)+lgamma(phi_survey*ss_survey(y))-lgamma(ss_survey(y)+phi_survey*ss_survey(y))+sum2(y);
  //     } // end acomp flag
  //   } // end y
  // } // end survey acomp fleets
  // 
  // 
  // vector<Type>sum3(tEnd);
  // vector<Type>sum4(tEnd);
  // //
  // sum3.setZero();
  // sum4.setZero();
  // 
  // // Likelihood: age comps in catches
  // // need to remove/change this phi_Catch thing and provide fish catch comps by fleet
  // for(int fish_flt =0;fish_flt<(nfleets_fish);fish_flt++){ 
  //   for(int y=1;y<tEnd;y++){ // Loop over available years
  //     if(catch_yf_obs(y,fish_flt)>0){ // only bother if we caught something
  //       if(flag_catch(y) == 1){ // Flag if  there was a measurement that year
  //         for(int a=0;a<(age_catch.rows()-1);a++){ // Loop over ages for catch comp (there are only 15 in obs)
  //           sum3(y) += lgamma(catch_yf_obs(y,fish_flt)*age_catch(a,y)+1);
  //           sum4(y) += lgamma(catch_yf_obs(y,fish_flt)*age_catch(a,y) + phi_catch*catch_yf_obs(y,fish_flt)*catch_acomp_f_est(y,a,fish_flt)) -
  //             lgamma(phi_catch*catch_yf_obs(y,fish_flt)*catch_acomp_f_est(y,a,fish_flt));
  //         } // end ages
  //         ans_catchcomp += lgamma(catch_yf_obs(y,fish_flt)+1)-sum3(y)+lgamma(phi_catch*catch_yf_obs(y,fish_flt))
  //           -lgamma(catch_yf_obs(y,fish_flt)+phi_catch*catch_yf_obs(y,fish_flt))+sum4(y);
  //       } // end catch flag 
  //     } // end if catches > 0
  //   } // end y
  // } // end fish fleets 
  // 
  // // Likelihood: SD Recruitment (hyperprior)
  // Type ans_SDR = 0.0;
  // for(int k=0;k<(nstocks);k++){
  //   for(int y=0;y<(tEnd-1);y++){ // Start y loop
  //     ans_SDR += Type(0.5)*(tildeR_yk(y,k)*tildeR_yk(y,k))/(SDR*SDR)+b(y)*log(SDR*SDR);
  //   }
  // }
  // 
  // // Likelihood: Error for Selectivity
  // Type ans_psel = 0.0;
  // for(int y=0;y<year_sel;y++){ // Start y loop
  //   for(int i=0;i<psel_fish.size();i++){ // Start y loop
  //     ans_psel += Type(0.5)*(PSEL(i,y)*PSEL(i,y))/(sigma_psel*sigma_psel);
  //   }
  // }
  // 
  // // Likelihood: Priors on h and M
  // Type ans_priors = 0.0;
  // for(int y=0;y<(nage);y++){ // Start y loop
  //   for(int i=0;i<(nspace);i++){
  //     for(int a=0;a<(nage-1);a++){ // needs to loop over all years of inits
  //       ans_priors += Type(0.5)*(Ninit_ai(a,i)*Ninit_ai(a,i))/(SDR*SDR);
  //     } // end ages
  //   } // end space
  // } // end y
  // 
  // // ans_priors += -dnorm(logh,log(Type(0.777)),Type(0.113),TRUE);
  // 
  // // Likelihood: Prior on h
  // // ans_priors += -dbeta(h,Bprior,Aprior,TRUE);
  // 
  // if(sum_zero == 1){
  //   ans_priors += ((Type(0.0)-sum(tildeR_yk))*(Type(0.0)-sum(tildeR_yk)))/Type(0.01);
  // }
  // 
  // 
  // // ans_priors += -dnorm(logMinit, log(Type(0.2)), Type(0.1), TRUE);
  // ans_priors += 0.5*pow(logMinit-log(Type(0.2)),2)/Type(0.01);
  // 
  // vector<Type>ans_tot(7);
  // ans_tot(0) = ans_SDR;
  // ans_tot(1) = ans_psel;
  // ans_tot(2) = ans_catch;
  // ans_tot(3) = ans_survey;
  // ans_tot(4) = ans_survcomp;
  // ans_tot(5) = ans_catchcomp;
  // ans_tot(6) = ans_priors;
  // 
  // // Likelihood: TOTAL
  // Type ans = ans_SDR+ans_psel+ans_catch+ans_survey-ans_survcomp-ans_catchcomp+ans_priors;

  // Later Fix F in the likelihood and age comp in catch
  // Type ans = 0.0;
  // Report calculations
  // ADREPORT(logF)
  // ADREPORT(R)
  // ADREPORT(Fyear)
  // ADREPORT(surveyselc)
  // ADREPORT(catchselec)
  // ADREPORT(age_catch)
  // ADREPORT(age_catch_est)
  // ADREPORT(age_survey)
  // ADREPORT(age_survey_est)
  // ADREPORT(ans_tot)
  REPORT(N_0ais)
  // REPORT(SSB_0k)
  //   REPORT(SSB_yk)
  //   REPORT(SSB_0i)
  //   REPORT(SSB_yi)
  //   REPORT(Ninit_ai)
  //   REPORT(Fyear)
  //   REPORT(R_yk)
  //   REPORT(R_yi)
  //   REPORT(R_0k)
  //   REPORT(logR_0k)
  //   REPORT(omega_0ij)
  //   REPORT(tildeR_yk)
  //   REPORT(tildeR_initk)
  //   
  //   REPORT(Catch_yaf_est)
  //   REPORT(CatchN_yaf)
  //   REPORT(ans_tot)
  //   REPORT(Zsave)
  //   REPORT(survey_acomp_f_est)
  //   REPORT(catch_acomp_f_est)
  //   REPORT(flag_sel)
  //   REPORT(PSEL.cols())
  //   REPORT(selectivity_save)
  //   REPORT(surveyselc)
  //   REPORT(Length_yai_beg)
  //   REPORT(LengthAge_alyi_beg)
  //   REPORT(Length_yai_mid)
  //   REPORT(N_yais_beg)
  //   REPORT(survey_bio_f_est)
  //   REPORT(survey_bio_f_obs)
  //   REPORT(N_yais_mid)
  //   REPORT(Nsamp_acomp_f)
  //   return ans;
}

