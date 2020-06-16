// M KAPUR mod N Jacobsen
// Summer 2020 kapurm@uw.edu
#include <TMB.hpp>
#include <iostream>

template <class Type>
vector<Type> cumsum(vector<Type> x) {
  int n = x.size();
  vector<Type> ans(n);
  ans[0] = x[0];
  for (int i = 1; i < n; i++) ans[i] = x[i] + ans[i-1];
  return ans;
}

// TO DO
// correct growth for early ages
// make selex fleet specific
// introduce selectivity estimation
// introduce tuning of F and use Ztuned in Natage
// need error on tau
// M vs myear -- at age or what?
// estimate SDR (currently parameter)
// on pause: make master flag_fleet matrix

template<class Type>
Type objective_function<Type>::operator() ()
{
  // DATA AND PARAMETERS BY CATEGORY //
  
  // structure //
  DATA_INTEGER(nspace); // number of subreas to track
  DATA_INTEGER(nstocks); // number of stocks (demography)
  DATA_INTEGER(nfleets_surv); // number of survey fleets
  DATA_INTEGER(nfleets_fish); //number of fishery fleets
  DATA_INTEGER(nfleets_acomp); // number of age comp fleets
  DATA_INTEGER(nfleets_lcomp); //number of len comp fleets
  DATA_INTEGER(tEnd); // number of years modeled
  
  DATA_ARRAY(phi_if_surv); // turn on/off subareas for survey fleets
  DATA_ARRAY(phi_if_fish); // turn on/off subareas for fishery fleets
  DATA_ARRAY(phi_ik); // 0/1 nesting of subareas i into stocks k (rows)
  DATA_IVECTOR(phi_ik2); // vector stating which subarea (col) belongs to each stock k (value)
  
  DATA_ARRAY(tau_ik); // downscaling from stocks to sub-areas
  
  // biology // 
  // PARAMETER_VECTOR(Rin); // Time varying stuff
  DATA_INTEGER(nage); // Plus group
  DATA_VECTOR(age); // ages
  DATA_VECTOR(Msel); // How mortality scales with age
  DATA_VECTOR(mat_age); // Maturity ogive
  PARAMETER(logMinit); // Natural mortality
  
  // biology storage
  array<Type> Ninit_ai(nage,nspace); // initial numbers at age in subarea, just once
  array<Type> N_0ai(nage, nspace); // numbers in year 0 at age in subarea
  vector<Type> SSB_0k(nstocks); // virgin spawnbio by stock
  vector<Type> SSB_0i(nspace); // virgin spawnbio by subarea
  
  array<Type> N_yai_beg( tEnd+1, nage, nspace); N_yai_beg.setZero(); 
  array<Type> N_yai_mid( tEnd+1, nage, nspace); N_yai_mid.setZero(); 
  array<Type> SSB_yk(tEnd,nstocks);
  array<Type> SSB_yi(tEnd,nspace);
  array<Type> survey_bio_f_est(tEnd,nfleets_surv); // this is actually predicted
  array<Type> Zsave(nage,tEnd);
  
  // movement //
  PARAMETER_VECTOR(omega_0ij); // estimated age-0 movment among areas (used upon recruit gen)
  DATA_ARRAY(omega_ai); // eigenvect of movement between subareas for ages > 0
  DATA_ARRAY(X_ija); // prob trans between subareas at age
  
  // growth //
  DATA_ARRAY(Linf_yk); // sex, stock, year specific
  DATA_ARRAY(L1_yk); // length at age 4 by stock; linear before this
  DATA_ARRAY(kappa_yk);
  DATA_ARRAY(sigmaG_yk); // perhaps turn to parameter later
  DATA_ARRAY(phi_ij); // matrix of whether i,j are from distinct stocks (0 otherwise)
  DATA_INTEGER(LBins); // maximum length bin in CM
  
  array<Type> Length_yai_beg(tEnd+1,nage,nspace); // placeholder for true lengths-at-age
  array<Type> Length_yai_mid(tEnd+1,nage,nspace); // placeholder for true lengths-at-age
  array<Type> LengthAge_alyi_beg(nage,LBins,tEnd+1,nspace); // placeholder for true age-length dist
  array<Type> LengthAge_alyi_mid(nage,LBins,tEnd+1,nspace); // placeholder for true age-length dist
  
  DATA_ARRAY(wage_ssb); // Weight in the beginning of the year
  DATA_ARRAY(wage_catch); // Weight in catch
  DATA_ARRAY(wage_survey); // Weight in survey
  DATA_ARRAY(wage_mid); // Weight in the middle of the year
  
  // repro //
  PARAMETER_VECTOR(logR_0k); // Recruitment at equil by stock
  PARAMETER(logh); // Steepness
  PARAMETER_VECTOR(logh_k); // Steepness by stock
  DATA_INTEGER(sum_zero); // should rec dev's sum to zero?
  DATA_SCALAR(logSDR); // Can it be estimated as a fixed effect?
  
  // repro storage
  vector<Type> R(tEnd);
  array<Type>  R_yk(tEnd,nstocks); // stock-level recruitment (bev-holt)
  array<Type>  R_yi(tEnd,nspace); // subarea-level recruitment (downscaled)
  
  // observations //
  DATA_INTEGER(year_sel);
  DATA_INTEGER(selYear);
  DATA_SCALAR(logQ);
  DATA_VECTOR(b); // bias adjustment factor
  DATA_VECTOR(years);
  DATA_VECTOR(flag_sel);
  DATA_INTEGER(age_maxage); // Last age included in age comps
  DATA_VECTOR(ss_catch); // age comp sample size
  DATA_VECTOR(flag_catch); // Years age was sampled
  DATA_ARRAY(age_catch); // Age comps in catch -- should be by fleet
  
  // Selectivity
  DATA_INTEGER(Smin);
  DATA_INTEGER(Smin_survey);
  DATA_INTEGER(Smax);
  DATA_INTEGER(Smax_survey);
  array<Type>selectivity_save(nage,tEnd);
  // Survey Biomass
  DATA_VECTOR(flag_surv_bio); // Flag if survey occured
  DATA_VECTOR(survey_err);
  
  DATA_VECTOR(survey); // Acoustic survey - vector of obs x year 
  DATA_ARRAY(survey_bio_f_obs); // year x fleets relbio
  PARAMETER(logSDsurv); // Survey uncertainty
  
  // Survey Comps
  DATA_VECTOR(ss_survey); // Age comp sample size
  DATA_VECTOR(flag_surv_acomp); // Were ages sampled this year
  DATA_ARRAY(age_survey); // Age compositions, age x year
  DATA_ARRAY(survey_acomp_f_obs); // Observed survey compositions, age x year x nfleets_acomp {right now just survnfleets}
  array<Type> survey_acomp_f_est(tEnd, age_maxage, nfleets_acomp); 
  vector<Type> Nsamp_acomp_f(nfleets_acomp); // placeholder for number sampled by comp survey (pre dirichlet weighting)
  
  // Survey Selex
  DATA_SCALAR(smul); // Multiplier for survey selectivity
  DATA_SCALAR(sigma_psel); // selectivity SD
  DATA_SCALAR(logphi_survey);
  PARAMETER_VECTOR(psel_fish);
  PARAMETER_VECTOR(psel_surv);
  PARAMETER_VECTOR(F0);
  PARAMETER_ARRAY(PSEL); // Time varying selectivity
  
  // Catches
  DATA_ARRAY(catch_yf_obs); // obs catch by year and fleet
  DATA_SCALAR(logSDcatch); // Error on catch, should be by fleet
  PARAMETER(logphi_catch);
  
  // Catch Comps
  array<Type> catch_acomp_f_est(tEnd, age_maxage, nfleets_fish); // estimated catch comps; uses derived quants
  
  // Catch storage
  array<Type> Catch_yaf_est(tEnd, nage, nfleets_fish);  // estimated catches at age by fleet
  array<Type> CatchN_yaf(tEnd,nage,nfleets_fish);
  array<Type> Catch_yf_est(tEnd,nfleets_fish); // estimated total catches by fleet
  array<Type> CatchN(tEnd,nfleets_fish);
  
  // F tuning storage
  array<Type> Ftuned_yf(tEnd,nfleets_fish); // final tuned fleet and yr specific F
  array<Type> Ftuned_ym(tEnd,nfleets_fish); // summation of Fs nested in mgmt regions
  array<Type> Ztuned_yf(tEnd, nage, nspace); // final tuned subarea and yr specific Z
  
  
  vector<Type> Fyear(tEnd);
  vector<Type> Freal(nage);
  vector<Type> Z(nage);
  vector<Type> pmax_catch_save(tEnd);
  vector<Type> psel_fish_zero = psel_fish;
  
  
  // Priors
  DATA_SCALAR(Bprior);
  DATA_SCALAR(Aprior);
  
  // PARAMETER(logSDR);
  // PARAMETER(logphi_survey);
  
  // Transform out of log space
  Type SDsurv = exp(logSDsurv);
  Type SDcatch = exp(logSDcatch);
  Type SDR = exp(logSDR);
  vector<Type> R_0k = exp(logR_0k);
  Type h = exp(logh);
  vector<Type> h_k = exp(logh_k);
  
  Type Minit = exp(logMinit);
  Type q = exp(logQ);
  Type phi_survey = exp(logphi_survey);
  Type phi_catch = exp(logphi_catch);
  
  //  Minor calculations
  vector<Type> M = Minit*Msel; // Natural mortality
  vector<Type> Myear = M*Msel; // Natural mortality (if we want to change it later)
  vector<Type> Zzero = M; // total mortality without fishing
  vector<Type> logF(tEnd);
  array<Type> tildeR_yk(tEnd,nstocks); // recdevs
  vector<Type> tildeR_initk(nstocks); // recdevs for init
 
  //  END DATA & PARS, BEGIN MODEL //
  
  // recdevs placeholder
  for(int k=0;k<(nstocks);k++){
    for(int time=0;time<(tEnd-1);time++){
      tildeR_yk(time,k) =0;
    }
    tildeR_yk(tEnd-1,k) =0;
  }
  for(int k=0;k<(nstocks);k++){
      tildeR_initk(k) =0;
  }
  
  /// LIKELY MOVE THIS OUT FOR DIFF SELEX HANDLING ----
  
  // selectivity
  // survey
  vector<Type>surveyselc(nage);
  Type pmax = sum(psel_surv);
  Type ptmp = 0;
  
  for(int j=0;j<nage;j++){ // Fix the survey selectivity
    if (age(j) < Smin_survey){
      surveyselc(j) = 0;
    }
    if (age(j) == Smin_survey){
      ptmp = 0;
      surveyselc(j) = exp(ptmp-pmax);
    }
    if ((age(j) > Smin_survey) & (age(j) <= Smax_survey)){
      ptmp = psel_surv(j-Smin_survey-1)+ptmp;
      surveyselc(j) = exp(ptmp-pmax);
    }
    if(age(j) > (Smax_survey)){
      surveyselc(j) = surveyselc(Smax_survey);
    }
  }
  // Fishing mortality
  vector<Type> catchselec(nage);
  Type pmax_catch = max(cumsum(psel_fish));
  Type ptmp_catch = 0;
  //
  for(int time=0;time<tEnd;time++){ // Start time loop
    for (int j=0;j<nage;j++){
      if (age(j) < Smin){
        catchselec(j) = 0;
      }
      if (age(j) == Smin){
        ptmp_catch = 0;
        catchselec(j) = exp(ptmp_catch-pmax_catch);
      }
      if ((age(j) > Smin) & (age(j) <= Smax)){
        ptmp_catch = psel_fish(j-Smin-1)+ptmp_catch;
        catchselec(j) = exp(ptmp_catch-pmax_catch);
      }
      if(age(j) > (Smax_survey)){
        catchselec(j) = catchselec(Smax);
      }
    }
  }
  
  
  // vector<Type> Nzero(nage); // Numbers with no fishing
  //vector<Type>Meq = cumsum(M);
  //vector<Type>Fpope(tEnd);
  // vector<Type>Bmid(nage); // Biomass at the middle of the year
  // // vector<Type>test(3);
  // // test = cumsum(test);
  // REPORT(test)
  //array<Type> PSEL_save(5,)
  
  // Equilibrium Unfished numbers-at-age, subarea (outside of time loop)
  // identical to Ninit except no recdevs
  N_0ai.setZero();   
  for(int k=0;k<(nstocks);k++){
    for(int i=0;i<(nspace);i++){ 
      for(int a=0;a<(nage-1);a++){
        N_0ai(a,i) = 0.5*omega_ai(a,i)*R_0k(k)*tau_ik(k,i)*exp(-(M(a)*age(a))); // compound multiply duh
      }
      // note the A+ group will be in slot A-1
      N_0ai(nage-1,i) = omega_ai(nage-1,i)* N_0ai(nage-2,i)*exp(-(M(nage-2)*age(nage-2))) 
        /(Type(1.0)-exp(-M(nage-1)*age(nage-1)));
    } // end subareas
  } // end stocks
  
  
  // Equilibrium Unfished SSB, stock (outside of time loop)
  for(int i=0;i<(nspace);i++){ 
    for(int a=0;a<nage;a++){ // Loop over ages
      SSB_0i(i) += mat_age(a)*N_0ai(a,i)*0.5;
      for(int k=0;k<(nstocks);k++){ 
        SSB_0k(k) += phi_ik(k,i)*mat_age(a)*N_0ai(a,i)*0.5;
      } // end stocks
    } // end ages
  } // end space
  
  
  // The first year of the simulation is initialized with the following age distribution 
  for(int k=0;k<(nstocks);k++){
    for(int i=0;i<(nspace);i++){
      for(int a=0;a<(nage-1);a++){
        Ninit_ai(a,i) = 0.5* omega_ai(a,i) * tau_ik(k,i) * R_0k(k)* exp(-(M(a)*age(a))) * exp(-0.5*SDR*SDR+tildeR_initk(k));
      } // end ages
      Ninit_ai(nage-1,i) = (omega_ai(nage-1,i) * Ninit_ai(nage-2,i) *
        exp(-M(nage-1)*age(nage-1)))/(Type(1.0)-exp(-(M(nage-1)*age(nage-1)))* exp(-0.5*SDR*SDR+tildeR_initk(k)));
    } // end space
  } // end stocks
  
  
  for(int time=0;time<(tEnd);time++){ // Start time loop
    
    pmax_catch_save(time) = pmax_catch;
    // Take care of selectivity
    
    
    if (flag_sel(time) == 1){
      
      for(int i=0;i<psel_fish.size();i++){
        psel_fish(i) = psel_fish_zero(i)+PSEL(i,time-selYear+1)*sigma_psel; // 27 is the number of years selectivity is calculated PSEL.cols()-1/ time-selYear-
      }
      
      pmax_catch = max(cumsum((psel_fish)));
      pmax_catch_save(time) = pmax_catch;
      
      for(int j=0;j<(nage);j++){ // Fix the Catch selectivity
        if (age(j) == Smin){
          ptmp_catch = 0;
          catchselec(j) = exp(ptmp_catch-pmax_catch);
        }
        if ((age(j) > Smin) & (age(j) <= Smax)){
          ptmp_catch = psel_fish(j-Smin-1)+ptmp_catch;
          catchselec(j) = exp(ptmp_catch-pmax_catch);
        }
        if(age(j) > (Smax_survey)){
          catchselec(j) = catchselec(Smax);
        }
      }
    }
    
    for(int fish_flt=0;fish_flt<(nfleets_fish);fish_flt++){
      Catch_yf_est(time,fish_flt) = 0;
    }
    
    Fyear(time) = F0(time);
    
    
    for(int a=0;a<(nage);a++){ // Loop over other ages
      Freal(a) = Fyear(time)*catchselec(a);
      Z(a) = Freal(a)+Myear(a);
      selectivity_save(a,time) = catchselec(a);
      Zsave(a,time) = Z(a);
    }
    
    
    // model year zero, use last year of Ninit_ai, and equil movement (omega) and downscaling (tau)
    // note we are assuming unfished here as the exponent is M only
    // add length at age initial here so that time loop can be removed
    if (time == 0){  
      for(int k=0;k<(nstocks);k++){
        for(int i=0;i<(nspace);i++){ 
          for(int j=0;j<(nspace);j++){ 
            for(int a=1;a<(nage-1);a++){ // we will fill recruits (a0) later
              Type pLeave = 0.0; Type NCome = 0.0; // reset for new age
              if(i != j){
                pLeave += X_ija(i,j,a); // will do 1-this for proportion which stay
                NCome += X_ija(j,i,a)*Ninit_ai(a,j); // actual numbers incoming
              }
              
              // this is the synthesis syntax; 10 is placeholder for LMIN
              // likely need a lower L1 at age stock-specific and linear before that age
              Length_yai_beg(time,a,i) = Linf_yk(0,phi_ik2(i))+(10-Linf_yk(0,phi_ik2(i)))*
                exp(-kappa_yk(0,phi_ik2(i))*a);
              Length_yai_mid(time,a,i) = Linf_yk(0,phi_ik2(i))+(10-Linf_yk(0,phi_ik2(i)))*
                exp(-0.5*kappa_yk(0,phi_ik2(i))*a);
              
              N_yai_beg(time,a,i) = ((1-pLeave)*Ninit_ai(a,i) + NCome)*exp(-M(a));
            } // end ages
            Type pLeave = 0.0; Type NCome = 0.0; // reset for plusgroup age
            // plus group includes those already at A AND age into A
            if(i != j){
              pLeave += X_ija(i,j,nage-1); 
              NCome += X_ija(j,i,nage-1)*(N_yai_beg(0,nage-1,j) + N_yai_beg(0,nage-2,j)) ; // if M becomes spatial use M_aj here
            }
            N_yai_beg(time,nage-1,i) =  ((1-pLeave)*(N_yai_beg(0,nage-1,i)+N_yai_beg(0,nage-2,i)) + NCome)*exp(-M(nage-1));
           
           Length_yai_beg(time,nage-1,i) = Linf_yk(0,phi_ik2(i))+(10-Linf_yk(0,phi_ik2(i)))*
             exp(-kappa_yk(0,phi_ik2(i))*nage-1);
           Length_yai_mid(time,nage-1,i) = Linf_yk(0,phi_ik2(i))+(10-Linf_yk(0,phi_ik2(i)))*
             exp(-0.5*kappa_yk(0,phi_ik2(i))*nage-1);
          } // end subareas j
        } // end subareas i
      } // end stocks
    } // end time == 0
    
    // calculate SSB using N at beginning of year

    for(int i=0;i<(nspace);i++){
      for(int a=0;a<nage;a++){ // Loop over ages
        SSB_yi(time,i) += N_yai_beg(time,a,i)*wage_ssb(a,time)*0.5; // for storage
        for(int k=0;k<(nstocks);k++){  
          SSB_yk(time,k) += phi_ik(k,i)*N_yai_beg(time,a,i)*wage_ssb(a,time)*0.5; // hat
        } // end stocks
      } // end ages
    } // end space
    
    
    // generate recruits (N age = 0) this year based on present SSB
    for(int i=0;i<(nspace);i++){
      for(int k=0;k<(nstocks);k++){  
        R_yk(time,k) += phi_ik(k,i)*(4*h_k(k)*R_0k(k)*SSB_yk(time,k)
                                       /(SSB_0k(k)*(1-h_k(k))+ 
          SSB_yk(time,k)*(5*h_k(k)-1)))*exp(-0.5*b(time)*SDR*SDR+tildeR_yk(time,k));
      } // end stocks
      R_yi(time,i) = R_yk(time,phi_ik2(i))*0.5; //tau_ik(k,i)*omega_0ij(i); // downscale to subarea including age-0 movement
      N_yai_beg(time,0,i) =  R_yi(time,i); // fill age-0 recruits
    } // end space
    
    // N- and Nominal Length - at-age for the middle of this year and beginning of next
      for(int i=0;i<(nspace);i++){ 
        for(int j=0;j<(nspace);j++){ 
          for(int a=1;a<(nage-1);a++){
            Type pLeave = 0.0; Type NCome = 0.0; 
            if(i != j){
              pLeave += X_ija(i,j,a); 
              NCome += X_ija(j,i,a)*N_yai_beg(time,a,j); 
            }
            N_yai_mid(time,a,i) = N_yai_beg(time,a,i)*exp(-0.4); // no movement?
            N_yai_beg(time+1,a,i) = ((1-pLeave)*N_yai_beg(time,a,i) + NCome)*exp(-0.4); // this exponent needs to be Ztuned eventually
           
           // as in document: next year A1 == this year A0 plus growth
            Length_yai_beg(time+1,a,i) = Length_yai_beg(time,a-1,i) + (Linf_yk(time,phi_ik2(i))-Length_yai_beg(time,a-1,i))*(1-exp(-kappa_yk(time,phi_ik2(i))));
            Length_yai_mid(time,a,i) = Length_yai_beg(time,a,i) + (Linf_yk(time,phi_ik2(i))-Length_yai_beg(time,a,i))*(1-exp(-0.5*kappa_yk(time,phi_ik2(i))));
          } // end ages
          // plus groups
          Type pLeave = 0.0; Type NCome = 0.0; 
          if(i != j){
            pLeave += X_ija(i,j,nage-1); 
            NCome += X_ija(j,i,nage-1)*(N_yai_beg(time,nage-1,j) + N_yai_beg(time,nage-2,j)); 
          }
          N_yai_mid(time,nage-1,i) = N_yai_beg(time,nage-1,i)*exp(-0.4);
          N_yai_beg(time+1,nage-1,i) =   ((1-pLeave)*(N_yai_beg(time,nage-1,i)+N_yai_beg(time,nage-2,i)) + NCome)*exp(-0.4);
          // plus group weighted average (we already have the numbers at age)
          Length_yai_beg(time,nage-1,i) = (N_yai_beg(time,nage-2,i)*
            (Length_yai_beg(time,nage-2,i)+(Linf_yk(time,phi_ik2(i))-Length_yai_beg(time,nage-2,i)*(1-exp(-kappa_yk(time,phi_ik2(i)))))) +
            N_yai_beg(time,nage-1,i)*
            (Length_yai_beg(time,nage-1,i)+(Linf_yk(time,phi_ik2(i))-Length_yai_beg(time,nage-1,i))*(1-exp(-kappa_yk(time,phi_ik2(i))))))/
              (N_yai_beg(time,nage-2,i) + N_yai_beg(time,nage-1,i));
          
          Length_yai_mid(time,nage-1,i) = (N_yai_mid(time,nage-2,i)*
            (Length_yai_beg(time,nage-2,i)+(Linf_yk(time,phi_ik2(i))-Length_yai_beg(time,nage-2,i)*(1-exp(-0.5*kappa_yk(time,phi_ik2(i)))))) +
            N_yai_mid(time,nage-1,i)*
            (Length_yai_beg(time,nage-1,i)+(Linf_yk(time,phi_ik2(i))-Length_yai_beg(time,nage-1,i))*(1-exp(-0.5*kappa_yk(time,phi_ik2(i))))))/
              (N_yai_mid(time,nage-2,i) + N_yai_mid(time,nage-1,i));
        } // end subareas j
      } // end subareas i
      
      // reweight length-at-age based on movement from other stocks
      for(int i=0;i<(nspace);i++){
        for(int a=1;a<(nage);a++){
          Type LCome = 0.0; Type NCome = 0.0;
          for(int j=0;j<(nspace);j++){
            // sum up other areas applicable to this subarea + age
            if(i != j){
              LCome += phi_ij(i,j)*N_yai_beg(time,a,j)*Length_yai_beg(time,a,j); // for numerator
              NCome += phi_ij(i,j)*N_yai_beg(time,a,j); // for denom
            }
          } // end subareas j
          Length_yai_beg(time+1,a,i) = (N_yai_beg(time,a,i)*Length_yai_beg(time,a,i) + LCome)/(N_yai_beg(time,a,i)+NCome);
        } // end ages
      } // end subareas i
      
        
    // prob of length-at-age
    for(int i=0;i<(nspace);i++){
      for(int a=1;a<(nage);a++){
        LengthAge_alyi_beg(a,0,time,i) = pnorm(Type(1.0),  Length_yai_beg(time,a,i), sigmaG_yk(time,phi_ik2(i)));
        LengthAge_alyi_mid(a,0,time,i) = pnorm(Type(1.0),  Length_yai_mid(time,a,i), sigmaG_yk(time,phi_ik2(i)));
        for(int l=1;l<(LBins-1);l++){
          LengthAge_alyi_beg(a,l,time,i) = pnorm(Type(l+1),  Length_yai_beg(time,a,i), sigmaG_yk(time,phi_ik2(i))) -
            pnorm(Type(l),  Length_yai_beg(time,a,i), sigmaG_yk(time,phi_ik2(i)));
          LengthAge_alyi_mid(a,l,time,i) = pnorm(Type(l+1),  Length_yai_mid(time,a,i), sigmaG_yk(time,phi_ik2(i))) -
            pnorm(Type(l),  Length_yai_mid(time,a,i), sigmaG_yk(time,phi_ik2(i)));
          } // end LBins
        LengthAge_alyi_beg(a,LBins-1,time,i) = 1-pnorm(Type(LBins-1), Length_yai_beg(time,a,i), sigmaG_yk(time,phi_ik2(i)));
        LengthAge_alyi_mid(a,LBins-1,time,i) = 1-pnorm(Type(LBins-1), Length_yai_mid(time,a,i), sigmaG_yk(time,phi_ik2(i)));
      } // end ages
    } // end nspace
    
    
    // Catch at beginning of year
    // Hybrid F tuning inputs
    // Type Fmax = 3.5;
    // Type v1=0.99; //corresponds to an Fmax of 3
    // Type v1=0.865; //corresponds to an Fmax of 2
    // Type v1=0.95;
    // Type v2=30;
    // Type F_no_inc=7;
    // Type term0 = 0.0;
    // Type term1 = 0.0;
    // Type term2 = 0.0;
    
    // array<Type> F1_yf(tEnd,nfleets_fish);
    
    
    // for(int i=0;i<(nspace);i++){
    //   for(int a=0;a<nage;a++){
    //     for(int fish_flt =0;fish_flt<(nfleets_fish);fish_flt++){
    //       // make an initial guess for F using obs catch - need to update selex
    //       F1_yf(time,fish_flt) = catch_yf_obs(time, fish_flt)/
    //         (phi_if_fish(fish_flt, i) * N_yai_beg(time,a,i)*wage_catch(a,time) *  selectivity_save(a,time) + catch_yf_obs(time, fish_flt));
    //       
    //       for(int fiter=0; fiter<10;fiter++){
    //         // modify the guess (overwrite)
    //         
    //         term0 = 1/(1+exp(v2*( F1_yf(time,fish_flt) - v1)));
    //         term1 = F1_yf(time,fish_flt)*term0;
    //         term2 = v1*(1-term0);
    //         F1_yf(time,fish_flt) = -log(1-(term1+term2));
    //         
    //       } // end hybrid F iterations
    //       Catch_yaf_est(time,a,fish_flt) = (Freal(a)/(Z(a)))*(1-exp(-Z(a)))* phi_if_fish(fish_flt, i)* N_yai_beg(time,a,i)*wage_catch(a,time); // do this by fleet with phi
    //       CatchN_yaf(time,a,fish_flt) = (Freal(a)/(Z(a)))*(1-exp(-Z(a)))* phi_if_fish(fish_flt, i)* N_yai_beg(time,a,i);// Calculate the catch in kg
    //       Catch_yf_est(time,fish_flt) += Catch_yaf_est(time,a,fish_flt); // sum over the current catch at age 
    //       CatchN(time,fish_flt) += CatchN_yaf(time,a,fish_flt);
    //     } // end fishery fleets
    //   } // end ages
    // } // end nspace
    
    
    // Estimate survey biomass at midyear
    for(int i=0;i<(nspace);i++){
      for(int a=0;a<nage;a++){   
        for(int sur_flt =0;sur_flt<(nfleets_surv);sur_flt++){
          survey_bio_f_est(time,sur_flt) += surveyselc(a)*wage_survey(a,time)*phi_if_surv(sur_flt,i)*N_yai_mid(time,a,i)*q; // need to include phi matrix to conditionally sum biomass over i 
          Nsamp_acomp_f(sur_flt) += surveyselc(a)*phi_if_surv(sur_flt,i)*N_yai_mid(time,a,i); // To use with age comps; may need to change phi to sum acomp surveys
        } // end surv fleets
      } // end ages
    } // end nspace
    
    
    // estimate age comps in surveys
    // need to include error here
    for(int surv_flt_acomp =0;surv_flt_acomp<(nfleets_acomp);surv_flt_acomp++){
      if(flag_surv_acomp(time) == 1){ // flag if there is an age measurement this year
        for(int i=0;i<(nspace);i++){
          for(int a=0;a<(nage-1);a++){ // Loop over other ages
            if(a< age_maxage){
              survey_acomp_f_est(time,a,surv_flt_acomp) = (surveyselc(a+1)*phi_if_surv(surv_flt_acomp,i)*N_yai_mid(time,a+1,i))/Nsamp_acomp_f(surv_flt_acomp); // estimated comps based on nbeg, should be fleet accrued
            }else{
              survey_acomp_f_est(time,age_maxage-1,surv_flt_acomp) += (surveyselc(a+1)*phi_if_surv(surv_flt_acomp,i)*N_yai_mid(time,a+1,i))/Nsamp_acomp_f(surv_flt_acomp); // placeholder note the indexing on ntot might be off
            } // end else
          } // end ages
        } // end nspace
      }  // end flag
    } // end acomp survey fleets
    
    // estimate age comps in catches
    // need to include error here
    if(flag_catch(time) == 1){ // Flag if  there was a measurement that year
      for(int fish_flt =0;fish_flt<(nfleets_fish);fish_flt++){
        for(int a=0;a<(nage-1);a++){ // Loop over ages for catch comp
          if(a<age_maxage){
            catch_acomp_f_est(time,a,fish_flt) = (CatchN_yaf(time,a+1,fish_flt)/CatchN(time,fish_flt)); // Catch comp (1 bc the data starts at age = 1)
          }else{
            catch_acomp_f_est(time,age_maxage-1,fish_flt) += (CatchN_yaf(time,a+1,fish_flt)/CatchN(time,fish_flt));
          } // end else
        } // end ages
      } // end fish_flt
    } // end flag
  } // END TIME LOOP
  
  
  // LIKELIHOODS //
  // using namespace density;
  Type ans_survey=0.0;
  ////Save the observation model estimates
  
  // Likelihood: survey biomass
  for(int surv_flt =0;surv_flt<(nfleets_surv);surv_flt++){
    for(int time=1;time<tEnd;time++){ // Survey Surveyobs
      if(flag_surv_bio(time) == 2){
        ans_survey += -dnorm(log(survey_bio_f_est(time,surv_flt)), 
                             log(survey_bio_f_obs(time,surv_flt)), 
                             SDsurv+survey_err(time), TRUE); // the err also needs to be by flt
      } // end survey flag
      
    } // end time
  } // end surv_flt
  
  
  // Likelihood: catches
  Type ans_catch = 0.0;
  for(int fish_flt =0;fish_flt<(nfleets_fish);fish_flt++){
    for(int time=0;time<tEnd;time++){ // Total Catches
      ans_catch += -dnorm(log(Catch_yf_est(time,fish_flt)+1e-6), log(catch_yf_obs(time,fish_flt)+1e-6), SDcatch, TRUE); 
      
    }
  }
  
  REPORT(ans_catch)
    
  // Likelihood: age comps in survey
  Type ans_survcomp = 0.0;
  Type ans_catchcomp = 0.0;
  
  
  vector<Type>sum1(tEnd);
  vector<Type>sum2(tEnd);
  
  sum1.setZero();
  sum2.setZero();
  
  for(int surv_flt_acomp =0;surv_flt_acomp<(nfleets_acomp);surv_flt_acomp++){
    for(int time=1;time<tEnd;time++){ // Loop over available years
      if(flag_surv_acomp(time) == 1){ // Flag if  there was a measurement that year
        for(int a=1;a<age_maxage;a++){ // Loop over other ages (first one is empty for survey)
          // NOTE THAT THE survey_acomp_f_obs ARE IN A X TIME X FLEET, which is not the typical ordering
          // need to replace SS_survey(time) with the number of samples from each year x fleet
          sum1(time) += lgamma(ss_survey(time)*survey_acomp_f_obs(a,time,surv_flt_acomp)+1);
          sum2(time) += lgamma(ss_survey(time)*survey_acomp_f_obs(a,time,surv_flt_acomp) + phi_survey*ss_survey(time)*survey_acomp_f_est(time,a,surv_flt_acomp)) -
            lgamma(phi_survey*ss_survey(time)*survey_acomp_f_est(time,a,surv_flt_acomp));
        } // end ages
        ans_survcomp += lgamma(ss_survey(time)+1)-sum1(time)+lgamma(phi_survey*ss_survey(time))-lgamma(ss_survey(time)+phi_survey*ss_survey(time))+sum2(time);
      } // end acomp flag
    } // end time
  } // end survey acomp fleets
  
  
  vector<Type>sum3(tEnd);
  vector<Type>sum4(tEnd);
  //
  sum3.setZero();
  sum4.setZero();
  
  // Likelihood: age comps in catches
  // need to remove/change this phi_Catch thing and provide fish catch comps by fleet
  for(int fish_flt =0;fish_flt<(nfleets_fish);fish_flt++){ 
    for(int time=1;time<tEnd;time++){ // Loop over available years
      if(catch_yf_obs(time,fish_flt)>0){ // only bother if we caught something
        if(flag_catch(time) == 1){ // Flag if  there was a measurement that year
          for(int a=0;a<(age_catch.rows()-1);a++){ // Loop over ages for catch comp (there are only 15 in obs)
            sum3(time) += lgamma(catch_yf_obs(time,fish_flt)*age_catch(a,time)+1);
            sum4(time) += lgamma(catch_yf_obs(time,fish_flt)*age_catch(a,time) + phi_catch*catch_yf_obs(time,fish_flt)*catch_acomp_f_est(time,a,fish_flt)) -
              lgamma(phi_catch*catch_yf_obs(time,fish_flt)*catch_acomp_f_est(time,a,fish_flt));
          } // end ages
          ans_catchcomp += lgamma(catch_yf_obs(time,fish_flt)+1)-sum3(time)+lgamma(phi_catch*catch_yf_obs(time,fish_flt))
            -lgamma(catch_yf_obs(time,fish_flt)+phi_catch*catch_yf_obs(time,fish_flt))+sum4(time);
        } // end catch flag 
      } // end if catches > 0
    } // end time
  } // end fish fleets 
  
  // Likelihood: SD Recruitment (hyperprior)
  Type ans_SDR = 0.0;
  for(int k=0;k<(nstocks);k++){
    for(int time=0;time<(tEnd-1);time++){ // Start time loop
      ans_SDR += Type(0.5)*(tildeR_yk(time,k)*tildeR_yk(time,k))/(SDR*SDR)+b(time)*log(SDR*SDR);
    }
  }

  // Likelihood: Error for Selectivity
  Type ans_psel = 0.0;
  for(int time=0;time<year_sel;time++){ // Start time loop
    for(int i=0;i<psel_fish.size();i++){ // Start time loop
      ans_psel += Type(0.5)*(PSEL(i,time)*PSEL(i,time))/(sigma_psel*sigma_psel);
    }
  }
  
  // Likelihood: Priors on h and M
  Type ans_priors = 0.0;
  for(int time=0;time<(nage);time++){ // Start time loop
    for(int i=0;i<(nspace);i++){
      for(int a=0;a<(nage-1);a++){ // needs to loop over all years of inits
        ans_priors += Type(0.5)*(Ninit_ai(a,i)*Ninit_ai(a,i))/(SDR*SDR);
      } // end ages
    } // end space
  } // end time
  
  // ans_priors += -dnorm(logh,log(Type(0.777)),Type(0.113),TRUE);
  
  // Likelihood: Prior on h
  ans_priors += -dbeta(h,Bprior,Aprior,TRUE);
  
  if(sum_zero == 1){
    ans_priors += ((Type(0.0)-sum(tildeR_yk))*(Type(0.0)-sum(tildeR_yk)))/Type(0.01);
  }
  
  
  // ans_priors += -dnorm(logMinit, log(Type(0.2)), Type(0.1), TRUE);
  ans_priors += 0.5*pow(logMinit-log(Type(0.2)),2)/Type(0.01);

  vector<Type>ans_tot(7);
  ans_tot(0) = ans_SDR;
  ans_tot(1) = ans_psel;
  ans_tot(2) = ans_catch;
  ans_tot(3) = ans_survey;
  ans_tot(4) = ans_survcomp;
  ans_tot(5) = ans_catchcomp;
  ans_tot(6) = ans_priors;
  
  // Likelihood: TOTAL
  Type ans = ans_SDR+ans_psel+ans_catch+ans_survey-ans_survcomp-ans_catchcomp+ans_priors;

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
   REPORT(SSB_0k)
    REPORT(SSB_yk)
    REPORT(SSB_0i)
    REPORT(SSB_yi)
    REPORT(Ninit_ai)
    REPORT(Fyear)
    REPORT(R_yk)
    REPORT(R_yi)
    REPORT(R_0k)
    REPORT(logR_0k)
    REPORT(omega_0ij)
    REPORT(tildeR_yk)
    REPORT(tildeR_initk)
    REPORT(N_0ai)
    REPORT(Catch_yaf_est)
    REPORT(CatchN_yaf)
    REPORT(ans_tot)
    REPORT(Zsave)
    REPORT(survey_acomp_f_est)
    REPORT(catch_acomp_f_est)
    REPORT(flag_sel)
    REPORT(PSEL.cols())
    REPORT(selectivity_save)
    REPORT(surveyselc)
    REPORT(Length_yai_beg)
    REPORT(LengthAge_alyi_beg)
    REPORT(Length_yai_mid)
    REPORT(N_yai_beg)
    REPORT(survey_bio_f_est)
    REPORT(survey_bio_f_obs)
    REPORT(N_yai_mid)
    REPORT(Nsamp_acomp_f)
    return ans;
}

