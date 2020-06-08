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
// when ready disable the original n_beg etc and remove the numbers
// add in growth and movement properly
// add spatial R init and R (likely need phi_ik)
// make selex fleet specific
// introduce tuning of F


template<class Type>
Type objective_function<Type>::operator() ()
{
  // Data input
  DATA_INTEGER(nspace); // number of subreas to track
  DATA_INTEGER(nstocks); // number of stocks (demography)
  
  DATA_INTEGER(nfleets_surv); // number of survey fleets
  DATA_INTEGER(nfleets_fish); //number of fishery fleets
  DATA_INTEGER(nfleets_acomp); // number of age comp fleets
  DATA_INTEGER(nfleets_lcomp); //number of len comp fleets
  
  DATA_ARRAY(wage_ssb); // Weight in the beginning of the year
  DATA_ARRAY(wage_catch); // Weight in catch
  DATA_ARRAY(wage_survey); // Weight in survey
  DATA_ARRAY(wage_mid); // Weight in the middle of the year
  // Age
  DATA_INTEGER(nage); // Plus group
  DATA_INTEGER(sum_zero); // should rec dev's sum to zero?
  DATA_VECTOR(age); // ages
  DATA_INTEGER(tEnd);
  DATA_INTEGER(year_sel);
  DATA_INTEGER(selYear);
  DATA_SCALAR(logQ);
  DATA_VECTOR(b); // bias adjustment factor
  DATA_VECTOR(years);
  DATA_VECTOR(flag_sel);
  // // Selectivity
  DATA_INTEGER(Smin);
  DATA_INTEGER(Smin_survey);
  DATA_INTEGER(Smax);
  DATA_INTEGER(Smax_survey);
  // // // Survey
  DATA_VECTOR(survey); // Acoustic survey - vector of obs x year 
  DATA_ARRAY(survey2); // year x fleets relbio
  
  DATA_ARRAY(phi_if_surv); // turn on/off subareas for survey fleets
  DATA_ARRAY(phi_if_fish); // turn on/off subareas for fishery fleets
  DATA_ARRAY(phi_ik); // nesting of subareas i into stocks k (rows)
  DATA_ARRAY(tau_ik); // downscaling from stocks to sub-areas
  
  DATA_VECTOR(flag_surv_bio); // Flag if survey occured
  DATA_VECTOR(survey_err);
  DATA_VECTOR(ss_survey); // Age comp sample size
  DATA_VECTOR(flag_surv_acomp); // Were ages sampled this year
  DATA_ARRAY(age_survey); // Age compositions, age x year
  DATA_ARRAY(age_survey2); // Age compositions, age x year x nfleets_acomp {right now just survnfleets}
  
  DATA_INTEGER(age_maxage); // Last age included in age comps
  DATA_SCALAR(smul); // Multiplier for survey selectivity
  // Catches
  DATA_VECTOR(Catchobs); // Total catch
  DATA_ARRAY(Catchobs2); // catch by fleet
  
  DATA_VECTOR(ss_catch); // age comp sample size
  DATA_VECTOR(flag_catch); // Years age was sampled
  DATA_ARRAY(age_catch); // Age comps
  
  DATA_SCALAR(logSDcatch); // Error on catch
  DATA_SCALAR(logSDR); // Can it be estimated as a fixed effect?
  DATA_SCALAR(sigma_psel); // selectivity SD
  // Mortality
  DATA_VECTOR(Msel); // How mortality scales with age
  DATA_VECTOR(Matsel); // Maturity ogive
  

  // Priors
  DATA_SCALAR(Bprior);
  DATA_SCALAR(Aprior);
  // Time parameters Parameter integers
  PARAMETER(logRinit); // Recruitment at
  PARAMETER(logh); // Steepness
  PARAMETER(logMinit); // Natural mortality
  PARAMETER(logSDsurv); // Survey uncertainty
  //PARAMETER(logSDR);
  PARAMETER(logphi_catch);
  // PARAMETER(logphi_survey);
  DATA_SCALAR(logphi_survey);
  PARAMETER_VECTOR(psel_fish);
  PARAMETER_VECTOR(psel_surv);
  PARAMETER_VECTOR(initN);
  PARAMETER_VECTOR(Rin); // Time varying stuff
  PARAMETER_ARRAY(PSEL); // Time varying selectivity
  PARAMETER_VECTOR(F0);
  // Transform out of log space
  Type SDsurv = exp(logSDsurv);
  Type SDcatch = exp(logSDcatch);
  Type SDR = exp(logSDR);
  Type Rinit = exp(logRinit);
  Type h = exp(logh);
  Type Minit = exp(logMinit);
  Type q = exp(logQ);
  Type phi_survey = exp(logphi_survey);
  Type phi_catch = exp(logphi_catch);
  //  Minor calculations
  vector<Type> M = Minit*Msel; // Natural mortality
  vector<Type> logF(tEnd);
  vector<Type> logR(tEnd);
  // Vectors for saving stuff
  vector<Type> R(tEnd);
  array<Type>  R_k(tEnd,nstocks); // stock-level recruitment (bev-holt)
  array<Type>  R_i(tEnd,nspace); // subarea-level recruitment (downscaled)
  
  
  array<Type> CatchAge(nage,tEnd); // original
  array<Type> CatchAge2(tEnd, nage, nfleets_fish);
  
  array<Type> CatchNAge(nage,tEnd);
  array<Type> CatchNAge2(tEnd,nage,nfleets_fish);
  
  
  for(int j=0;j<(tEnd-1);j++){
    logR(j)=Rin(j);
  }
  logR(tEnd-1) = 0;
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
  array<Type>   Nzero3(nage, nspace); Nzero3.setZero();   // grant's approach
  vector<matrix<Type> > Nzero2(nspace); // create this many matrices within a vector
  vector<Type> Nzero(nage); // Numbers with no fishing
  //vector<Type>Meq = cumsum(M);
  
  for(int i=0;i<(nspace);i++){ // there are nspace+1 slots, the last one is for total
    Nzero(0) = Rinit;
    Nzero3(0,i) = Rinit;
    for(int a=1;a<(nage-1);a++){
      Nzero(a) = Rinit * exp(-(M(a)*age(a)));
      Nzero3(a,i) = Rinit * exp(-(M(a)*age(a)));
      
    }
    Nzero(nage-1) = (Rinit*exp(-(M(nage-2)*age(nage-1))))/(Type(1.0)-exp(-M(nage-1))); // note the A+ will be in slot A-1
    Nzero3(nage-1,i) = (Rinit*exp(-(M(nage-2)*age(nage-1))))/(Type(1.0)-exp(-M(nage-1))); // note the A+ will be in slot A-1
    
    Nzero2(i) = Nzero;
  } // end subareas
  
  
  array<Type> SSBage(nage);
  array<Type> Catchinit(nage);
  array<Type>selectivity_save(nage,tEnd);
  Type SSBzero = 0;
  vector<Type> SSBzero2(nspace);
  vector<Type> Zzero = M;
  
  for(int i=0;i<(nspace);i++){ 
    for(int a=0;a<nage;a++){ // Loop over ages
      SSBzero += Matsel(a)*Nzero(a)*0.5;// original
      SSBzero2(i) += Matsel(a)*Nzero3(a,i)*0.5;
    } 
  }
  // Run the initial distribution
  REPORT(SSBzero);
  REPORT(SSBzero2);
  
  // Type SSBinit = 0;
  //
  // // for(int i=0;i<(nage);i++){ // Loop over other ages
  // //     Ninit(i,0) = Nzero(i);
  // //   }
  // Ninit(0) = Rinit;
  // for(int i=1;i<(nage-1);i++){
  //     Ninit(i) = Rinit * exp(-(M(i)*age(i)))*exp(-0.5*0*SDR*SDR+initN(i-1));
  //   }
  // //
  // Ninit(nage-1) = Rinit*exp(-(M(nage-1)*age(nage-1)))/(Type(1.0)-exp(-M(nage-1)))*exp(-0.5*0*SDR*SDR+initN(nage-2));
  //
  // for(int i=0;i<(nage);i++){ // Loop over other ages
  //   SSBinit += Ninit(i)*Matsel(i)*0.5;
  // }
  
  
  array<Type>Catch(tEnd,nfleets_fish);
  // vector<Type>Catch(tEnd); //original
  array<Type>CatchN(tEnd,nfleets_fish);
  // vector<Type>CatchN(tEnd); // original
  
  matrix<Type> N_mid(nage,tEnd+1);//previously array
  matrix<Type> N_beg(nage,tEnd+1); //previously array
  vector<matrix<Type> > N_beg2(nspace); 
  vector<matrix<Type> > N_mid2(nspace); 
  array<Type>   N_beg3( tEnd+1, nage, nspace); N_beg3.setZero(); 
  array<Type>   N_mid3( tEnd+1, nage, nspace); N_mid3.setZero(); 
  
  // N_mid2.setZero();
  // }
  // // Run the model over time
  array<Type> SSB(tEnd);
  array<Type> SSB2(tEnd,nspace);
  array<Type> Surveyobs(tEnd); // Survey observed Surveyobs
  array<Type> surv_pred(tEnd,nfleets_surv); // this is actually predicted
  
  // array<Type>Surveyobs_tot(tEnd); // Total Surveyobs over age 2
  // array<Type>surv_pred_tot(tEnd); // Total Surveyobs over age 2
  
  array<Type>age_survey_est(age_maxage,tEnd);
  array<Type>age_survey_est2(tEnd, age_maxage, nfleets_acomp); 
  
  array<Type>age_catch_est(age_maxage,tEnd);
  array<Type>age_catch_est2(tEnd, age_maxage, nspace);
  
  array<Type>Zsave(nage,tEnd);
  //
  age_survey_est.setZero();
  age_catch_est.setZero();
  Catch.setZero();
  CatchN.setZero();
  vector<Type> Myear = M*Msel; // Natural mortality (if we want to change it later)
  //
  vector<Type> Fyear(tEnd);
  vector<Type> Freal(nage);
  vector<Type> Z(nage);
  vector<Type>pmax_catch_save(tEnd);
  vector<Type>psel_fish_zero = psel_fish;
  vector<Type>Catchsave(tEnd);
  //vector<Type>Fpope(tEnd);
  
  // vector<Type>Bmid(nage); // Biomass at the middle of the year
  // // vector<Type>test(3);
  // // test = cumsum(test);
  // REPORT(test)
  //array<Type> PSEL_save(5,)
  
  for(int time=0;time<(tEnd);time++){ // Start time loop
    
    Type Ntot_survey = 0;
    vector<Type>Ntot_survey2(nfleets_acomp);
    pmax_catch_save(time) = pmax_catch;
    // Take care of selectivity
    REPORT(flag_sel)
      REPORT(PSEL.cols())
      
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
        Catch(time,fish_flt) = 0;
      }
      Fyear(time) = F0(time);
      
      
      if (time == 0){ // YEAR ZERO
        for(int i=0;i<(nspace);i++){ 
          N_beg.setZero();
          N_beg2(i).setZero();
          for(int a=1;a<(nage-1);a++){
            N_beg(a,time) = Rinit * exp(-0.5*0*SDR*SDR+initN(a-1))*exp(-Myear(a)*age(a));
            N_beg3(time,a,i) = Rinit * exp(-0.5*0*SDR*SDR+initN(a-1))*exp(-Myear(a)*age(a));
          }
          N_beg(nage-1, time) = Rinit * exp(-0.5*0*SDR*SDR+initN(nage-2)) * exp(-Myear(nage-1) * age(nage-1)) / (1 - exp(-Myear(nage-1)));
          N_beg2(i) = N_beg;
          N_beg3(time,nage-1,i) =  Rinit * exp(-0.5*0*SDR*SDR+initN(nage-2)) * exp(-Myear(nage-1) * age(nage-1)) / (1 - exp(-Myear(nage-1)));
          
        } // end subareas
      } // end time == 0
      
      for(int i=0;i<(nspace);i++){ 
        for(int a=0;a<nage;a++){ // Loop over ages
          SSB(time) += N_beg(a,time)*wage_ssb(a,time)*0.5; // hat
          SSB2(time,i) += N_beg3(time,a,i)*wage_ssb(a,time)*0.5; // hat
          // SSB2(time,i) += N_beg2(i)(a,time)*wage_ssb(a,time)*0.5; // hat
        }
      }
      
      for(int a=0;a<(nage);a++){ // Loop over other ages
        Freal(a) = Fyear(time)*catchselec(a);
        Z(a) = Freal(a)+Myear(a);
        selectivity_save(a,time) = catchselec(a);
        Zsave(a,time) = Z(a);
      }
      
      for(int i=0;i<(nspace);i++){
        R(time) = (4*h*Rinit*SSB(time)/(SSBzero*(1-h)+ SSB(time)*(5*h-1)))*exp(-0.5*b(time)*SDR*SDR+logR(time));
        for(int k=0;k<(nstocks);k++){
          R_k(time,k) += phi_ik(k,i)*(4*h*Rinit*SSB2(time,i)/(SSBzero2(i)*(1-h)+ SSB2(time,i)*(5*h-1)))*exp(-0.5*b(time)*SDR*SDR+logR(time));
          R_i(time,i) = R_k(time,k)*tau_ik(k,i); // downscale to subarea
        } // end stocks
        N_beg2(i)(0,time) = R(time);
        N_beg3(time,0,i) =  R_i(time,i);
      } // end space
      //Type smul = Type(0.58);
      for(int i=0;i<(nspace);i++){
        // Catch(time,i) = 0;
        for(int a=0;a<(nage-1);a++){ // Loop over other ages
          N_mid(a,time) =  N_beg2(i)(a,time)*exp(-Z(a)*smul);
          N_mid3(time,a,i) = N_beg3(time,a,i)*exp(-Z(a)*smul);
          N_beg2(i)(a+1,time+1) =  N_beg2(i)(a,time)*exp(-Z(a));
          N_beg3(time+1,a+1,i) =  N_beg3(time,a,i)*exp(-Z(a));
          
        }
        N_mid2(i) = N_mid;
        // Plus group
        N_mid2(i)(nage-1, time) =  N_beg2(i)(nage-2,time)*exp(-Z(nage-2)*0.5)+ N_beg2(i)(nage-1,time)*exp(-Z(nage-1)*smul);
        N_mid3(time,nage-1,i) =  N_beg3(time,nage-2,i)*exp(-Z(nage-2)*0.5)+ N_beg3(time,nage-1,i)*exp(-Z(nage-1)*smul);
        N_beg2(i)(nage-1, time+1) =  N_beg2(i)(nage-2,time)*exp(-Z(nage-2))+ N_beg2(i)(nage-1,time)*exp(-Z(nage-1));
        N_beg3(time+1,nage-1,i) =  N_beg3(time,nage-2,i)*exp(-Z(nage-2))+ N_beg3(time,nage-1,i)*exp(-Z(nage-1));
        
      }
      
      for(int i=0;i<(nspace);i++){
        
        for(int a=0;a<nage;a++){
          
          for(int fish_flt =0;fish_flt<(nfleets_fish);fish_flt++){
            
            CatchAge(a,time)= (Freal(a)/(Z(a)))*(1-exp(-Z(a)))* N_beg2(i)(a,time)*wage_catch(a,time);// Calculate the catch in kg
            CatchAge2(time,a,fish_flt) = (Freal(a)/(Z(a)))*(1-exp(-Z(a)))* phi_if_fish(fish_flt, i)* N_beg3(time,a,i)*wage_catch(a,time); // do this by fleet with phi
            
            CatchNAge(a,time) = (Freal(a)/(Z(a)))*(1-exp(-Z(a)))* N_beg2(i)(a,time);// Calculate the catch in kg
            CatchNAge2(time,a,fish_flt) = (Freal(a)/(Z(a)))*(1-exp(-Z(a)))* phi_if_fish(fish_flt, i)* N_beg3(time,a,i);// Calculate the catch in kg
            
            // Catch(time,i) += CatchAge(a,time); // sum over the current catch at age
            Catch(time,fish_flt) += CatchAge2(time,a,fish_flt); // sum over the current catch at age 
            
            // CatchN(time,i) += CatchNAge(a,time);
            CatchN(time,fish_flt) += CatchNAge2(time,a,fish_flt);
          }
          
          for(int sur_flt =0;sur_flt<(nfleets_surv);sur_flt++){
            
            Surveyobs(time) += surveyselc(a)*wage_survey(a,time)*N_mid(a,time)*q; 
            surv_pred(time,sur_flt) += surveyselc(a)*wage_survey(a,time)*phi_if_surv(sur_flt,i)*N_mid3(time,a,i)*q; // need to include phi matrix to conditionally sum biomass over i 
            Ntot_survey += surveyselc(a)*N_mid(a,time); // To use with age comps
            Ntot_survey2(sur_flt) += surveyselc(a)*phi_if_surv(sur_flt,i)*N_mid3(time,a,i); // To use with age comps; may need to change phi to sum acomp surveys
            
          } // end fleets
        } // end ages
      } // end nspace
      
      REPORT(Ntot_survey2)
        
        for(int surv_flt_acomp =0;surv_flt_acomp<(nfleets_acomp);surv_flt_acomp++){
          if(flag_surv_acomp(time) == 1){ // flag if there is an age measurement this year
            for(int i=0;i<(nspace);i++){
              for(int a=0;a<(nage-1);a++){ // Loop over other ages
                if(a< age_maxage){
                  age_survey_est(a,time) = (surveyselc(a+1)*N_mid(a+1,time))/Ntot_survey;
                  age_survey_est2(time,a,surv_flt_acomp) = (surveyselc(a+1)*phi_if_surv(surv_flt_acomp,i)*N_mid3(time,a+1,i))/Ntot_survey2(surv_flt_acomp); // estimated comps based on nbeg, should be fleet accrued
                  
                }else{
                  age_survey_est(age_maxage-1,time) += (surveyselc(i+1)*N_mid(i+1,time))/Ntot_survey;
                  age_survey_est2(time,age_maxage-1,surv_flt_acomp) += (surveyselc(a+1)*phi_if_surv(surv_flt_acomp,i)*N_mid3(time,a+1,i))/Ntot_survey2(surv_flt_acomp); // placeholder note the indexing on ntot might be off
                  
                } // end else
              } // end ages
            } // end nspace
          }  // end flag
        } // end acomp survey fleets
        
        if(flag_catch(time) == 1){ // Flag if  there was a measurement that year
          
          for(int i=0;i<(nspace);i++){
            for(int a=0;a<(nage-1);a++){ // Loop over ages for catch comp
              if(a<age_maxage){
                age_catch_est(a,time) = (CatchNAge(a+1,time)/CatchN(time,i)); // Catch comp (1 bc the data starts at age = 1)
                age_catch_est2(time,a,i) = (CatchNAge2(time,a+1,i)/CatchN(time,i)); // Catch comp (1 bc the data starts at age = 1)
                
              }else{
                age_catch_est(age_maxage-1,time) += (CatchNAge(i+1,time)/CatchN(time,i));
                age_catch_est2(time,age_maxage-1,i) += (CatchNAge2(time,a+1,i)/CatchN(time,i));
              } // end else
            } // end ages
          } // end nspace
        } // end flag
  } // END TIME LOOP
  
  
  // OBSERVATION MODEL //
  // using namespace density;
  Type ans_survey=0.0;
  ////Save the observation model estimates
  for(int surv_flt =0;surv_flt<(nfleets_surv);surv_flt++){
    for(int time=1;time<tEnd;time++){ // Survey Surveyobs
      if(flag_surv_bio(time) == 2){
        // ans_survey += -dnorm(log(Surveyobs(time)), log(survey(time)), SDsurv+survey_err(time), TRUE);
        ans_survey += -dnorm(log(surv_pred(time,surv_flt)), log(survey2(time,surv_flt)), SDsurv+survey_err(time), TRUE); // the err also needs to be by flt
      } // end survey flag
      
    } // end time
  } // end surv_flt
  
  Type ans_catch = 0.0;
  for(int fish_flt =0;fish_flt<(nfleets_fish);fish_flt++){
    // for(int i=0;i<(nspace);i++){
    for(int time=0;time<tEnd;time++){ // Total Catches
      // ans_catch += -dnorm(log(Catch(time,i)+1e-6), log(Catchobs(time)+1e-6), SDcatch, TRUE); // this likelihood needs to be by fleet, not space
      ans_catch += -dnorm(log(Catch(time,fish_flt)+1e-6), log(Catchobs2(time,fish_flt)+1e-6), SDcatch, TRUE); // this likelihood needs to be by fleet, not space
      
    }
  }
  
  REPORT(ans_catch)
    
    // LIKELIHOOD: age comp in survey
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
          // sum1(time) += lgamma(ss_survey(time)*age_survey(a,time)+1);
          // sum2(time) += lgamma(ss_survey(time)*age_survey(a,time) + phi_survey*ss_survey(time)*age_survey_est(a,time)) -
          //   lgamma(phi_survey*ss_survey(time)*age_survey_est(a,time));
          // NOTE THAT THE age_survey2 OBS ARE IN A X TIME X FLEET, which is not the typical ordering
          sum1(time) += lgamma(ss_survey(time)*age_survey2(a,time,surv_flt_acomp)+1);
          sum2(time) += lgamma(ss_survey(time)*age_survey2(a,time,surv_flt_acomp) + phi_survey*ss_survey(time)*age_survey_est2(time,a,surv_flt_acomp)) -
            lgamma(phi_survey*ss_survey(time)*age_survey_est2(time,a,surv_flt_acomp));
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
  
  // for(int time=1;time<tEnd;time++){ // Loop over available years
  //   if(Catch(time)>0){
  //     
  //     if(flag_catch(time) == 1){ // Flag if  there was a measurement that year
  //       for(int i=0;i<age_maxage;i++){ // Loop over other ages (first one is empty for survey)
  //         sum3(time) += lgamma(ss_catch(time)*age_catch(i,time)+1);
  //         sum4(time) += lgamma(ss_catch(time)*age_catch(i,time) + phi_catch*ss_catch(time)*age_catch_est(i,time)) - lgamma(phi_catch*ss_catch(time)*age_catch_est(i,time));
  //       }
  //       ans_catchcomp += lgamma(ss_catch(time)+1)-sum3(time)+lgamma(phi_catch*ss_catch(time))-lgamma(ss_catch(time)+phi_catch*ss_catch(time))+sum4(time);
  //     }
  //   }
  // }
  
  
  Type ans_SDR = 0.0;
  
  for(int time=0;time<(tEnd-1);time++){ // Start time loop
    ans_SDR += Type(0.5)*(logR(time)*logR(time))/(SDR*SDR)+b(time)*log(SDR*SDR);
  }
  
  
  
  
  
  // Error for Selectivity
  Type ans_psel = 0.0;
  //
  for(int time=0;time<year_sel;time++){ // Start time loop
    for(int i=0;i<psel_fish.size();i++){ // Start time loop
      ans_psel += Type(0.5)*(PSEL(i,time)*PSEL(i,time))/(sigma_psel*sigma_psel);
    }
  }
  
  // Priors on h and M
  Type ans_priors = 0.0;
  
  for(int time=0;time<(nage-1);time++){ // Start time loop
    ans_priors += Type(0.5)*(initN(time)*initN(time))/(SDR*SDR);
  }
  
  // ans_priors += -dnorm(logh,log(Type(0.777)),Type(0.113),TRUE);
  
  // Prior on h
  ans_priors += -dbeta(h,Bprior,Aprior,TRUE);
  
  if(sum_zero == 1){
    ans_priors += ((Type(0.0)-sum(logR))*(Type(0.0)-sum(logR)))/Type(0.01);
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
  
  Type ans = ans_SDR+ans_psel+ans_catch+ans_survey-ans_survcomp-ans_catchcomp+ans_priors;
  //
  
  
  // Later Fix F in the likelihood and age comp in catch
  // Type ans = 0.0;
  // Report calculations
  ADREPORT(SSB)
    //ADREPORT(N)
    ADREPORT(Catch)
    ADREPORT(logF)
    ADREPORT(R)
    ADREPORT(Surveyobs)
    ADREPORT(Fyear)
    ADREPORT(surveyselc)
    ADREPORT(catchselec)
    ADREPORT(age_catch)
    ADREPORT(age_catch_est)
    ADREPORT(age_survey)
    ADREPORT(age_survey_est)
    ADREPORT(ans_tot)
    ADREPORT(SSBzero)
    
    // REPORT(N_beg2)
    REPORT(SSB)
    REPORT(SSB2)
    REPORT(Fyear)
    REPORT(Catch)
    REPORT(R)
    REPORT(R_k)
    REPORT(R_i)
    REPORT(Nzero2)
    REPORT(Nzero3)
    REPORT(CatchAge2)
    REPORT(CatchNAge2)
    REPORT(ans_tot)
    REPORT(Zsave)
    REPORT(age_survey_est)
    REPORT(age_catch_est)
    REPORT(age_survey_est2)
    REPORT(age_catch_est2)
    
    REPORT(CatchN)
    REPORT(selectivity_save)
    REPORT(surveyselc)
    REPORT(N_beg3)
    REPORT(N_beg2)
    REPORT(N_beg)
    REPORT(N_mid)
    REPORT(N_mid2)
    REPORT(surv_pred)
    REPORT(survey2)
    REPORT(N_mid3)
    REPORT(Surveyobs)
    
    return ans;
}
