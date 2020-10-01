#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  
  // DATA INPUT //
  // structure //
  DATA_INTEGER(nspace); // number of subreas to track
  DATA_INTEGER(nstocks); // number of stocks (demography)
  DATA_INTEGER(nage); // Plus group
  DATA_VECTOR(age); // ages
  DATA_INTEGER(tEnd); // number of years modeled
  DATA_VECTOR(years); // number of years modeled
  int nyear = years.size();
  // int nsex = 2;
  
  DATA_INTEGER(nfleets_surv); // number of survey fleets
  DATA_INTEGER(nfleets_fish); //number of fishery fleets
  DATA_INTEGER(nfleets_acomp); // number of age comp fleets
  DATA_INTEGER(nfleets_lcomp); //number of len comp fleets
  DATA_INTEGER(nmgmt_reg); // mgmt regions (3)
  
  DATA_IMATRIX(phi_if_surv); // turn on/off subareas for survey fleets
  DATA_IMATRIX(phi_if_fish); // turn on/off subareas for fishery fleets
  DATA_IMATRIX(phi_fm); //  fleets to mgmt areas
  DATA_IMATRIX(phi_fm_acomp); // acomp fleets to mgmt areas
  DATA_IMATRIX(phi_ki); // 0/1 nesting of subareas i into stocks k (rows)
  DATA_IMATRIX(phi_ik2); // vector stating which subarea (col) belongs to each stock k (value)
  DATA_MATRIX(tau_ki); // downscaling recruits from stocks to sub-areas
  DATA_IMATRIX(phi_fm_acomp2); //  fleets to mgmt areas
  // DATA_MATRIX(phi_lcomp_fm); //  fleets to mgmt areas
  
  // // DEMOGRAPHY // 
  DATA_VECTOR(mat_age); // natural mortality at age
  
  // // movement //
  DATA_ARRAY(omega_ais); // eigenvect of movement between subareas for ages > 0
  DATA_ARRAY(X_ijas); // prob trans between subareas at age
  // 
  // // growth //
  DATA_ARRAY(unfished_ALK_F); // for use in SSB0 calcs
  DATA_ARRAY(wtatlen_kab); // aL^b values by stock
  DATA_ARRAY(Linf_yk); // sex, stock, year specific
  DATA_ARRAY(L1_yk); // length at age 4 by stock; linear before this
  DATA_ARRAY(kappa_yk);
  DATA_ARRAY(sigmaG_yk); // perhaps turn to parameter later
  DATA_ARRAY(phi_ij); // matrix of whether i,j are from distinct stocks (0 otherwise)
  DATA_INTEGER(LBins); // maximum length bin in CM
  DATA_ARRAY(mla_yais); // for early or late periods, most likely length at age(for selex)
  DATA_ARRAY(mat_ak); // maturity at age for stock
  // 
  // // Selectivity
  DATA_IVECTOR(selType_fish); // 0 == AGESEL, 1= LENSEL
  DATA_IVECTOR(selShape_fish);
  DATA_IVECTOR(selType_surv); // 0 == AGESEL, 1 = LENSEL
  DATA_IVECTOR(selShape_surv);
  
  // // Time varying parameter blocks (indexed as h) - each vector contains the terminal years of
  // // each time block. Used for both selectivity and catchability
  DATA_IMATRIX(fsh_blks);        // fishery
  DATA_IMATRIX(srv_blks);       // survey
  // 
  // // Survey Biomass
  DATA_ARRAY(surv_yf_obs);
  // DATA_VECTOR(survey_err);
  array<Type> survey_yf_pred(nyear, nfleets_surv);
  
  // // Age Comps
  DATA_MATRIX(age_error); // nmgmt_reg x 100 ages
  DATA_MATRIX(age_error_SD); // nmgmt_reg x 100 ages
  DATA_IMATRIX(acomp_flt_type); // 0 for commercial, 1 for survey
  // // STORAGE ///
  // // Catches
  DATA_ARRAY(catch_yf_obs); // obs catch by year and fleet
  array<Type> catch_yaf_pred(tEnd, nage, nfleets_fish);  // estimated catches at age by fleet
  array<Type> catch_yf_pred(tEnd,nfleets_fish);
  array<Type> catch_yfi_pred(tEnd,nfleets_fish,nspace);
  array<Type> catch_yaif_pred(tEnd,nage,nspace,nfleets_fish);
  array<Type> CatchN_yaf(tEnd,nage,nfleets_fish);
  array<Type> N_avail_yf(tEnd, nfleets_fish);
  array<Type> N_weight_yfi(tEnd, nfleets_fish,nspace);
  // Switch for selectivity type: 0 = a50, a95 logistic; 1 = a50, slope logistic
  // Predicted selectivity
  array<Type> fsh_slx_yafs(nyear, LBins, nfleets_fish,2);           // Fishery selectivity-at-age by sex (on natural scale)
  array<Type> srv_slx_yafs(nyear, LBins, nfleets_surv+nfleets_acomp,2);  //
  // F tuning
  int niter = 50;
  array<Type> F1_yf(tEnd,nfleets_fish+1, niter+1); // intermediate f guess storage
  array<Type> F2_yf(tEnd,nfleets_fish+1, niter+1); // intermediate f guess storage
  array<Type> Freal_yf(tEnd,nfleets_fish); // final tuned fleet and yr specific F
  array<Type> Zreal_ya(tEnd,nage); // temp tuned fleet Z by y and age
  array<Type> Zreal_yai(tEnd,nage,nspace); // temp tuned fleet Z by y and age and area
  array<Type> F_area_yfi(tEnd,nfleets_fish,nspace); // temp tuned fleet Z by y and age
  array<Type> F_ym(tEnd,nmgmt_reg); //dodo
  array<Type> F_ydm(tEnd,nfleets_fish,nspace); //dodo
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
  // Recruits
  vector<Type> R(tEnd);
  array<Type>  R_yk(tEnd,nstocks); // stock-level recruitment (bev-holt)
  array<Type>  R_yi(tEnd,nspace); // subarea-level recruitment (downscaled)
  // Length at age
  array<Type> Length_yais_beg(tEnd+1,nage,nspace,2); // placeholder for true lengths-at-age
  array<Type> Length_yais_mid(tEnd+1,nage,nspace,2); // placeholder for true lengths-at-age
  array<Type> Length_yais_end(tEnd+1,nage,nspace,2); // placeholder for true lengths-at-age
  array<Type> LengthAge_alyis_beg(nage,LBins,tEnd+1,nspace,2); // placeholder for true age-length dist
  array<Type> LengthAge_alyis_mid(nage,LBins,tEnd+1,nspace,2); // placeholder for true age-length dist
  array<Type> LengthAge_alyis_end(nage,LBins,tEnd+1,nspace,2); // placeholder for true age-length dist
  // age comps
  array<Type> acomp_yaf_temp(tEnd, nage, nfleets_acomp); // predicted acomps from commercial fisheries
  array<Type> comm_acomp_yafs_pred(tEnd, nage, 2, 2); // predicted acomps from commercial fisheries
  array<Type> surv_acomp_yafs_pred(tEnd, nage, 6, 2); // predicted acomps from surveys (without biomass)
  array<Type> Nsamp_acomp_yf(tEnd, nfleets_acomp); // placeholder for number sampled by comp survey (pre dirichlet weighting)
  // // PARAMETERS //
  PARAMETER_VECTOR(logh_k); // Steepness by stock
  PARAMETER_VECTOR(logR_0k); // Recruitment at equil by stock
  PARAMETER_VECTOR(omega_0ij); // estimated age-0 movment among areas (used upon recruit gen)
  PARAMETER_VECTOR(logq_f); // Q by survey fleet
  PARAMETER_VECTOR(b); // bias adjustment factor
  PARAMETER(logSDR);
  PARAMETER_ARRAY(log_fsh_slx_pars);       // Fishery selectivity (selShape controls parameterization)
  PARAMETER_ARRAY(log_srv_slx_pars);       // Survey selectivity (selShape controls parameterization)
  // // Transform out of log space
  // // Type SDsurv = exp(logSDsurv);
  // // Type SDcatch = exp(logSDcatch);
  Type SDR = exp(logSDR);
  vector<Type> R_0k = exp(logR_0k);
  vector<Type> h_k = exp(logh_k);
  vector<Type> q_f = exp(logq_f);
  array<Type> tildeR_yk(tEnd,nstocks); // recdevs
  vector<Type> tildeR_initk(nstocks); // recdevs for init

  // Fishery selectivity
  // Number of parameters in the chosen selectivity type: 
  int npar_slx = log_fsh_slx_pars.dim(1); // dim = array dimensions; 1 = # columns in array = # params in slx_type
  // Preliminary calcs to bring parameters out of log space
  vector<int> a1_dim =log_fsh_slx_pars.dim;
  array<Type> fsh_slx_pars(a1_dim);
  fsh_slx_pars.setZero();
  for (int fish_flt = 0; fish_flt < nfleets_fish; fish_flt++) {
    for (int n = 0; n < npar_slx; n++) { // loop over alpha and beta
      // for (int h = 0; h < fsh_blks.size(); h++) { // loop time blocks
        for (int s = 0; s < 2; s++) { // loop sexes
          fsh_slx_pars(fish_flt,n,0,s) = exp(log_fsh_slx_pars(fish_flt,n,0,s));
        } // end sex
      // } // end blocks
    } // end alpha, beta
  } // end fish fleets
  // Notes on the following syntax: the do while allows you to estimate parameters within a y block. It
  // "does" the looping over year and age "while" within the y block, then
  // iterates to the next block. Year is not in a for loop because it is
  // iterated by the do statement.

  // The switch for slx_shape allows you to change parameterization SHAPE. This could
  // easily be expanded to accomodate any selectivity type (the fsh_slx_pars
  // allows for a flexible number of parameters and y blocks)
  // slx pars setup is fleet x alpha, beta x time block (1 for now) x sex 
  for(int fish_flt =0;fish_flt<(nfleets_fish);fish_flt++){ // loop fleets
    int i = 0;
    for(int y = 0; y < nyear; y++){ // loop years
      do{
        
        //       switch (selType_fish(fish_flt)) { // 0 is age, 1 is leng
        //       case 0: // enter age based sel
        //         for (int s = 0; s < 2; s++) { // loop sexes
        //           // Selectivity switch (case 0 or 1 references the value of slx_type)
        //           switch (selShape_fish(fish_flt)) { // age sel
        //           case 0: // Logistic with a50 and a95, where  fsh_slx_pars(fish_flt,0,0,s) = a50 and  fsh_slx_pars(fish_flt,1,0,s) = a95
        //             for (int a= 0; a < nage; a++){
        //               fsh_slx_yafs(i,a,fish_flt,s) = Type(1.0) / ( Type(1.0) + exp(-log(Type(19)) *
        //                 (a -  fsh_slx_pars(fish_flt,0,0,s)) / ( fsh_slx_pars(fish_flt,1,0,s) - 
        //                 fsh_slx_pars(fish_flt,0,0,s))));
        //             } // end ages
        //             break;
        //           // case 1: // Logistic with a50 and slope, where  fsh_slx_pars(fish_flt,0,0,s) = a50 and  fsh_slx_pars(fish_flt,1,0,s) = slope.
        //           //   //  *This is the preferred logistic parameterization b/c it reduces parameter correlation*
        //           //   for (int a= 0; a < nage; a++){
        //           //     fsh_slx_yafs(i,a,fish_flt,s)  = Type(1.0) / ( Type(1.0) + exp( Type(-1.0) *
        //           //        fsh_slx_pars(fish_flt,1,0,s) * (a -  fsh_slx_pars(fish_flt,0,0,s)) ) );
        //           //   } // end ages
        //           //   break;
        //           // case 2: // Dome Normal with alpha (mean) and beta (sd)
        //           //   for (int a= 0; a < nage; a++){
        //           //     fsh_slx_yafs(i,a,fish_flt,s)  = exp(-(0.5 * (a - fsh_slx_pars(y,2,s))/pow(fsh_slx_pars(y,3,s),2)));
        //           //   } // end ages
        //           //   break;
        //           // case 3: // Dome Gamma with alpha (mean) and beta (sd)
        //           //   vector<Type>selG(nage);
        //           //   for (int a= 0; a < nage; a++) {
        //           //     selG(a)= pow(a, (fsh_slx_pars(y,2,s) - 1)) * exp(-a/fsh_slx_pars(y,3,s));
        //           //   } // end ages
        //           //   for (int a= 0;a < nage; a++) {
        //           //     fsh_slx_yafs(i,a,fish_flt,s) = selG(a) / max(selG);
        //           //   } // end ages
        //             break;
        //           } // end sex
        //         } // end switch selShape
        //         break; // break age sel
        //       // case 1: // enter length based sel
        //       //   for (int s = 0; s < 2; s++) {
        //       //     switch (selShape_fish(fish_flt)) {
        //       //     case 0: // Logistic with a50 and a95, where  fsh_slx_pars(fish_flt,0,0,s) = a50 and  fsh_slx_pars(fish_flt,1,0,s) = a95
        //       //       for (int l = 0; l < LBins; l++){
        //       //         fsh_slx_yafs(i,l,fish_flt,s) = Type(1.0) / ( Type(1.0) + exp(-log(Type(19)) *
        //       //           (l -  fsh_slx_pars(fish_flt,0,0,s)) / ( fsh_slx_pars(fish_flt,1,0,s) -  fsh_slx_pars(fish_flt,0,0,s))));
        //       //       } // end ages
        //       //       break;
        //       //     case 1: // Logistic with a50 and slope, where  fsh_slx_pars(fish_flt,0,0,s) = a50 and  fsh_slx_pars(fish_flt,1,0,s) = slope.
        //       //       //  *This is the preferred logistic parameterization b/c it reduces parameter correlation*
        //       //       for (int l = 0; l < LBins; l++){
        //       //         fsh_slx_yafs(i,l,fish_flt,s)  = Type(1.0) / ( Type(1.0) + exp( Type(-1.0) *
        //       //            fsh_slx_pars(fish_flt,1,0,s) * (l -  fsh_slx_pars(fish_flt,0,0,s)) ) );
        //       //       } // end len
        //       //       break;
        //       //     case 2: // Dome Normal with alpha (mean) and beta (sd)
        //       //       for (int l = 0; l < LBins; l++){
        //       //         fsh_slx_yafs(i,l,fish_flt,s)  = exp(-(0.5 * (l - fsh_slx_pars(y,2,s))/pow(fsh_slx_pars(y,3,s),2)));
        //       //       } // end len
        //       //       break;
        //       //     case 3: // Dome Gamma with alpha (mean) and beta (sd)
        //       //       vector<Type>selG(LBins);
        //       //       for (int l = 0; l < LBins; l++){
        //       //         selG(l)= pow(l, (fsh_slx_pars(y,2,s) - 1)) * exp(-l/fsh_slx_pars(y,3,s));
        //       //       } // end len
        //       //       for (int l = 0; l < LBins; l++){
        //       //         fsh_slx_yafs(i,l,fish_flt,s) = selG(l) / max(selG);
        //       //       } // end len
        //       //       break;
        //       //     } // end switch selShape
        //       //   } // end sex
        //         break;
        //       } // end switch selType
        // 
        i++;
      } while (i <= fsh_blks(y,fish_flt));
    } // end y blocks
  } // end fish_flt
  
  Type ans = 1.0;
  return ans;
}
