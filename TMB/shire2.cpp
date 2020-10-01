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
  // // PARAMETERS //
  PARAMETER_VECTOR(logh_k); // Steepness by stock
  PARAMETER_VECTOR(logR_0k); // Recruitment at equil by stock
  PARAMETER_VECTOR(omega_0ij); // estimated age-0 movment among areas (used upon recruit gen)
  PARAMETER_VECTOR(logq_f); // Q by survey fleet
  PARAMETER_VECTOR(b); // bias adjustment factor
  PARAMETER(logSDR);
  PARAMETER_ARRAY(log_fsh_slx_pars);       // Fishery selectivity (selShape controls parameterization)
  PARAMETER_ARRAY(log_srv_slx_pars);       // Survey selectivity (selShape controls parameterization)
  
  Type ans = 1.0;
  return ans;
}
