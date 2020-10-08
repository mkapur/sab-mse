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
  int nsex = 2;

  DATA_INTEGER(nfleets_surv); // number of survey fleets
  DATA_INTEGER(nfleets_fish); //number of fishery fleets
  DATA_INTEGER(nfleets_acomp); // number of age comp fleets
  // DATA_INTEGER(nfleets_lcomp); //number of len comp fleets
  DATA_INTEGER(nmgmt_reg); // mgmt regions (3)

  DATA_IMATRIX(phi_if_surv); // turn on/off subareas for survey fleets
  DATA_IMATRIX(phi_if_fish); // turn on/off subareas for fishery fleets
  DATA_IMATRIX(phi_if_acomp); // turn on/off subareas for fishery fleets
  DATA_IMATRIX(phi_fm); //  fleets to mgmt areas
  DATA_IMATRIX(phi_ff_acomp); // where the acomp fleets are in fish/surv fleet vectors
  DATA_IMATRIX(phi_fm_acomp); // acomp fleets to mgmt areas
  DATA_IMATRIX(phi_im); // 0/1 nesting of subareas i into stocks k (rows)
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
  DATA_IARRAY(mla_yais); // most likely length bin in age bin
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
  DATA_ARRAY(surv_yf_err);
  array<Type> surv_yf_pred(nyear, nfleets_surv);

  // Age Comps
  DATA_MATRIX(age_error); // nmgmt_reg x 100 ages
  DATA_MATRIX(age_error_SD); // nmgmt_reg x 100 ages
  DATA_IMATRIX(acomp_flt_type); // 0 for commercial, 1 for survey
  DATA_ARRAY(acomp_yafs_obs);
  // // STORAGE ///
  // // Catches
  DATA_ARRAY(catch_yf_obs); // obs catch by year and fleet
  DATA_ARRAY(catch_yf_error); // right now just 0.1 for all fleets and all years
  array<Type> catch_yaf_pred(tEnd, nage, nfleets_fish);  // estimated catches at age by fleet
  array<Type> catch_yf_pred(tEnd,nfleets_fish);
  array<Type> catch_yfi_pred(tEnd,nfleets_fish,nspace);
  array<Type> catch_yaif_pred(tEnd,nage,nspace,nfleets_fish);
  array<Type> CatchN_yaf(tEnd,nage,nfleets_fish);
  array<Type> N_avail_yf(tEnd, nfleets_fish);
  array<Type> N_weight_yfi(tEnd, nfleets_fish,nspace);
  // Switch for selectivity type: 0 = a50, a95 logistic; 1 = a50, slope logistic
  // Predicted selectivity
  array<Type> fsh_slx_yafs(nyear, LBins, nfleets_fish, nsex);           // Fishery selectivity-at-age by sex (on natural scale)
  array<Type> srv_slx_yafs(nyear, LBins, nfleets_surv+(nfleets_acomp-5),nsex);  // five of the acomp fleets are surveys; the other two are fsh
  vector<Type>selG(nage);
  vector<Type>selGL(LBins);

  // F tuning
  int niter = 5;
  vector<Type>Z_a_TEMP(nage);
  vector<Type>Z_a_TEMP2(nage);
  array<Type> catch_afk_TEMP(nage, nfleets_fish, niter+1);
  array<Type> F1_yf(tEnd,nfleets_fish, niter+1); // intermediate f guess storage
  array<Type> F2_yf(tEnd,nfleets_fish, niter+1); // intermediate f guess storage
  array<Type> Freal_yf(tEnd,nfleets_fish); // final tuned fleet and yr specific F
  array<Type> Zreal_ya(tEnd,nage); // temp tuned fleet Z by y and age
  array<Type> Zreal_yai(tEnd,nage,nspace); // temp tuned fleet Z by y and age and area
  array<Type> F_area_yfi(tEnd,nfleets_fish,nspace); // temp tuned fleet Z by y and age
  array<Type> F_ym(tEnd,nmgmt_reg); //dodo
  array<Type> F_ydm(tEnd,nfleets_fish,nspace); //dodo
  // biology storage
  array<Type> Ninit_ais(nage,nspace,nsex); // initial numbers at age in subarea, just once
  array<Type> N_0ais(nage, nspace,nsex); // numbers in year 0 at age in subarea
  vector<Type> SSB_0k(nstocks); // virgin spawnbio by stock
  vector<Type> SSB_0i(nspace); // virgin spawnbio by subarea
  array<Type> N_yais_beg( tEnd+1, nage, nspace,nsex); N_yais_beg.setZero();
  array<Type> N_yais_mid( tEnd+1, nage, nspace,nsex); N_yais_mid.setZero();
  array<Type> N_yais_end( tEnd+1, nage, nspace,nsex); N_yais_end.setZero();
  array<Type> SSB_yk(tEnd,nstocks);
  array<Type> SSB_yi(tEnd,nspace);
  array<Type> SSB_ym(tEnd,nmgmt_reg);
  // Recruits
  array<Type>  R_yk(tEnd,nstocks); // stock-level recruitment (bev-holt)
  array<Type>  R_yi(tEnd,nspace); // subarea-level recruitment (downscaled)
  array<Type>  R_ym(tEnd,nmgmt_reg); // subarea-level recruitment (downscaled)

  // Length at age
  array<Type> Length_yais_beg(tEnd+1,nage,nspace,nsex); // placeholder for true lengths-at-age
  array<Type> Length_yais_mid(tEnd+1,nage,nspace,nsex); // placeholder for true lengths-at-age
  array<Type> Length_yais_end(tEnd+1,nage,nspace,nsex); // placeholder for true lengths-at-age
  array<Type> LengthAge_alyis_beg(nage,LBins,tEnd+1,nspace,nsex); // placeholder for true age-length dist
  array<Type> LengthAge_alyis_mid(nage,LBins,tEnd+1,nspace,nsex); // placeholder for true age-length dist
  array<Type> LengthAge_alyis_end(nage,LBins,tEnd+1,nspace,nsex); // placeholder for true age-length dist
  // age comps
  array<Type> acomp_yaf_temp(tEnd, nage, nfleets_acomp); // placeholder multiplier for all acomp fleets
  array<Type> comm_acomp_yafs_pred(tEnd, nage, 5, nsex); // predicted acomps from commercial fisheries
  array<Type> surv_acomp_yafs_pred(tEnd, nage, nfleets_acomp-5, nsex); // predicted acomps from surveys (without biomass)
  array<Type> Nsamp_acomp_yf(tEnd, nfleets_surv+nfleets_acomp); // placeholder for number sampled by comp survey (pre dirichlet weighting)
  // // PARAMETERS //
  PARAMETER_VECTOR(logh_k); // Steepness by stock
  PARAMETER_VECTOR(logR_0k); // Recruitment at equil by stock
  PARAMETER_VECTOR(omega_0ij); // estimated age-0 movment among areas (used upon recruit gen)
  PARAMETER_VECTOR(logq_f); // Q by survey fleet
  PARAMETER_VECTOR(b); // bias adjustment factor
  PARAMETER_VECTOR(logpi_acomp); // dirichlet scalar for acomp sampling
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
  vector<Type> pi_acomp = exp(logpi_acomp);
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
      for (int s = 0; s < nsex; s++) { // loop sexes
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
    for(int y = 0; y < nyear; y++){ // loop years; this should really loop over the # of blocks and replace the fixed zero
      do{
        switch (selType_fish(fish_flt)) { // 0 is age, 1 is leng
        case 0: // enter age based sel
          for (int s = 0; s < nsex; s++) { // loop sexes
            // Selectivity switch (case 0 or 1 references the value of slx_type)
            switch (selShape_fish(fish_flt)) { // age sel
            case 0: // Logistic with a50 and a95, where  fsh_slx_pars(fish_flt,0,0,s) = a50 and  fsh_slx_pars(fish_flt,1,0,s) = a95
              for (int a= 0; a < nage; a++){
                fsh_slx_yafs(i,a,fish_flt,s) = Type(1.0) / ( Type(1.0) + exp(-log(Type(19)) *
                  (a -  fsh_slx_pars(fish_flt,0,0,s)) / ( fsh_slx_pars(fish_flt,1,0,s) -
                  fsh_slx_pars(fish_flt,0,0,s))));
              } // end ages
              break;
            case 1: // Logistic with a50 and slope, where  fsh_slx_pars(fish_flt,0,0,s) = a50 and  fsh_slx_pars(fish_flt,1,0,s) = slope.
              //  *This is the preferred logistic parameterization b/c it reduces parameter correlation*
              for (int a= 0; a < nage; a++){
                fsh_slx_yafs(i,a,fish_flt,s)  = Type(1.0) / ( Type(1.0) + exp( Type(-1.0) *
                  fsh_slx_pars(fish_flt,1,0,s) * (a -  fsh_slx_pars(fish_flt,0,0,s)) ) );
              } // end ages
              break;
            case 2: // Dome Normal with alpha (mean) and beta (sd)
              for (int a= 0; a < nage; a++){
                fsh_slx_yafs(i,a,fish_flt,s)  = exp(-(0.5 * (a -    fsh_slx_pars(fish_flt,0,0,s))/pow(   fsh_slx_pars(fish_flt,1,0,s),2)));
              } // end ages
              break;
            case 3: // Dome Gamma with alpha (mean) and beta (sd)
              selG.setZero();
              for (int a= 0; a < nage; a++) {
                selG(a)= pow(a, (   fsh_slx_pars(fish_flt,0,0,s) - 1)) * exp(-a/   fsh_slx_pars(fish_flt,1,0,s));
              } // end ages
              for (int a= 0;a < nage; a++) {
                fsh_slx_yafs(i,a,fish_flt,s) = selG(a) / max(selG);
              } // end ages
              break;
            } // end switch selShape
          } // end sex
          break; // break age sel
        case 1: // enter length based sel
          for (int s = 0; s < nsex; s++) {
            switch (selShape_fish(fish_flt)) {
            case 0: // Logistic with a50 and a95, where  fsh_slx_pars(fish_flt,0,0,s) = a50 and  fsh_slx_pars(fish_flt,1,0,s) = a95
              for (int l = 0; l < LBins; l++){
                fsh_slx_yafs(i,l,fish_flt,s) = Type(1.0) / ( Type(1.0) + exp(-log(Type(19)) *
                  (l -  fsh_slx_pars(fish_flt,0,0,s)) / ( fsh_slx_pars(fish_flt,1,0,s) -  fsh_slx_pars(fish_flt,0,0,s))));
              } // end ages
              break;
            case 1: // Logistic with a50 and slope, where  fsh_slx_pars(fish_flt,0,0,s) = a50 and  fsh_slx_pars(fish_flt,1,0,s) = slope.
              //  *This is the preferred logistic parameterization b/c it reduces parameter correlation*
              for (int l = 0; l < LBins; l++){
                fsh_slx_yafs(i,l,fish_flt,s)  = Type(1.0) / ( Type(1.0) + exp( Type(-1.0) *
                  fsh_slx_pars(fish_flt,1,0,s) * (l -  fsh_slx_pars(fish_flt,0,0,s)) ) );
              } // end len
              break;
            case 2: // Dome Normal with alpha (mean) and beta (sd)
              for (int l = 0; l < LBins; l++){
                fsh_slx_yafs(i,l,fish_flt,s)  = exp(-(0.5 * (l -    fsh_slx_pars(fish_flt,0,0,s))/pow(   fsh_slx_pars(fish_flt,1,0,s),2)));
              } // end len
              break;
            case 3: // Dome Gamma with alpha (mean) and beta (sd)
              selGL.setZero();
              for (int l = 0; l < LBins; l++){
                selGL(l)= pow(l, (   fsh_slx_pars(fish_flt,0,0,s) - 1)) * exp(-l/   fsh_slx_pars(fish_flt,1,0,s));
              } // end len
              for (int l = 0; l < LBins; l++){
                fsh_slx_yafs(i,l,fish_flt,s) = selGL(l) / max(selGL);
              } // end len
              break;
            } // end switch selShape
          } // end sex for case 1
          break;
        } // end switch selType
        i++;
      } while (i <= fsh_blks(y,fish_flt)); // bracket i estimation for years designated by this block
    } // end y blocks
  } // end fish_flt

  vector<int> a2_dim = log_srv_slx_pars.dim;
  array<Type> srv_slx_pars(a2_dim);
  srv_slx_pars.setZero();
  for (int srv_flt = 0; srv_flt < a2_dim(0); srv_flt++) {
    for (int n = 0; n < npar_slx; n++) { // loop over alpha and beta
      // for (int h = 0; h < srv_blks.size(); h++) { // loop time blocks
      for (int s = 0; s < nsex; s++) { // loop sexes
        srv_slx_pars(srv_flt,n,0,s) = exp(log_srv_slx_pars(srv_flt,n,0,s));
      } // end sex
      // } // end blocks
    } // end alpha, beta
  } // end srv fleets
  // doing five of these to account for five surveys w acomp
  for(int srv_flt =0;srv_flt<(nfleets_surv+(nfleets_acomp-5));srv_flt++){ // loop fleets
    int i = 0; // re-set i to 0
    for(int y = 0; y < nyear; y++){ // loop years; this should really loop over the # of blocks and replace the fixed zero
      do{
        switch (selType_surv(srv_flt)) { // 0 is age, 1 is leng
        case 0: // enter age based sel
          for (int s = 0; s < nsex; s++) { // loop sexes
            // Selectivity switch (case 0 or 1 references the value of slx_type)
            switch (selShape_surv(srv_flt)) { // age sel
            case 0: // Logistic with a50 and a95, where  srv_slx_pars(srv_flt,0,0,s) = a50 and  srv_slx_pars(srv_flt,1,0,s) = a95
              for (int a= 0; a < nage; a++){
                srv_slx_yafs(i,a,srv_flt,s) = Type(1.0) / ( Type(1.0) + exp(-log(Type(19)) *
                  (a -  srv_slx_pars(srv_flt,0,0,s)) / ( srv_slx_pars(srv_flt,1,0,s) -
                  srv_slx_pars(srv_flt,0,0,s))));
              } // end ages
              break;
            case 1: // Logistic with a50 and slope, where  srv_slx_pars(srv_flt,0,0,s) = a50 and  srv_slx_pars(srv_flt,1,0,s) = slope.
              //  *This is the preferred logistic parameterization b/c it reduces parameter correlation*
              for (int a= 0; a < nage; a++){
                srv_slx_yafs(i,a,srv_flt,s)  = Type(1.0) / ( Type(1.0) + exp( Type(-1.0) *
                  srv_slx_pars(srv_flt,1,0,s) * (a -  srv_slx_pars(srv_flt,0,0,s)) ) );
              } // end ages
              break;
            case 2: // Dome Normal with alpha (mean) and beta (sd)
              for (int a= 0; a < nage; a++){
                srv_slx_yafs(i,a,srv_flt,s)  = exp(-(0.5 * (a -    srv_slx_pars(srv_flt,0,0,s))/pow(   srv_slx_pars(srv_flt,1,0,s),2)));
              } // end ages
              break;
            case 3: // Dome Gamma with alpha (mean) and beta (sd)
              selG.setZero();
              for (int a= 0; a < nage; a++) {
                selG(a)= pow(a, (   srv_slx_pars(srv_flt,0,0,s) - 1)) * exp(-a/   srv_slx_pars(srv_flt,1,0,s));
              } // end ages
              for (int a= 0;a < nage; a++) {
                srv_slx_yafs(i,a,srv_flt,s) = selG(a) / max(selG);
              } // end ages
              break;
            } // end switch selShape
          } // end sex
          break; // break age sel
        case 1: // enter length based sel
          for (int s = 0; s < nsex; s++) {
            switch (selShape_surv(srv_flt)) {
            case 0: // Logistic with a50 and a95, where  srv_slx_pars(srv_flt,0,0,s) = a50 and  srv_slx_pars(srv_flt,1,0,s) = a95
              for (int l = 0; l < LBins; l++){
                srv_slx_yafs(i,l,srv_flt,s) = Type(1.0) / ( Type(1.0) + exp(-log(Type(19)) *
                  (l -  srv_slx_pars(srv_flt,0,0,s)) / ( srv_slx_pars(srv_flt,1,0,s) -  srv_slx_pars(srv_flt,0,0,s))));
              } // end ages
              break;
            case 1: // Logistic with a50 and slope, where  srv_slx_pars(srv_flt,0,0,s) = a50 and  srv_slx_pars(srv_flt,1,0,s) = slope.
              //  *This is the preferred logistic parameterization b/c it reduces parameter correlation*
              for (int l = 0; l < LBins; l++){
                srv_slx_yafs(i,l,srv_flt,s)  = Type(1.0) / ( Type(1.0) + exp( Type(-1.0) *
                  srv_slx_pars(srv_flt,1,0,s) * (l -  srv_slx_pars(srv_flt,0,0,s)) ) );
              } // end len
              break;
            case 2: // Dome Normal with alpha (mean) and beta (sd)
              for (int l = 0; l < LBins; l++){
                srv_slx_yafs(i,l,srv_flt,s)  = exp(-(0.5 * (l -    srv_slx_pars(srv_flt,0,0,s))/pow(   srv_slx_pars(srv_flt,1,0,s),2)));
              } // end len
              break;
            case 3: // Dome Gamma with alpha (mean) and beta (sd)
              selGL.setZero();
              for (int l = 0; l < LBins; l++){
                selG(l)= pow(l, (   srv_slx_pars(srv_flt,0,0,s) - 1)) * exp(-l/   srv_slx_pars(srv_flt,1,0,s));
              } // end len
              for (int l = 0; l < LBins; l++){
                srv_slx_yafs(i,l,srv_flt,s) = selG(l) / max(selG);
              } // end len
              break;
            } // end switch selShape
          } // end sex for case 1
          break;
        } // end switch selType
        i++;
      } while (i <= srv_blks(y,srv_flt)); // bracket i estimation for years designated by this block
    } // end y blocks
  } // end srv
  REPORT(fsh_slx_yafs);
  REPORT(srv_slx_yafs);
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
  for(int s=0;s<nsex;s++){
    for(int k=0;k<(nstocks);k++){
      for(int i=0;i<(nspace);i++){
        for(int a=0;a<(nage-1);a++){
          N_0ais(a,i,s) = 0.5*R_0k(k)*tau_ki(k,i)*exp(-(mat_age(a)*age(a))); // compound multiply duh
        }  // note the A+ group will be in slot A-1
        N_0ais(nage-1,i,s) =  N_0ais(nage-2,i,s)*exp(-sum(mat_age))
          /(Type(1.0)-exp(-mat_age(nage-1)));
      } // end subareas
    } // end stocks
  } // end sex

  // Equilibrium Unfished SSB, stock (outside of y loop)
  for(int i=0;i<(nspace);i++){
    for(int a=0;a<nage;a++){ // Loop over ages
      SSB_0i(i) += mat_age(a)*
        N_0ais(a,i,0)*
        wtatlen_kab(phi_ik2(i),0)*
        pow(unfished_ALK_F(a,i),wtatlen_kab(phi_ik2(i),1))*
        mat_ak(a,phi_ik2(i));
      for(int k=0;k<(nstocks);k++){
        SSB_0k(k) += phi_ki(k,i)*SSB_0i(i);
      } // end stocks
    } // end ages
  } // end space
  //
  // // The first year of the simulation is initialized with the following age distribution
  Ninit_ais.setZero();
  for(int y=0;y<(10*nage);y++){
    for(int s=0;s<nsex;s++){
      for(int i=0;i<(nspace);i++){
        for(int a=0;a<(nage-1);a++){
          Ninit_ais(a,i,s) += 0.5* omega_ais(a,i,s) *
            tau_ki(phi_ik2(i),i) *
            R_0k(phi_ik2(i))* exp(-(mat_age(a)*age(a))) *
            exp(-0.5*SDR*SDR+tildeR_initk(phi_ik2(i)));
        } // end ages
        Ninit_ais(nage-1,i,s) += (omega_ais(nage-1,i,s) * Ninit_ais(nage-2,i,s) *
          exp(-mat_age(nage-1)*age(nage-1)))/(Type(1.0)-exp(-sum(mat_age))*
          exp(-0.5*SDR*SDR+tildeR_initk(phi_ik2(i))));
      } // end space
    } // end sex
  } // end yinit
  // std::cout << "Done" << std::endl;

  // std::cout << " Here" << "\n";
  // for(int y=0;y<(tEnd);y++){ // Start y loop
    for(int y=0;y<3;y++){ // Start y loop
      
    // model year zero, use last year of Ninit_ai, and equil movement (omega) and downscaling (tau)
    // note we are assuming unfished here as the exponent is M only
    // note that in tmb, plus group is in slot nage-1
    // so the incoming plus-groupers will be in slots nage-1 or nage-2 in prior year
    std::cout << y << " start year loop" << "\n";
    if (y == 0){
      for(int i=0;i<(nspace);i++){
        for(int s=0;s<nsex;s++){
          Length_yais_beg(0,0,i,s) = L1_yk(y,phi_ik2(i),s);
          N_yais_beg(0,0,i,s) = Ninit_ais(0,i,s);
          N_yais_mid(0,0,i,s) = N_yais_beg(0,0,i,s)*exp(-mat_age(0)/2);
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
              NCome += X_ijas(j,i,nage-1,s)*(Ninit_ais(nage-1,j,s) + Ninit_ais(nage-2,j,s)) ;
            } // end i != j
          } // end subareas j
          N_yais_beg(y,nage-1,i,s) =  ((1-pLeave)*(Ninit_ais(nage-1,i,s)+Ninit_ais(nage-2,i,s)) + NCome)*
            exp(-mat_age(nage-1)/2);
          Length_yais_beg(y,nage-1,i,s) = Linf_yk(0,phi_ik2(i),s)+(L1_yk(0,phi_ik2(i),s)-Linf_yk(0,phi_ik2(i),s))*
            exp(-kappa_yk(0,phi_ik2(i),s)*(nage-1));
          Length_yais_mid(y,nage-1,i,s) =  Linf_yk(0,phi_ik2(i),s)+(L1_yk(0,phi_ik2(i),s)-Linf_yk(0,phi_ik2(i),s))*
            exp(-0.5*kappa_yk(0,phi_ik2(i),s)*(nage-1));
        } // end sexes
      } // end subareas i
    } // end y == 0
    std::cout << y << " did year zero" << "\n";
    
    Type lenstep = 0.0; Type lenslope = 0.0;
    // N- and Nominal Length - at-age for the middle of this year 
    for(int s=0;s<nsex;s++){
      for(int i=0;i<(nspace);i++){
        N_yais_mid(y,0,i,s) = N_yais_beg(y,0,i,s)*exp(-mat_age(0)/2);
        // linear growth below A4 as in synthesis
        if(L1_yk(y,phi_ik2(i),s) < 3){
          lenstep = L1_yk(y,phi_ik2(i),s);
          lenslope = (L1_yk(y,phi_ik2(i),s) - lenstep) / 3;
        } else if(L1_yk(y,phi_ik2(i),s) >= 3){
          lenstep = 3.0;
          lenslope = (L1_yk(y,phi_ik2(i),s) - lenstep) / 3;
        }
        for(int a=0;a<4;a++){
          Length_yais_beg(y,a,i,s) = lenstep+lenslope*a;
        } // end linear age
        Length_yais_beg(y,4,i,s) =  L1_yk(y,phi_ik2(i),s);
        // for(int a=0;a<4;a++){
          // Length_yais_mid(y,a,i,s) = Length_yais_beg(y,a,i,s) + (Linf_yk(y,phi_ik2(i),s)-Length_yais_beg(y,a,i,s)*
          //   (1-exp(-0.5*kappa_yk(y,phi_ik2(i),s))));
        // } // end linear age
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
          Length_yais_beg(y,a,i,s) =  Linf_yk(y,phi_ik2(i),s)+(L1_yk(y,phi_ik2(i),s)-Linf_yk(y,phi_ik2(i),s))*
            exp(-kappa_yk(y,phi_ik2(i),s)*a);
          // Length_yais_beg(y+1,a,i,s)  = Length_yais_beg(y,a-1,i,s) +
          //   (Linf_yk(y,phi_ik2(i),s)-Length_yais_beg(y,a-1,i,s))*
          //   (1-exp(-kappa_yk(y,phi_ik2(i),s)));
          Length_yais_mid(y,a,i,s)= Length_yais_beg(y,a,i,s) +
            (Linf_yk(y,phi_ik2(i),s)-Length_yais_beg(y,a,i,s))*
            (1-exp(-0.5*kappa_yk(y,phi_ik2(i),s)));
        } // end nonlinear growth ages
        // plus groups
        Type pLeave = 0.0; Type NCome = 0.0;
        for(int j=0;j<(nspace);j++){
          if(i != j){
            pLeave += X_ijas(i,j,nage-1,s);
            NCome += X_ijas(j,i,nage-1,s)*(N_yais_beg(y,nage-1,j,s) + N_yais_beg(y,nage-2,j,s));
          } // end i != j
        } // end subareas j
        N_yais_mid(y,nage-1,i,s) =((1-pLeave)*N_yais_beg(y,nage-1,i,s) + NCome)*exp(-mat_age(nage-1)/2);
        // plus group weighted average (we already have the numbers at age)
        Length_yais_beg(y,nage-1,i,s) = (N_yais_beg(y,nage-2,i,s)*
          (Length_yais_beg(y,nage-2,i,s)+
          (Linf_yk(y,phi_ik2(i),s))-
          Length_yais_beg(y,nage-2,i,s)*
          (1-exp(-kappa_yk(y,phi_ik2(i),s)))) +
          N_yais_beg(y,nage-1,i,s)  *
          (Length_yais_beg(y,nage-1,i,s) +
          (Linf_yk(y,phi_ik2(i),s) -
          Length_yais_beg(y,nage-1,i,s))*(1-exp(-kappa_yk(y,phi_ik2(i),s)))))/
            (N_yais_beg(y,nage-2,i,s) + N_yais_beg(y,nage-1,i,s));
        Length_yais_mid(y,nage-1,i,s) = (N_yais_mid(y,nage-2,i,s)*
          (Length_yais_beg(y,nage-2,i,s)+(Linf_yk(y,phi_ik2(i),s)-
          Length_yais_beg(y,nage-2,i,s)*(1-exp(-0.5*kappa_yk(y,phi_ik2(i),s))))) +
          N_yais_mid(y,nage-1,i,s)*
          (Length_yais_beg(y,nage-1,i,s)+(Linf_yk(y,phi_ik2(i),s)-
          Length_yais_beg(y,nage-1,i,s))*(1-exp(-0.5*kappa_yk(y,phi_ik2(i),s)))))/
            (N_yais_mid(y,nage-2,i,s) + N_yais_mid(y,nage-1,i,s));
      } // end subareas i
    } // end sexes
    std::cout << y << " before prob LAA" << "\n";
    // prob of length-at-age
    for(int s=0;s<nsex;s++){
      for(int i=0;i<(nspace);i++){
        for(int a=0;a<(nage);a++){
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
    std::cout << y << " after LengthAge_alyis_mid" << "\n";
    // Catch at beginning of year
    // Hybrid F tuning inputs & temp storage
    Type v1 = 0.7; Type v2 = 30; Type Fmax = 1.5;
    // for(int fish_flt =0;fish_flt<(nfleets_fish);fish_flt++){}
      // std::cout << y << 	" " << fish_flt << 	" " << catch_yf_obs(y,fish_flt+1) << std::endl;
      for(int fish_flt =0;fish_flt<(nfleets_fish);fish_flt++){
        if(catch_yf_obs(y,fish_flt+1) != -1){
          // std::cout << fish_flt << " F TUNING" << "\n";
          catch_yaf_pred.setZero();
          catch_yf_pred.setZero();
          catch_yfi_pred.setZero();
          catch_yaif_pred.setZero();
          catch_afk_TEMP.setZero();
          Type denom = 0;
          for(int s=0;s<nsex;s++){
            for(int i=0;i<(nspace);i++){
              for(int a=0;a<(nage);a++){
                switch(selType_fish(fish_flt)){
                case 0: // age sel
                    denom += phi_if_fish(fish_flt,i)*
                      fsh_slx_yafs(y,a,fish_flt,s)*
                      N_yais_mid(y,a,i,s)*
                      wtatlen_kab(phi_ik2(i),0)*
                      pow( mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1))+
                      catch_yf_obs(y,fish_flt+1);
                    break;
                  case 1: // length sel
                    denom += phi_if_fish(fish_flt,i)*
                      fsh_slx_yafs(y, mla_yais(y,a,i,s),fish_flt,s)*
                      N_yais_mid(y,a,i,s)*
                      mla_yais(y,a,i,s)*
                      wtatlen_kab(phi_ik2(i),0)*
                      pow( mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1))+
                      catch_yf_obs(y,fish_flt+1);
                  break;
                } // end selType_fish
              } // end age
            } // end space
          } // end sex
          F1_yf(y,fish_flt,0) = catch_yf_obs(y, fish_flt)/denom;
          Type latest_guess = F1_yf(y,fish_flt,0);
          std::cout << y  <<"\t  pre-iter   denom  = " <<  denom  << "\n";
          std::cout << y  <<"\t  pre-iter   F1_yf(y,fish_flt,k)  = " <<  latest_guess  << "\n";
          // k iterations
          for(int k=1;k<(niter+1);k++){
            std::cout << y << "\t" << k <<  "\t" << fish_flt <<  "\t" << "doing k iters" << "\n";
            // modify the guess Eq 20
            Type term0 = 1/(1+exp(v2*( latest_guess - v1)));
            Type term1 = latest_guess*term0;
            Type term2 = v1*(1-term0);
            F1_yf(y,fish_flt,k) = -log(1-(term1+term2));
            std::cout << y << "\t" << k << "\t     F1_yf(y,fish_flt,k)  = " <<   F1_yf(y,fish_flt,k)  << "\n";
            Z_a_TEMP.setZero();
            for(int i=0;i<(nspace);i++){
              for(int a=0;a<(nage);a++){
                for(int s=0;s<nsex;s++){
                  switch(selType_fish(fish_flt)){
                  case 0: // age sel
                    for(int s=0;s<nsex;s++){
                      Z_a_TEMP[a] += fsh_slx_yafs(y, a, fish_flt, s)*F1_yf(y,fish_flt,k) + mat_age(a);
                    } // end sex for z a temp
                    catch_afk_TEMP(a,fish_flt,k) +=
                      F1_yf(y,fish_flt,k)/Z_a_TEMP[a]*
                      (1-exp(-Z_a_TEMP[a]))*
                      phi_if_fish(fish_flt,i)*
                      fsh_slx_yafs(y,a,fish_flt,s)*
                      N_yais_mid(y,a,i,s) *
                      wtatlen_kab(phi_ik2(i),0) *
                      pow(mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1));
                    break;
                  case 1: // length sel
                    for(int s=0;s<nsex;s++){
                      Z_a_TEMP[a] += fsh_slx_yafs(y, mla_yais(y,a,i,s), fish_flt, s)*F1_yf(y,fish_flt,k) + mat_age(a);
                    } // end sex for z a temp
                    // for(int l=0;l<(LBins);l++){
                    catch_afk_TEMP(a,fish_flt,k) +=
                      F1_yf(y,fish_flt,k)/Z_a_TEMP[a]*
                      (1-exp(-Z_a_TEMP[a]))*
                      phi_if_fish(fish_flt,i)*
                      fsh_slx_yafs(y,mla_yais(y,a,i,s),fish_flt,s)*
                      N_yais_mid(y,a,i,s)*
                      mla_yais(y,a,i,s)*
                      wtatlen_kab(phi_ik2(i),0)*
                      pow(mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1));
                    break;
                  } // end selType_fish
                } // end sex
                std::cout << y << "\t" << k << "\t" << "\t  catch_afk_TEMP[a]  = " << catch_afk_TEMP(a,fish_flt,k) << "\n";
              } // end age
            } // end space
            // std::cout << y << "\t" << k << "\t" << "\t  Z_a_TEMP[0]  = " << Z_a_TEMP(0) << "\n";
            // std::cout << y << "\t" << k << "\t" << "\t  catch_afk_TEMP[a+]  = " << catch_afk_TEMP(nage-1,fish_flt,k) << "\n";
            // std::cout << y << "\t" << k << "\t" << "\t  catch_afk_TEMP[a0]  = " << catch_afk_TEMP(0,fish_flt,k) << "\n";
            std::cout << y << "\t" << k << "\t" << "\t  catch_yf_obs(y,fish_flt+1) = " << catch_yf_obs(y,fish_flt+1) << "\n";
            
            vector<Type>Adj(niter+1);
            for(int a=0;a<(nage);a++){
              Adj(k) += catch_yf_obs(y,fish_flt+1)/catch_afk_TEMP(a,fish_flt,k);
            }
            std::cout << y << "\t" << k << "\t" << "\t  Adjk  = " <<   Adj(k)  << "\n";
            
            // Get new Z given ADJ - need to add discard here and switch selex
            Z_a_TEMP2.setZero();
            for(int a=0;a<(nage);a++){
              for(int s=0;s<nsex;s++){
                Z_a_TEMP2(a) += Adj(k)  *
                  fsh_slx_yafs(y, a, fish_flt, s) * F1_yf(y, fish_flt, k) +
                  mat_age(a);
              } // end sex
            } // end age
            std::cout << y << "\t" << k << "\t" << "\t  Z_a_TEMP2[0]  = " <<   Z_a_TEMP2(0)  << "\n";
            
            // Second Guess for F (EQ 24)
            Type denom = 0;
            for(int s=0;s<nsex;s++){
              for(int i=0;i<(nspace);i++){
                for(int a=0;a<(nage);a++){
                  switch(selType_fish(fish_flt)){
                  case 0: // age sel
                    denom += phi_if_fish(fish_flt,i)*
                      fsh_slx_yafs(y,a,fish_flt,s)*
                      N_yais_mid(y,a,i,s)*
                      wtatlen_kab(phi_ik2(i),0)*
                      pow(mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1))*
                      (1-exp(-Z_a_TEMP2(a))) * (F1_yf(y,fish_flt,k)/(Z_a_TEMP2(a)));
                    break;
                  case 1: // length sel
                    denom += phi_if_fish(fish_flt,i)*
                      fsh_slx_yafs(y,mla_yais(y,a,i,s),fish_flt,s)*
                      N_yais_mid(y,a,i,s)*
                      mla_yais(y,a,i,s)*
                      wtatlen_kab(phi_ik2(i),0)*
                      pow(mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1))*
                      (1-exp(-Z_a_TEMP2(a))) * (F1_yf(y,fish_flt,k)/(Z_a_TEMP2(a)));
                    break;
                  } // end selType_fish
                } // end age
                std::cout << y << "\t" << k << "\t" << i << "\t  space accumdenom for F2  = " <<   denom  << "\n";
              } // end space
            } // end sex
            F2_yf(y, fish_flt, k) = catch_yf_obs(y, fish_flt+1)/denom;
            // std::cout << y << "\t" << k << "\t" << "\t  denom for F2  = " <<   denom  << "\n";
            
            std::cout << y << "\t" << k << "\t    F2_yf(y, fish_flt, k) = " <<    F2_yf(y, fish_flt, k)   << "\n";
            // Modify the guess again Eq 25
            term0 = 1/(1+exp(v2*( F2_yf(y,fish_flt,k )- v1*Fmax)));
            term1 = F2_yf(y,fish_flt,k)*term0;
            term2 = v1*(1-term0);
            F2_yf(y, fish_flt, k) = -log(1-(term1+term2));
            latest_guess =    F2_yf(y, fish_flt, k);
            std::cout << y << "\t" << k << "\t latest guess (mod   F2_yf(y, fish_flt, k) )= " << latest_guess << "\n";
          } // end k hybrid F iterations
          // std::cout << y << "\t" << "end k iters" << "\n";
          // Define F, Z and predicted catches
          Freal_yf(y, fish_flt) = F2_yf(y, fish_flt, niter); //final as Freal_yf
          // annoying multi-loops for F in area  get total N exploitable by this fleet
          for(int a=0;a<(nage);a++){
            for(int i=0;i<(nspace);i++){
              for(int s=0;s<nsex;s++){
                N_avail_yf(y,fish_flt) += phi_if_fish(fish_flt, i)*N_yais_mid(y,a,i,s);
              } // end sex
            } // end nspace
          } // end nage
          std::cout << y << "\t" << "end navail" << "\n";
          // get ratio of N in area & reweight F; will just return Freal and 1 for single-area fisheries
          for(int i=0;i<(nspace);i++){
            for(int a=0;a<(nage);a++){
              for(int s=0;s<nsex;s++){
                N_weight_yfi(y,fish_flt, i) = (phi_if_fish(fish_flt, i)* N_yais_mid(y,a,i,s)) /N_avail_yf(y,fish_flt);
              } // end sex
            } // end age
            F_area_yfi(y,fish_flt,i) = Freal_yf(y, fish_flt) * N_weight_yfi(y,fish_flt, i);
          } // end space
          std::cout << y << "\t" << "end N_weight_yfi" << "\n";
          // add together for mgmt regions
          for(int m=1;m<(nmgmt_reg);m++){
            F_ym(y,m) += phi_fm(fish_flt,m)*Freal_yf(y, fish_flt);
          } // end mgmt regions
        // generate predicted catches
        std::cout << y << "\t" << "end F_ym" << "\n";
        for(int i=0;i<(nspace);i++){
          switch(selType_fish(fish_flt)){
          case 0: // age sel
            for(int a=0;a<(nage);a++){
              Zreal_ya(y,a) += Freal_yf(y, fish_flt) + mat_age(a)/2;
              Zreal_yai(y,a,i) += F_area_yfi(y, fish_flt,i) + mat_age(a)/2;

              for(int s=0;s<nsex;s++){
                catch_yaf_pred(y,a,fish_flt) +=
                  Freal_yf(y, fish_flt)/ Zreal_ya(y,a) *
                  (1-exp(- Zreal_ya(y,a) ))*
                  phi_if_fish(fish_flt,i)*
                  fsh_slx_yafs(y,a,fish_flt,s)*
                  N_yais_mid(y,a,i,s)*
                  wtatlen_kab(phi_ik2(i),0)*
                  pow(mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1));

                catch_yaif_pred(y,a,i,fish_flt) += (F_area_yfi(y,fish_flt,i)/
                  ( Zreal_yai(y,a,i)))*(1-exp(- Zreal_yai(y,a,i)  ))*
                    phi_if_fish(fish_flt, i)*
                    fsh_slx_yafs(y,a,fish_flt,s)*
                    N_yais_mid(y,a,i,s)*
                    wtatlen_kab(phi_ik2(i),0)*
                    pow(mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1));
              } // end sex
            } // end age
            break;
          case 1: // length sel
            for(int a=0;a<(nage);a++){
              Zreal_ya(y,a) += Freal_yf(y, fish_flt) + mat_age(a)/2;
              Zreal_yai(y,a,i) += F_area_yfi(y, fish_flt,i) + mat_age(a)/2;
            } // end age for Z
            // for(int l=0;l<(LBins);l++){
              for(int a=0;a<(nage);a++){
                for(int s=0;s<nsex;s++){
                  catch_yaf_pred(y,a,fish_flt) +=
                    Freal_yf(y, fish_flt)/ Zreal_ya(y,a) *
                    (1-exp(- Zreal_ya(y,a) ))*
                    phi_if_fish(fish_flt,i)*
                    fsh_slx_yafs(y,mla_yais(y,a,i,s),fish_flt,s)*
                    N_yais_mid(y,a,i,s)*
                    mla_yais(y,a,i,s)*
                    wtatlen_kab(phi_ik2(i),0)*
                    pow(mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1));

                  catch_yaif_pred(y,a,i,fish_flt) += (F_area_yfi(y,fish_flt,i)/
                    ( Zreal_yai(y,a,i)))*(1-exp(- Zreal_yai(y,a,i)  ))*
                      phi_if_fish(fish_flt,i)*
                      fsh_slx_yafs(y, mla_yais(y,a,i,s),fish_flt,s)*
                      N_yais_mid(y,a,i,s)*
                      mla_yais(y,a,i,s)*
                      wtatlen_kab(phi_ik2(i),0)*
                      pow(mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1));
                } // end sex
              } // end age
            break;
          } // end selType_fish
          for(int a=0;a<(nage);a++){
            catch_yfi_pred(y,fish_flt,i) += catch_yaif_pred(y,a,i,fish_flt);
            catch_yf_pred(y,fish_flt) += catch_yaf_pred(y,a,fish_flt);
          } // end age
        } // end space
      } // end -1 NA trap
    }// end nfleets_fish
      
    std::cout << y << "END OF NFLEETS FISH F TUNING" << "\n";
           
    // N_yais_end ----
    //fill EOY and beginning of next year using Ztuned
    //this will populate ages 2:nage using the end-of year biomass, which accounts for the remaineder
    //of the mortality and the tuned F extraction.
    
    for(int s=0;s<nsex;s++){
      for(int i=0;i<(nspace);i++){
        for(int a=0;a<(nage);a++){
          N_yais_end(y,a,i,s) = N_yais_mid(y,a,i,s)*exp(-(0.2));
          // N_yais_end(y,a,i,s) = N_yais_mid(y,a,i,s)*exp(-(Zreal_yai(y,a,i)));
        }
        for(int a=1;a<(nage-1);a++){
          N_yais_beg(y+1,a,i,s) = N_yais_end(y,a-1,i,s);
        //   std::cout << "filling N for year " << y+1 << "\t space" << i << "\t age" <<  a <<  N_yais_beg(y+1,a,i,s)  << "\n";
        }
        N_yais_beg(y+1,(nage-1),i,s)= N_yais_end(y,nage-1,i,s) + N_yais_end(y,nage-2,i,s);
        // std::cout << "filling N for year " << y+1 << "\t space" << i << "\t" << "\n";
      } // end subareas i
    } // end sexes
    std::cout << y << " N yais end" << "\n";
    // //reweight length-at-age given movement
    for(int s=0;s<nsex;s++){
      for(int i=0;i<(nspace);i++){
        for(int a=0;a<(nage);a++){
          Type LCome = 0.0; Type NCome = 0.0;
          for(int j=0;j<(nspace);j++){
            if(i != j){
              LCome = LCome + phi_ij(i,j)*N_yais_end(y,a,j,s)*Length_yais_mid(y,a,j,s); // for numerator
              NCome = NCome + phi_ij(i,j)*N_yais_end(y,a,j,s); // for denom
            }
          } // end subareas j
          Length_yais_end(y,a,i,s) =    (N_yais_end(y,a,i,s)*Length_yais_mid(y,a,i,s) + LCome)/
            (N_yais_end(y,a,i,s)+NCome);
          Length_yais_beg(y+1,a,i,s) = (N_yais_end(y,a,i,s)*Length_yais_mid(y,a,i,s) + LCome)/
            (N_yais_end(y,a,i,s)+NCome);

        } // end ages
      } // end subareas i
    } // end sexes
        std::cout << y << " reweight length-at-age given movement" << "\n";
    // // SSB_yi, SSB_yk
    for(int i=0;i<(nspace);i++){
      for(int a=0;a<(nage);a++){
        SSB_yi(y,i) += N_yais_end(y,a,i,0)*
          wtatlen_kab(phi_ik2(i),0)*
          pow(mla_yais(y,a,i,0),wtatlen_kab(phi_ik2(i),1))*
          mat_ak(a,phi_ik2(i));
      } // end ages
    } // end space
    std::cout << y << "\t" << "end SSB_yi" << "\n";
    for(int k=0;k<(nstocks);k++){
      for(int i=0;i<(nspace);i++){
        SSB_yk(y,k) +=  phi_ki(k,i)*SSB_yi(y,i);
      } // end stocks
    } // end space
    std::cout << y << "\t" << "end SSB_yk" << "\n";
    for(int m=0;m<(nmgmt_reg);m++){
      for(int i=0;i<(nspace);i++){
        SSB_ym(y,m) += phi_im(i,m)*SSB_yi(y,i);
      } // end space
    } //end mgmt
    std::cout << y << "\t" << "end SSB_ym" << "\n";
    // R_yi, R_yk
    for(int k=0;k<(nstocks);k++){
      // SSB_yk already has summation
      R_yk(y,k) = (4*h_k(k)*R_0k(k)*SSB_yk(y,k))/
        (SSB_0k(k)*(1-h_k(k))+
          SSB_yk(y,k)*(5*h_k(k)-1))*exp(-0.5*b[y]*SDR*SDR+tildeR_yk(y,k));
    }  // end stocks
    std::cout << y << "\t" << "end R_yk" << "\n";
    for(int i=0;i<(nspace);i++){
      R_yi(y,i) = R_yk(y,phi_ik2(i))*tau_ki(phi_ik2(i),i);//*omega_0ij(i); /// downscale to subarea including age-0 movement
      N_yais_beg(y+1,0,i,0) = 0.5*R_yi(y,i);
      N_yais_beg(y+1,0,i,1) = 0.5*R_yi(y,i);
    } /// end space
    std::cout << y << "\t" << "end R_yi" << "\n";
    for(int m=0;m<(nmgmt_reg);m++){
      for(int i=0;i<(nspace);i++){
        R_ym(y,m) += phi_im(i,m)*R_yi(y,i);
      } // end space
    } //end mgmt
    std::cout << y << "\t" << "end R_ym" << "\n";
    // Estimate survey biomass at midyear
    for(int i=0;i<(nspace);i++){
      for(int s=0;s<nsex;s++){
        for(int surv_flt =0;surv_flt<(nfleets_surv);surv_flt++){
          if(surv_yf_obs(y,surv_flt) != -1){
            switch(selType_surv(surv_flt)){
            case 0: // age sel
              for(int a=0;a<nage;a++){
                surv_yf_pred(y,surv_flt) += q_f(surv_flt)*
                  srv_slx_yafs(y,a,surv_flt,s)*
                  phi_if_surv(surv_flt,i)*
                  N_yais_mid(y,a,i,s)*
                  wtatlen_kab(phi_ik2(i),0)*
                  pow(mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1));

                Nsamp_acomp_yf(y,surv_flt) +=  srv_slx_yafs(y,a,surv_flt,s)*
                  phi_if_surv(surv_flt,i)*
                  N_yais_mid(y,a,i,s);
              } // end ages
              break;
            case 1:
                for(int a=0;a<(nage);a++){
                  surv_yf_pred(y,surv_flt) +=  q_f(surv_flt)*
                    srv_slx_yafs(y,mla_yais(y,a,i,s),surv_flt,s)*
                    phi_if_surv(surv_flt,i)*
                    N_yais_mid(y,a,i,s)*
                    mla_yais(y,a,i,s)*
                    wtatlen_kab(phi_ik2(i),0)*
                    pow(mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1));

                  Nsamp_acomp_yf(y,surv_flt) +=  srv_slx_yafs(y,mla_yais(y,a,i,s),surv_flt,s)*
                    phi_if_surv(surv_flt,i)*
                    N_yais_mid(y,a,i,s);
              }
              break;
            } // end selType_fish
          } // end check that it's not an NA year
        } // end survey fleets
      } // end sexes
    } // end nspace
    std::cout << y << "\t" << "end surv_yf_pred" << "\n";
    // predicted age comps, given error
    for(int acomp_flt = 0;acomp_flt<(nfleets_acomp);acomp_flt++){
      // age 0
      acomp_yaf_temp(y,0,acomp_flt) = pnorm(age(0), age_error(phi_fm_acomp2(acomp_flt),0), age_error_SD(phi_fm_acomp2(acomp_flt),0));
      // Loop over ages
      for(int a=0;a<(nage-1);a++){
        acomp_yaf_temp(y,a,acomp_flt) =
          pnorm(Type(a+1),   age_error(phi_fm_acomp2(acomp_flt),a),  age_error_SD(phi_fm_acomp2(acomp_flt),a)) -
          pnorm(age(a),   age_error(phi_fm_acomp2(acomp_flt),a),  age_error_SD(phi_fm_acomp2(acomp_flt),a));
      } // end ages
      acomp_yaf_temp(y,nage-1,acomp_flt) = Type(1.0) - pnorm(Type(nage-1),
                     age_error(phi_fm_acomp2(acomp_flt),nage-1),
                     age_error_SD(phi_fm_acomp2(acomp_flt),nage-1));

      // calculate nsamp within this loop
      for(int i=0;i< nspace;i++){
        for(int s=0;s<nsex;s++){
          for(int a=0;a<(nage-1);a++){
            if(acomp_flt_type(acomp_flt) == 0){
              switch(selType_fish(phi_ff_acomp(acomp_flt,0))){
              case 0: // age sel fish fleet
                Nsamp_acomp_yf(y,phi_ff_acomp(acomp_flt,2)) +=
                  fsh_slx_yafs(y,a,phi_ff_acomp(acomp_flt,0),s)*
                  phi_if_acomp(acomp_flt,i)*
                  N_yais_mid(y,a,i,s);
                break;
              case 1: // len sel fish fleet
                  Nsamp_acomp_yf(y,phi_ff_acomp(acomp_flt,2)) +=
                    fsh_slx_yafs(y,mla_yais(y,a,i,s),phi_ff_acomp(acomp_flt,0),s)*
                    phi_if_acomp(acomp_flt,i)*
                    mla_yais(y,a,i,s)*
                    N_yais_mid(y,a,i,s);
                break;
              } //end selType switch for comms
            }else{
              if(selType_surv(phi_ff_acomp(acomp_flt,1)) == 0){
                Nsamp_acomp_yf(y,phi_ff_acomp(acomp_flt,2)) +=
                  srv_slx_yafs(y,a,phi_ff_acomp(acomp_flt,1),s)*
                  phi_if_acomp(acomp_flt,i)*
                  N_yais_mid(y,a,i,s);
              }else{
                  Nsamp_acomp_yf(y,phi_ff_acomp(acomp_flt,2)) +=
                    srv_slx_yafs(y,mla_yais(y,a,i,s),phi_ff_acomp(acomp_flt,1),s)*
                    phi_if_acomp(acomp_flt,i)*
                    mla_yais(y,a,i,s)*
                    N_yais_mid(y,a,i,s);
              } // end selType switch for survs
            } // end fltType switch
          } // end ages for nsamp
        } // end sex loop for nsamp
      } // end nspace for nsamp
      for(int a=0;a<(nage);a++){
        for(int i=0;i<(nspace);i++){
          for(int s=0;s<nsex;s++){
            if(acomp_flt_type(acomp_flt) == 0){
              switch(selType_fish(phi_ff_acomp(acomp_flt,0))){
              case 0: // age sel fish fleet
                comm_acomp_yafs_pred(y,a,phi_ff_acomp(acomp_flt,3),s) +=
                  acomp_yaf_temp(y,a,acomp_flt)*
                  fsh_slx_yafs(y,a,acomp_flt,s)*
                  phi_if_acomp(acomp_flt,i)*
                  N_yais_mid(y,a,i,s)/  Nsamp_acomp_yf(y,phi_ff_acomp(acomp_flt,2));
                break;
              case 1: // len sel fish fleet
                  comm_acomp_yafs_pred(y,a,phi_ff_acomp(acomp_flt,3),s) +=
                    acomp_yaf_temp(y,a,acomp_flt)*
                    fsh_slx_yafs(y,mla_yais(y,a,i,s),phi_ff_acomp(acomp_flt,0),s)*
                    phi_if_acomp(acomp_flt,i)*
                    mla_yais(y,a,i,s)*
                    N_yais_mid(y,a,i,s)/  Nsamp_acomp_yf(y,phi_ff_acomp(acomp_flt,2));
                break;
              } //end selType switch for comms
            }else{
              if(selType_surv(phi_ff_acomp(acomp_flt,1)) == 0){
                surv_acomp_yafs_pred(y,a,phi_ff_acomp(acomp_flt,4),s) +=
                  acomp_yaf_temp(y,a,acomp_flt)*
                  srv_slx_yafs(y,a,phi_ff_acomp(acomp_flt,1),s)*
                  phi_if_acomp(acomp_flt,i)*
                  N_yais_mid(y,a,i,s)/  Nsamp_acomp_yf(y,phi_ff_acomp(acomp_flt,2));
              }else{
                for(int l=1;l<(LBins);l++){
                  surv_acomp_yafs_pred(y,a,phi_ff_acomp(acomp_flt,4),s) +=
                    acomp_yaf_temp(y,a,acomp_flt)*
                    srv_slx_yafs(y,mla_yais(y,a,i,s),phi_ff_acomp(acomp_flt,1),s)*
                    phi_if_acomp(acomp_flt,i)*
                    mla_yais(y,a,i,s)*
                    N_yais_mid(y,a,i,s)/  Nsamp_acomp_yf(y,phi_ff_acomp(acomp_flt,2));
                } // end lbins
              } // end seltype surv ifelse
            } // end acomp fleet type
          } // end sex
        } // end space
      } // end age
    } // end acomp fleets
  } // END YEARS; END MODEL RUN

  
  // // LIKELIHOODS //
  // Likelihood: survey biomass
  Type ans_survey=0.0;
  for(int surv_flt =0;surv_flt<(nfleets_surv);surv_flt++){
    for(int y=0;y<tEnd;y++){ // Survey Surveyobs
      if(surv_yf_obs(surv_flt) != -1){
        ans_survey -= dnorm(log(surv_yf_pred(y,surv_flt)),
                            log(surv_yf_obs(y,surv_flt)),
                            surv_yf_err(y,surv_flt), TRUE); 
      } // end survey for neg 1
    } // end y
  } // end surv_flt
  
  // Likelihood: catches
  Type ans_catch = 0.0;
  for(int fish_flt =0;fish_flt<(nfleets_fish);fish_flt++){
    for(int y=0;y<tEnd;y++){
      if(catch_yf_obs(fish_flt) != -1){
        ans_catch -= dnorm(log(catch_yf_pred(y,fish_flt)),
                           log(catch_yf_obs(y,fish_flt)),
                           catch_yf_error(y,fish_flt), TRUE); 
      }
    } // end y
  } // end fish_flt
  
  // Likelihood: age comps in surveys & catches
  Type ans_survcomp = 0.0;
  Type ans_catchcomp = 0.0;
  vector<Type>sum1(tEnd); // survey comp likelihood
  vector<Type>sum2(tEnd); // fishery comp likelihood
  sum1.setZero();
  sum2.setZero();
  for(int acomp_flt = 0;acomp_flt<(nfleets_acomp);acomp_flt++){
    for(int y=1;y<tEnd;y++){ // Loop over available years      
      for(int s=0;s<nsex;s++){
        for(int a=0;a<nage;a++){ // Loop over other ages (first one is empty for survey)
          if(acomp_yafs_obs(y,a,acomp_flt,s) == -1){ // Flag if  there was a measurement that year
            sum1(y) += lgamma(Nsamp_acomp_yf(y,phi_ff_acomp(acomp_flt,2))*acomp_yafs_obs(y,a,acomp_flt,s)+1);
            if(acomp_flt_type(acomp_flt) == 0){
              sum2(y) += lgamma(Nsamp_acomp_yf(y,phi_ff_acomp(acomp_flt,2))*
                acomp_yafs_obs(y,a,acomp_flt,s) +
                pi_acomp(acomp_flt)*
                Nsamp_acomp_yf(y,phi_ff_acomp(acomp_flt,2))*
                comm_acomp_yafs_pred(y,a,phi_ff_acomp(acomp_flt,3),s)) -
                - lgamma(pi_acomp(acomp_flt)*
                Nsamp_acomp_yf(y,phi_ff_acomp(acomp_flt,2))*
                comm_acomp_yafs_pred(y,a,phi_ff_acomp(acomp_flt,3),s));
              
              ans_catchcomp += lgamma(Nsamp_acomp_yf(y,phi_ff_acomp(acomp_flt,2))+1)-
                sum1(y)+
                lgamma(pi_acomp(acomp_flt)*Nsamp_acomp_yf(y,phi_ff_acomp(acomp_flt,2)))-
                lgamma(Nsamp_acomp_yf(y,phi_ff_acomp(acomp_flt,2))+
                pi_acomp(acomp_flt)*
                Nsamp_acomp_yf(y,phi_ff_acomp(acomp_flt,2)))+
                sum2(y);
              
            } else{
              sum2(y) += lgamma(Nsamp_acomp_yf(y,phi_ff_acomp(acomp_flt,2))*
                acomp_yafs_obs(y,a,acomp_flt,s) +
                pi_acomp(acomp_flt)*
                Nsamp_acomp_yf(y,phi_ff_acomp(acomp_flt,2))*
                surv_acomp_yafs_pred(y,a,phi_ff_acomp(acomp_flt,4),s))-
                lgamma(pi_acomp(acomp_flt)*
                Nsamp_acomp_yf(y,phi_ff_acomp(acomp_flt,2))*
                surv_acomp_yafs_pred(y,a,phi_ff_acomp(acomp_flt,4),s));
              
              ans_survcomp += lgamma(Nsamp_acomp_yf(y,phi_ff_acomp(acomp_flt,2))+1)-
                sum1(y)+
                lgamma(pi_acomp(acomp_flt)*Nsamp_acomp_yf(y,phi_ff_acomp(acomp_flt,2)))-
                lgamma(Nsamp_acomp_yf(y,phi_ff_acomp(acomp_flt,2))+
                pi_acomp(acomp_flt)*
                Nsamp_acomp_yf(y,phi_ff_acomp(acomp_flt,2)))+
                sum2(y);
            } // end switch for comm or surv type
          } // end acomp flag
        } // end age
      } // end sex
    } // end y
  } // end acomp fleets
  
  // Likelihood: SD Recruitment (hyperprior)
  Type ans_SDR = 0.0;
  for(int k=0;k<(nstocks);k++){
    for(int y=0;y<(tEnd-1);y++){ // Start y loop
      ans_SDR += Type(0.5)*(tildeR_yk(y,k)*tildeR_yk(y,k))/(SDR*SDR)+b(y)*log(SDR*SDR);
    }
  }
  // Likelihood: Error for Selectivity
  Type ans_psel = 0.0;
  // for(int y=0;y<year_sel;y++){ // Start y loop
  //   for(int i=0;i<psel_fish.size();i++){ // Start y loop
  //     ans_psel += Type(0.5)*(PSEL(i,y)*PSEL(i,y))/(sigma_psel*sigma_psel);
  //   }
  // }
  
  // // Likelihood: Priors on h and M
  Type ans_priors = 0.0;
  // for(int i=0;i<(nspace);i++){
  //   for(int a=0;a<(nage-1);a++){ // needs to loop over all years of inits
  //     for(int s=0;s<nsex;s++){
  //       ans_priors += Type(0.5)*(Ninit_ais(a,i,s)*Ninit_ai(a,i,s))/(SDR*SDR);
  //     } // end sex
  //   } // end ages
  // } // end space
  
  // Likelihood: Prior on h
  // ans_priors += -dnorm(logh,log(Type(0.777)),Type(0.113),TRUE);
  // ans_priors += -dbeta(h,Bprior,Aprior,TRUE);
  // if(sum_zero == 1){
  //   ans_priors += ((Type(0.0)-sum(tildeR_yk))*(Type(0.0)-sum(tildeR_yk)))/Type(0.01);
  // }
  // ans_priors += -dnorm(logMinit, log(Type(0.2)), Type(0.1), TRUE);
  // ans_priors += 0.5*pow(logMinit-log(Type(0.2)),2)/Type(0.01);
  // 
  vector<Type>ans_tot(7);
  ans_tot(0) = ans_SDR;
  ans_tot(1) = ans_psel;
  ans_tot(2) = ans_catch;
  ans_tot(3) = ans_survey;
  ans_tot(4) = ans_survcomp;
  ans_tot(5) = ans_catchcomp;
  ans_tot(6) = ans_priors;
  // 
  // // Likelihood: TOTAL
  Type ans = 1.0;
  // ans_SDR+
  // ans_psel+
  // ans_catch;//+
  // ans_survey-
  // ans_survcomp-
  // ans_catchcomp+
  // ans_priors;
  // Type ans = 0.0;
  // Report calculations
  
  // numbers @ age
  REPORT(Ninit_ais);
  REPORT(N_0ais);
  REPORT(N_yais_beg);
  REPORT(N_yais_mid);
  REPORT(N_yais_end);
  
  // len at age
  REPORT(Length_yais_beg);
  REPORT(Length_yais_mid);
  REPORT(Length_yais_end);
  REPORT(LengthAge_alyis_beg);
  REPORT(LengthAge_alyis_mid);
  REPORT(LengthAge_alyis_end);
  
  // SSB and recruits
  REPORT(SSB_yi);
  REPORT(SSB_ym);
  REPORT(SSB_yk);
  REPORT(SSB_0i);
  REPORT(SSB_0k);
  REPORT(R_yi);
  REPORT(R_ym);
  REPORT(R_yk);
  REPORT(R_0k);
  
  // catches and tuning
  REPORT(catch_yaf_pred);  
  REPORT(catch_yf_pred);  
  REPORT(catch_yfi_pred);  
  REPORT(catch_yaif_pred);  
  REPORT(Freal_yf);
  REPORT(F1_yf)
  REPORT(F2_yf);
  REPORT(Zreal_yai);
  REPORT(F_area_yfi);
  
  // survey biomass
  REPORT(surv_yf_pred);
  
  // age comps
  REPORT(comm_acomp_yafs_pred);
  REPORT(surv_acomp_yafs_pred);
  REPORT(Nsamp_acomp_yf);
  
  // other stuff

  
  // REPORT PARS
  ADREPORT(logR_0k);
  ADREPORT(omega_0ij);
  ADREPORT(logh_k);
  ADREPORT(logq_f);
  REPORT(tildeR_yk);
  REPORT(tildeR_initk);
  REPORT(ans_tot);
  return ans;
}
