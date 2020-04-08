
# Make parameter list for TMB
build_pars <- function(
  # data,
  # nsex,
  # inits,
  # rec_devs_inits,
  # rinit_devs_inits,
  # Fdevs_inits,
  ...) {
  
  # Parameter starting values
  parameters <- list(
    
    dummy = 0,   # Used for troubleshooting model               
    
    log_M = log(0.1),  # M_type = 0 is fixed, 1 is estimated
    
    # Fishery selectivity - starting values developed using NOAA selectivity
    # curves see data/NOAA_sablefish_selectivities_2017_jys.xlxs
    log_fsh_slx_pars = 
      # Logistic with a50 and a95, data$slx_type = 0, single sex model
      if(data$slx_type == 0 & nsex == 1) {
        array(data = c(log(4.05), log(3.99), # Sexes combined
                       log(5.30), log(5.20)),
              dim = c(length(data$fsh_blks), 2, nsex)) # 2 = npar for this slx_type
        
        # Logistic with a50 and a95, data$slx_type = 0, sex-structured model
      } else if (data$slx_type == 0 & nsex == 2) {
        array(data = c(log(4.19), log(5.12), # Male
                       log(5.50), log(6.30),
                       log(3.91), log(2.87), # Female
                       log(5.20), log(4.15)),
              dim = c(length(data$fsh_blks), 2, nsex)) # 2 = npar for this slx_type
        
        # Logistic with a50 and slope, data$slx_type = 1, single sex model
      } else if (data$slx_type == 1 & nsex == 1) {
        array(data = c(log(4.05), log(3.99),
                       log(2.29), log(2.43)),
              dim = c(length(data$fsh_blks), 2, nsex)) # 2 = npar for this slx_type
        
      } else {  # Logistic with a50 and slope, data$slx_type = 1, sex-structured model
        array(data = c(log(5.12), log(4.22), # male
                       log(2.57), log(2.61),
                       log(2.87), log(3.86), # female
                       log(2.29), log(2.61)),
              dim = c(length(data$fsh_blks), 2, nsex)) }, # 2 = npar for this slx_type
    
    # Survey selectivity - starting values developed using NOAA selectivity curves
    log_srv_slx_pars = 
      # Logistic with a50 and a95, data$slx_type = 0, single sex model
      if(data$slx_type == 0 & nsex == 1) {
        array(data = c(rep(log(3.74), length(data$srv_blks)),
                       rep(log(5.20), length(data$srv_blks))),
              dim = c(length(data$srv_blks), 2, nsex)) # 2 = npar for this slx_type 
        
        # Logistic with a50 and a95, data$slx_type = 0, sex-structured model
      } else if (data$slx_type == 0 & nsex == 2) {
        array(data = c(rep(log(3.73), length(data$srv_blks)), # male
                       rep(log(5.20), length(data$srv_blks)),
                       rep(log(3.74), length(data$srv_blks)), # female
                       rep(log(5.20), length(data$srv_blks))),
              dim = c(length(data$srv_blks), 2, nsex)) # 2 = npar for this slx_type 
        
        # Logistic with a50 and slope, data$slx_type = 1, single sex model
      } else if (data$slx_type == 1 & nsex == 1) {
        array(data = c(rep(log(3.74), length(data$srv_blks)),
                       rep(log(1.96), length(data$srv_blks))),
              dim = c(length(data$srv_blks), 2, nsex)) # 2 = npar for this slx_type 
        
        # Logistic with a50 and slope, data$slx_type = 1, sex-structured model
      } else { 
        array(data = c(rep(log(3.72), length(data$srv_blks)), # male
                       rep(log(2.21), length(data$srv_blks)),
                       rep(log(3.75), length(data$srv_blks)), # female
                       rep(log(2.21), length(data$srv_blks))),
              dim = c(length(data$srv_blks), 2, nsex)) }, # 2 = npar for this slx_type
    
    # Catchability
    fsh_logq = inits %>% filter(grepl("fsh_logq", Parameter)) %>% pull(Estimate), 
    srv_logq = inits %>% filter(grepl("srv_logq", Parameter)) %>% pull(Estimate),
    mr_logq = inits %>% filter(grepl("mr_logq", Parameter)) %>% pull(Estimate),
    
    # Log mean recruitment and deviations (nyr)
    log_rbar = inits %>% filter(Parameter == "log_rbar") %>% pull(Estimate), 
    log_rec_devs = rec_devs_inits,
    
    # Log mean initial numbers-at-age and deviations (nage-2)
    log_rinit = inits %>% filter(Parameter == "log_rinit") %>% pull(Estimate), 
    log_rinit_devs = rinit_devs_inits,
    
    # Variability in rec_devs and rinit_devs
    log_sigma_r = log(0.6), #log(1.2), # Federal value of 1.2 on log scale
    
    # Fishing mortality
    log_Fbar = inits %>% filter(Parameter == "log_Fbar") %>% pull(Estimate),
    log_F_devs = Fdevs_inits,
    
    # SPR-based fishing mortality rates, i.e. the F at which the spawning biomass
    # per recruit is reduced to xx% of its value in an unfished stock
    log_spr_Fxx = inits %>% filter(grepl("spr_Fxx", Parameter)) %>% pull(Estimate), # F35, F40, F50, F60, F70
    
    # Parameter related to effective sample size for Dirichlet-multinomial
    # likelihood used for composition data. Default of 10 taken from LIME model by
    # M. Rudd. Estimated in log-space b/c it can only be positive.
    log_fsh_theta = log(10),   
    log_srv_theta = log(10)
  )
  
  return(parameters)
}

# Build random vars list for estimating random effects

build_random_vars <- function(
  # data
  ...) {
  
  # If you have a single sigma_r that governs the rinits and the rec_devs, in
  # MakeADFun() the random = c("rinits", "rec_devs") not random = "sigma_r". When
  # you're building the map for phases, it's sigma_r that gets muted as an "NA" if
  # it's not estimated as a random effect
  
  # Setup random effects
  random_vars <- c()
  if (data$random_rec == 1) {
    random_vars <- c("log_rec_devs", "log_rinit_devs")
  }
  
  # Fix parameter if sigma_r is not estimated via random effects
  if(data$random_rec == 0) {
    random_vars <- NULL 
  }
  
  return(random_vars)
}
