## Functions.R
## M S KAPUR adapted from J Sullivan seak_sablefish & Grant Adams
## preset data gen, post-processing and plotting functions for use in *_Master.R scripts

# Make data object (TMB-friendly)
build_dat <- function(...) {
  data <- list(
    nareas = nareas,
    nyr = nyr,
    nage = nage,
    nsex = nsex,
    nlenbin = nlenbin,
    lenbin = unique(len$length_bin), 
    
    # Switch recruitment estimation: 0 = penalized likelihood (fixed sigma_r), 1 =
    # random effects
    random_rec = rec_type,
    
    # Switch for selectivity type: 0 = a50, a95 logistic; 1 = a50, slope logistic
    slx_type = slx_type,
    
    # Swtich for age composition type (hopefully one day length comps too): 0 =
    # multinomial; 1 = Dirichlet-multinomial
    comp_type = comp_type,
    
    # Switch for assumption on SPR equilibrium recruitment. 0 = arithmetic mean
    # (same as Federal assessment), 1 = geometric mean, 2 = median (2 not coded
    # yet)
    spr_rec_type = spr_rec_type,
    
    # Natural mortality
    M_type = M_type,  # Switch for natural mortality: fixed = 0, estimated with prior = 1. 
    p_log_M = log(0.1), # Priors for natural mortality (same as 2016-2019 Federal assessment
    p_sigma_M = 0.1,
    
    # Time varying parameters - each vector contains the terminal years of each time block
    fsh_blks = c(14, max(ts$index)), #  fishery selectivity: limited entry in 1985, EQS in 1994 = c(5, 14, max(ts$year))
    srv_blks = c(max(ts$index)), # no breaks survey selectivity
    
    # Discard mortality rate in the directed fishery (currently either 0 or 0.16,
    # borrowed from the halibut fishery)
    dmr = array(data = ifelse(include_discards == TRUE, 0.16, 0), dim = c(nyr, nage, nsex)),
    
    # Probability of retaining a fish, sex- and age-based
    retention = 
      # 100% retention (assuming no discards)
      if(include_discards == FALSE & nsex == 1) {
        array(data = 1,
              # Number of rows could = time blocks but currently doesn't
              dim = c(1, nage, nsex))
        # Discards, single sex model
      } else if (include_discards == TRUE & nsex == 1) {
        array(data = filter(retention, Sex == "Combined") %>% pull(p),
              dim = c(1, nage, nsex))
      } else { # Discards, sex-structured
        array(data = filter(retention, Sex %in% c("Female","Male")) %>%
                mutate(sex = ifelse(Sex == "Male", 1, 2)) %>% 
                arrange(sex) %>% pull(p), dim = c(1, nage, nsex))
      },
    
    
    # Fxx levels that correspond with log_spr_Fxx in Parameter section
    Fxx_levels = c(0.35, 0.40, 0.50, 0.60, 0.70),
    
    # Priors ("p_" denotes prior)
    p_fsh_q = c(exp(-16), exp(-16)),
    sigma_fsh_q = c(1, 1),
    p_srv_q = exp(-17), 
    sigma_srv_q = 1,
    p_mr_q = 1.0,
    sigma_mr_q = 0.01,
    
    # Weights on likelihood components ("wt_" denotes weight) based on weights in Federal model
    wt_catch = 1.0,
    wt_fsh_cpue = 1.0,
    wt_srv_cpue = 1.0,
    wt_mr = 1.0,
    wt_fsh_age = 1.0,
    wt_srv_age = 1.0,
    wt_fsh_len = 1.0,
    wt_srv_len = 1.0,
    wt_rec_like = 2.0,
    wt_fpen = 0.1,
    wt_spr = 100,
    
    # Catch
    data_catch = ts$catch,
    sigma_catch = pull(ts, sigma_catch),
    
    # Mark-recapture estimates
    nyr_mr = n_distinct(mr, mr),
    yrs_mr = mr %>% distinct(index) %>% pull(),
    data_mr = pull(mr, mr),
    sigma_mr = rep(0.05, n_distinct(mr, mr)), #mr %>% pull(sigma_mr), 
    
    # Fishery CPUE
    nyr_fsh_cpue = fsh_cpue %>% n_distinct(fsh_cpue),
    yrs_fsh_cpue = fsh_cpue %>% distinct(index) %>% pull(),
    data_fsh_cpue = pull(fsh_cpue, fsh_cpue),
    sigma_fsh_cpue = pull(fsh_cpue, sigma_fsh_cpue),
    
    # Survey CPUE 
    nyr_srv_cpue = srv_cpue %>% n_distinct(srv_cpue),
    yrs_srv_cpue = srv_cpue %>% distinct(index) %>% pull(),
    data_srv_cpue = pull(srv_cpue, srv_cpue),
    sigma_srv_cpue = pull(srv_cpue, sigma_srv_cpue),
    
    # Timing in month fractions
    spawn_month = 2/12, # Feb
    srv_month = 7/12,   # Jul
    fsh_month = 8/12,   # Aug
    
    # Proportion mature-at-age - flexible to vary over time if you wanted, where
    # rows would be time blocks or annual variations
    prop_mature = matrix(data = rep(bio$prop_mature, 1),
                         ncol = nage, byrow = TRUE),
    
    # Vector of prop_fem *Need both prop_fem and sex_ratio to accommodate single
    # sex and sex-structured models*. In sex-structured version, the N matrix is
    # already split so you don't need prop_fem in spawning biomass calculation (it
    # will be a vector of 1's).
    prop_fem = if (nsex == 1) {bio$prop_fem} else { rep(1, nage) } ,
    
    # Sex ratio in the survey (matrix of 1's if single sex model so that N matrix
    # doesn't get split up by sex ratio)
    sex_ratio = if (nsex == 1) {matrix(data = 1, ncol = nage)
    } else {matrix(data =  c(c(1 - bio$prop_fem), # Proprtion male
                             bio$prop_fem), # Proportion female
                   ncol = nage, byrow = TRUE)},
    
    # Weight-at-age: currently weight-at-age is averaged over all years. If for
    # whatever reason you want to include multiple time periods for weight-at-age,
    # you would change the number of rows from 1 to the number of time periods or
    # years
    data_fsh_waa =
      if (nsex == 1) { # Single sex model
        filter(waa, Source == "Fishery (sexes combined)") %>%
          pull(weight) %>%  matrix(ncol = nage, nrow = nsex)  %>%
          array(dim = c(1, nage, nsex))} else {
            # Sex-structured - males in first matrix, females in second
            filter(waa, Source %in% c("Fishery (males)", "Fishery (females)")) %>%
              arrange(desc(Sex)) %>% pull(weight) %>% matrix(ncol = nage, nrow = nsex)  %>%
              array(dim = c(1, nage, nsex))},
    
    # Survey weight-at-age used for everything except predicting catch
    data_srv_waa = 
      if (nsex == 1) { # Single sex model
        filter(waa, Source == "Survey (sexes combined)") %>% 
          pull(weight) %>%  matrix(ncol = nage, nrow = nsex)  %>% 
          array(dim = c(1, nage, nsex))} else {
            # Sex-structured - males in first matrix, females in second
            filter(waa, Source %in% c("Survey (males)", "Survey (females)")) %>% 
              arrange(desc(Sex)) %>% pull(weight) %>% matrix(ncol = nage, nrow = nsex)  %>% 
              array(dim = c(1, nage, nsex))},
    
    # Fishery age comps
    nyr_fsh_age = fsh_age %>% distinct(year) %>% nrow(),
    yrs_fsh_age = fsh_age %>% distinct(index) %>% pull(),
    data_fsh_age = fsh_age %>% select(-c(year, index, Source, n, effn)) %>% as.matrix(),
    n_fsh_age = pull(fsh_age, n),        # total sample size
    effn_fsh_age = pull(fsh_age, effn),  # effective sample size, currently sqrt(n_fsh_age)
    
    # Survey age comps
    nyr_srv_age = srv_age %>% distinct(year) %>% nrow(),
    yrs_srv_age = srv_age %>% distinct(index) %>% pull(),
    data_srv_age = srv_age %>% select(-c(year, index, Source, n, effn)) %>% as.matrix(),
    n_srv_age = pull(srv_age, n),        # total sample size
    effn_srv_age = pull(srv_age, effn),  # effective sample size, currently sqrt(n_srv_age)
    
    # Fishery length comps
    nyr_fsh_len = length(unique(fsh_len$year)),
    yrs_fsh_len = fsh_len %>% distinct(index) %>% pull(),
    data_fsh_len =
      if (nsex == 1) { # Single sex model
        array(data = fsh_len %>% filter(Sex == "Sex combined") %>% pull(proportion),
              dim = c(length(unique(fsh_len$year)), nlenbin, nsex)) } else {
                # Sex-structured (make sure males are first)
                array(data = fsh_len %>% filter(Sex != "Sex combined") %>%
                        arrange(desc(Sex), year) %>% pull(proportion),
                      dim = c(length(unique(fsh_len$year)), nlenbin, nsex))},
    
    n_fsh_len =
      if (nsex == 1) { # Single sex model
        array(data = fsh_len %>% filter(Sex == "Sex combined") %>% distinct(year, Sex, n) %>% pull(n),
              dim = c(length(unique(fsh_len$year)), 1, nsex)) } else {
                # Sex-structured (make sure males are first)
                array(data = fsh_len %>% filter(Sex != "Sex combined") %>% 
                        distinct(year, Sex, n) %>%
                        arrange(desc(Sex)) %>% pull(n),
                      dim = c(length(unique(fsh_len$year)), 1, nsex))},
    
    effn_fsh_len =
      if (nsex == 1) { # Single sex model
        array(data = fsh_len %>% filter(Sex == "Sex combined") %>% distinct(year, Sex, effn) %>% pull(effn),
              dim = c(length(unique(fsh_len$year)), 1, nsex)) } else {
                # Sex-structured (make sure males are first)
                array(data = fsh_len %>% filter(Sex != "Sex combined") %>%
                        distinct(year, Sex, effn) %>% 
                        arrange(desc(Sex)) %>% pull(effn),
                      dim = c(length(unique(fsh_len$year)), 1, nsex))},
    
    # Survey length comps
    nyr_srv_len = length(unique(srv_len$year)),
    yrs_srv_len = srv_len %>% distinct(index) %>% pull(),
    data_srv_len =
      if (nsex == 1) { # Single sex model
        array(data = srv_len %>% filter(Sex == "Sex combined") %>% pull(proportion),                              
              dim = c(length(unique(srv_len$year)), nlenbin, nsex)) } else {
                # Sex-structured (make sure males are first)
                array(data = srv_len %>% filter(Sex != "Sex combined") %>%
                        arrange(desc(Sex)) %>% pull(proportion),
                      dim = c(length(unique(srv_len$year)), nlenbin, nsex))},
    n_srv_len =
      if (nsex == 1) { # Single sex model
        array(data = srv_len %>% filter(Sex == "Sex combined") %>% distinct(year, Sex, n) %>% pull(n),
              dim = c(length(unique(srv_len$year)), 1, nsex)) } else {
                # Sex-structured (make sure males are first)
                array(data = srv_len %>% filter(Sex != "Sex combined") %>%
                        distinct(year, Sex, n) %>% 
                        arrange(desc(Sex)) %>% pull(n),
                      dim = c(length(unique(srv_len$year)), 1, nsex))},
    effn_srv_len =
      if (nsex == 1) { # Single sex model
        array(data = srv_len %>% filter(Sex == "Sex combined") %>%  distinct(year, Sex, effn) %>% pull(effn),
              dim = c(length(unique(srv_len$year)), 1, nsex)) } else {
                # Sex-structured (make sure males are first)
                array(data = srv_len %>% filter(Sex != "Sex combined") %>%
                        distinct(year, Sex, effn) %>% 
                        arrange(desc(Sex)) %>% pull(effn),
                      dim = c(length(unique(srv_len$year)), 1, nsex))},
    
    # Ageing error matrix
    ageing_error = ageing_error,
    
    # Fishery age-length transition matrix *FLAG* currently only using the
    # survey-based matrices from the Feds. Use these as a placeholder for now.
    agelen_key_fsh =
      if (nsex == 1) { # Single sex model
        array(data = c(agelen_key_m), # only have sex-specific, use male curve as placeholder
              dim = c(nage, nage, nsex))} else {
                # Sex-structured (make sure males are first)
                array(data = c(agelen_key_m, agelen_key_f),
                      dim = c(nage, nage, nsex))},
    
    # Survey age-length transition matrix
    agelen_key_srv =
      if (nsex == 1) { # Single sex model
        array(data = c(agelen_key_m), # only have sex-specific, use male curve as placeholder
              dim = c(nage, nage, nsex))} else {
                # Sex-structured (make sure males are first)
                array(data = c(agelen_key_m, agelen_key_f),
                      dim = c(nage, nage, nsex))}
  )
  
  return(data)
}

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

# Build parameter bounds
# Original code by Grant Adams, adapted for use with the sablefish model

build_bounds <- function(param_list = NULL, data_list){
  
  # Debug function
  # param_list <- parameters
  # data_list <- data
  
  upper_bnd <- param_list
  lower_bnd <- param_list
  
  # General bounds
  for(i in 1:length(param_list)){
    upper_bnd[[i]] <- replace(upper_bnd[[i]], values = rep(Inf, length(upper_bnd[[i]])))
    lower_bnd[[i]] <- replace(lower_bnd[[i]], values = rep(-Inf, length(lower_bnd[[i]])))
  }
  
  # Natural mortality
  lower_bnd$log_M <- replace(lower_bnd$log_M, values = -3) 
  upper_bnd$log_M <- replace(upper_bnd$log_M, values = -1) 
  
  # Fishery selectivity
  lower_bnd$log_fsh_slx_pars[,,] <- replace(lower_bnd$log_fsh_slx_pars[,,], values = -2) 
  upper_bnd$log_fsh_slx_pars[,,] <- replace(upper_bnd$log_fsh_slx_pars[,,], values = 2) 
  
  # Survey selectivity
  lower_bnd$log_srv_slx_pars[,,] <- replace(lower_bnd$log_srv_slx_pars[,,], values = -2) 
  upper_bnd$log_srv_slx_pars[,,] <- replace(upper_bnd$log_srv_slx_pars[,,], values = 2) 
  
  # Fishery catchability
  lower_bnd$fsh_logq <- replace(lower_bnd$fsh_logq, values = rep(-30, length(lower_bnd$fsh_logq)))
  upper_bnd$fsh_logq <- replace(upper_bnd$fsh_logq, values = rep(5, length(upper_bnd$fsh_logq)))
  
  # Survey catchability
  lower_bnd$srv_logq <- replace(lower_bnd$srv_logq, values = rep(-30, length(lower_bnd$srv_logq)))
  upper_bnd$srv_logq <- replace(upper_bnd$srv_logq, values = rep(5, length(upper_bnd$srv_logq)))
  
  # Mark-recapture catchability
  lower_bnd$mr_logq <- replace(lower_bnd$mr_logq, values = rep(-20, length(lower_bnd$mr_logq)))
  upper_bnd$mr_logq <- replace(upper_bnd$mr_logq, values = rep(1, length(upper_bnd$mr_logq)))
  
  # Recruitment devs
  lower_bnd$log_rec_devs <- replace(lower_bnd$log_rec_devs, values = rep(-10, length(lower_bnd$log_rec_devs)))
  upper_bnd$log_rec_devs <- replace(upper_bnd$log_rec_devs, values = rep(10, length(upper_bnd$log_rec_devs)))
  
  # Initial numbers-at-age devs
  lower_bnd$log_rinit_devs <- replace(lower_bnd$log_rinit_devs, values = rep(-10, length(lower_bnd$log_rinit_devs)))
  upper_bnd$log_rinit_devs <- replace(upper_bnd$log_rinit_devs, values = rep(10, length(upper_bnd$log_rinit_devs)))
  
  # F devs
  lower_bnd$log_F_devs <- replace(lower_bnd$log_F_devs, values = rep(-15, length(lower_bnd$log_F_devs)))
  upper_bnd$log_F_devs <- replace(upper_bnd$log_F_devs, values = rep(15, length(upper_bnd$log_F_devs)))
  
  # SPR F rates
  lower_bnd$log_spr_Fxx <- replace(lower_bnd$log_spr_Fxx, values = rep(-5, length(lower_bnd$log_spr_Fxx)))
  upper_bnd$log_spr_Fxx <- replace(upper_bnd$log_spr_Fxx, values = rep(5, length(upper_bnd$log_spr_Fxx)))
  
  # Fishery age comp Dirichlet-multinomial theta
  lower_bnd$log_fsh_theta <- replace(lower_bnd$log_fsh_theta, values = rep(-5, length(lower_bnd$log_fsh_theta)))
  upper_bnd$log_fsh_theta <- replace(upper_bnd$log_fsh_theta, values = rep(15, length(upper_bnd$log_fsh_theta)))
  
  # Fishery age comp Dirichlet-multinomial theta
  lower_bnd$log_srv_theta <- replace(lower_bnd$log_srv_theta, values = rep(-5, length(lower_bnd$log_srv_theta)))
  upper_bnd$log_srv_theta <- replace(upper_bnd$log_srv_theta, values = rep(15, length(upper_bnd$log_srv_theta)))
  
  # Put bounds together
  bounds <- list(upper= upper_bnd, lower = lower_bnd)
  
  # Make sure inits are within bounds
  if( sum(unlist(bounds$upper) < as.numeric(unlist(param_list))) > 0 | sum(as.numeric(unlist(param_list)) < unlist(bounds$lower)) > 0 ){
    lower_check <- param_list
    upper_check <- param_list
    param_check <- data.frame(matrix(NA, nrow = length(param_list), ncol = 3))
    colnames(param_check) <- c("Parameter", "Lower", "Upper")
    param_check$Parameter <- names(param_list)
    
    for(i in 1:length(param_list)){
      lower_check[[i]] <- param_list[[i]] < lower_bnd[[i]]
      upper_check[[i]] <- param_list[[i]] > upper_bnd[[i]]
      param_check$Lower[i] <- sum(lower_check[[i]])
      param_check$Upper[i] <- sum(upper_check[[i]])
    }
    
    print("Non-zero value indicates error in initial value")
    print(param_check)
    stop("Initial parameter values are not within bounds")
  }
  
  return(bounds)
}