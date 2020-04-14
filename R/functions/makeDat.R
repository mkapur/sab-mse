
# Make data object (TMB-friendly)
makeDat <- function(...) {
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
