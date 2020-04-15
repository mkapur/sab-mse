## OM_Master.R
## M S KAPUR 
## Generate data, condition and export based on OM config, calling functions.R and TMB model XX

## Load packages and functions
source('./R/presets.R') 
mapply(source, list.files("./R/functions", pattern = ".R", full.names=TRUE))

## Load real data and pass to build_dat, build_pars to get in TMB-friendly format
surv <- read.csv("input/cleaned/clean_survey.csv")
NsurveyFleets <- length(unique(surv$fleet))

landings <- read.csv("./input/cleaned/clean_landings.csv") ## made in dataprep.R
discard <- read.csv("./input/cleaned/clean_discard.csv") ## made in dataprep.R
NfisheryFleets <- length(unique(landings$fleet))


agecomps <-  read.csv("./input/cleaned/clean_agecomps.csv") ## made in dataprep.R
lencomps <-  read.csv("./input/cleaned/clean_lencomps.csv") ## made in dataprep.R

NFleets <- sum(NsurveyFleets,NfisheryFleets)

## Structure
Nareas <- 7 ## number of modeled areas
startyr <- min(ts$year)               
endyr <- YEAR <- max(ts$year)           
nyears <- length(syr:lyr)               
rec_age <- min(waa$age)               # recruitment age                  
plus_group <- max(waa$age)            # plus group age
nage <- length(rec_age:plus_group)    # number of ages
nlenbin <- length(unique(len$length_bin)) # number of length bins - poss vary by area
nsex <- 2                 # single sex or sex-structured
nproj <- 1                # projection years *FLAG* eventually add to cpp file, currently just for graphics
include_discards <- TRUE  # include discard mortality, TRUE or FALSE
tmp_debug <- TRUE         # Temporary debug flag, shut off estimation of selectivity pars

# Model switches
rec_type <- 1     # Recruitment: 0 = penalized likelihood (fixed sigma_r), 1 = random effects
slx_type <- 1     # Selectivity: 0 = a50, a95 logistic; 1 = a50, slope logistic
comp_type <- 0    # Age comp likelihood (not currently developed for len comps): 0 = multinomial, 1 = Dirichlet-multinomial
spr_rec_type <- 1 # SPR equilbrium recruitment: 0 = arithmetic mean, 1 = geometric mean, 2 = median (not coded yet)
M_type <- 0       # Natural mortality: 0 = fixed, 1 = estimated with a prior


## Create dataframes, similar to VAST 
data <- makeDat()
parameters <- build_parameters()
random_vars <- build_random_vars()

TMBphase(data, parameters, random = random_vars, 
         model_name = "mod", phase = FALSE, 
         debug = FALSE)
