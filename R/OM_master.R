## OM_Master.R
## M S KAPUR 
## Generate data, condition and export based on OM config, calling functions.R and TMB model XX

## Load packages and functions
source('./R/presets.R') 
mapply(source, list.files("./R/functions", pattern = ".R", full.names=TRUE))

## Load real data and pass to build_dat, build_pars to get in TMB-friendly format
surv <- read.csv("input/Indices_SS3_2020-01-23v3.csv")  %>% filter(Fleet != "AllAreas") ## VAST stdization
NsurveyFleets <- length(unique(surv$Fleet))
landings <- read.csv("clean_landings.csv") ## made in dataprep.R
NfisheryFleets <- length(unique(landings$Fleet))
agecomps <-  read.csv("clean_agecomps.csv") ## made in dataprep.R
lencomps <-  read.csv("clean_agecomps.csv") ## made in dataprep.R
discard <- read.csv("clean_discard.csv") ## made in dataprep.R
NFleets <- sum(NsurveyFleets,NfisheryFleets)

## Structure
Nareas <- 7 ## number of modeled areas
startyr <- min(ts$year)                   # model start year
endyr <- YEAR <- max(ts$year)           # end year
nyears <- length(syr:lyr)                # number of years        
rec_age <- min(waa$age)               # recruitment age                  
plus_group <- max(waa$age)            # plus group age
nage <- length(rec_age:plus_group)    # number of ages
nlenbin <- length(unique(len$length_bin)) # number of length bins
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
data <- build_data()
parameters <- build_parameters()
random_vars <- build_random_vars()


