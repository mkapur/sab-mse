## OM_Master.R
## M S Kapur 
## Inspiration  from J Sullivan N Jacobsen Summer 2020
## kapurm@uw.edu
# file.copy("C:/Users/public/shire.cpp",here("TMB","shire.cpp"),overwrite = TRUE)

library(TMB)
library(dplyr)
library(reshape2)
library(ggplot2)
library(r4ss)
library(here)
library(ggsidekick)
compile(here("TMB","shireAEP.cpp"))
dyn.load(dynlib(here("TMB","shireAEP")))

source(here("R","functions",'load_files_OM.R'))
df <- load_data_OM(nspace = 6, move = TRUE) ## data that works with OM

# mappy <- list(
#   logh_k = factor(rep(NA, 4)),
#   logR_0k = factor(rep(NA, 4)), ## sum wc = 12
#   omega_0ij = factor(rep(NA,6)),
#   logq_f = factor(rep(NA, 5)),
#   b =  factor(rep(NA, 60)),  
#   logpi_acomp = factor(rep(NA,df$nfleets_acomp)),
#   logSDR = factor(NA),
#   ## structure is fleet x alpha, beta x time block (1 for now)x sex 
#   log_fsh_slx_pars = factor(array(NA, dim = c(df$nfleets_fish,2,1,2)))
#   # log_srv_slx_pars =  factor(array(NA, dim = c( df$nfleets_surv+(df$nfleets_acomp-5),2,1,2)))
# )

## ~90s with full years
system.time(obj <- MakeADFun(df,
                 parameters = df$parms,
                 # map = mappy, ## fix everything; for testing eigen fails
                 checkParameterOrder = TRUE)) 
## up to 30s
system.time(rep1 <- obj$report()) ## one off caclulation using start pars

likes <- rep1$ans_tot %>% matrix(., ncol = length(.)) %>% data.frame()
names(likes) = c("SDR","CATCH","SURVEY","SURVCOMP","CATCHCOMP","PRIORS")
likes

## about 65s for 22 years
system.time(opt <- TMBhelper::fit_tmb(obj)$opt) ## estimate; can repreat for stability)
# for (k in 1:2)  opt <- nlminb(model$env$last.par.best, obj$fn, obj$gr) 
best <- obj$env$last.par.best ## update object with the best parameters
## 81 s
system.time(rep <- sdreport(obj, par = best)) ## re-run & return values at best pars
dat <- obj$report(par = best)

likes <- rep$ans_tot %>% matrix(., ncol = length(.)) %>% data.frame()
names(likes) = c("SDR","PSEL","CATCH","SURVEY","SURVCOMP","CATCHCOMP","PRIORS")
likes

source(here("R","plotting","plotOM_outputs.R"))

opt$opt$par
opt$opt$objective
opt$opt$time_for_MLE
opt$opt$Convergence_check
opt$opt$AIC

steep <- exp(opt$par[1:4])
names(steep) <- paste0("h","_R",4:1)
logR_0 <- opt$par[5:8]
names(logR_0) <- paste0("logR_0","_R",4:1)

# # 
rep1$N_yais_beg[1:7,c(0:4,71),,1]
rep1$N_yais_mid[1:7,c(0:4,71),,1]
rep1$N_yais_end[1:6,c(0:4,71),,1]
rowSums(rep1$N_yais_end[1:25,,,1])
# rep1$SSB_yi[1:7,]
rep1$SSB_yk[1:25,]
rep1$R_yk[1:25,]
# rep1$R_yi[1:7,]
# # 
# rep1$catch_afk_TEMP[,8:9,]
# rep1$catch_yaf_pred[1:5,,8]
rep1$catch_yf_pred[1:10,]
# rep1$Zreal_yai[1:3,c(0:4,71),] ## no fleets here!
rep1$Freal_yf[1:20,]
# rep1$F1_yf[1:3,,]
# rep1$F2_yf[1:3,,]
# 
# rep1$Length_yais_beg[1:3,,,1]
# 

# rep1$fsh_slx_yafs[4:6,,5,1]

# compile("C:/Users/public/shire.cpp")
# dyn.load(dynlib("C:/Users/public/shire"))
# lower <- obj$par-Inf
# upper <- obj$par+Inf
# lower[names(lower) == 'logh_k'] <- log(0.0001)
# upper[names(upper) == 'psel_fish' ] <- 5
# upper[names(upper) == 'PSEL'] <- 9
# upper[names(upper) == 'logh'] <- log(0.999)
# upper[names(upper) == 'F0'] <- 2
