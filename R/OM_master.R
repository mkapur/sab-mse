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
dllUSE = c("shireAEP1010",'shire_v2')[2]

# compile(here("TMB","shireAEP.cpp"))
# dyn.load(dynlib(here("TMB","shireAEP")))
# compile(here("TMB","shireAEP1010.cpp"))
dyn.load(dynlib(here("TMB",dllUSE)))
# compile(here("TMB","shire_v2.cpp"))
# dyn.load(dynlib(here("TMB",dllUSE)))

source(here("R","functions",'load_files_OM.R'))
df <- load_data_OM(nspace = 6, move = TRUE) ## data that works with OM

# df$v1 = 0.99; df$Fmax = 3;
df$v1 = 0.7;  df$Fmax = 1.5;
# df$v1 = 0.65; df$Fmax = 1.15;
df$niter = 7
df$yRun = 5;# df$yRun = df$tEnd-1

mappy <- list(
  # logh_k = factor(rep(NA, 4)),
  # logR_0k = factor(rep(NA, 4)), ## sum wc = 12
  # omega_0ij = factor(matrix(NA, nrow = 6, ncol = 6)),
  # logq_f = factor(rep(NA, 5)),
  b =  factor(rep(NA, 60))
  # logpi_acomp = factor(rep(NA,df$nfleets_acomp)),
  # logSDR = factor(NA),
  ## structure is fleet x alpha, beta x time block (1 for now)x sex
  # log_fsh_slx_pars = factor(array(NA, dim = c(df$nfleets_fish,2,1,2)))
  # log_srv_slx_pars =  factor(array(NA, dim = c( df$nfleets_surv+(df$nfleets_acomp-5),2,1,2)))
)

## ~90s with full years
system.time(obj <- MakeADFun(df,
                 parameters = df$parms,
                 dll =dllUSE,
                 map = mappy, ## fix everything for testing eigen fails
                 checkParameterOrder = TRUE)) 
## up to 30s
# system.time(rep1 <- obj$report()) ## one off caclulation using start pars
# head(rep1$catch_yf_pred,10)
# likes <- rep1$ans_tot %>% matrix(., ncol = length(.)) %>% data.frame()
# names(likes) = c("SDR","CATCH","SURVEY","SURVCOMP","CATCHCOMP","PRIORS")
# likes
source(here('R','functions','boundPars.R')) ## bound selex etc
## about 65s for 22 years
## 7 hours for 44 yrs
## 3hrs for all years with bounds (just catch like)
system.time(opt <-
              TMBhelper::fit_tmb(
                obj,
                lower = lower,
                upper = upper,
                dll = dllUSE,
                control = list(eval.max = 1e6,
                               iter.max = 1e6)
              )$opt) ## estimate; can repreat for stability)
# for (k in 1:2)  opt <- nlminb(model$env$last.par.best, obj$fn, obj$gr) 
best <- obj$env$last.par.best ## update object with the best parameters
## 81 s
dat <- obj$report(par = best)
head(dat$catch_yf_pred,10)
head(df$catch_yf_obs,10)

system.time(rep <- sdreport(obj, par = best)) ## re-run & return values at best pars

likes <- dat$ans_tot %>% matrix(., ncol = length(.)) %>% data.frame()
names(likes) = c("SDR","CATCH","SURVEY","SURVCOMP","CATCHCOMP","PRIORS")
likes
## save everything and plot
writeOM(dat=dat,obj = obj, opt = opt, cppname = 'v2',
        runname = "-ltop3y1010")

opt2$par
opt2$objective
opt$time_for_MLE
opt2$Convergence_check
opt$AIC

steep <- exp(opt2$par[1:4])
names(steep) <- paste0("h","_R",4:1)
logR_0 <- opt$par[5:8]
names(logR_0) <- paste0("logR_0","_R",4:1)

# # 
dat$N_yais_beg[1:7,c(0:4,71),,1]
dat$N_yais_mid[1:7,c(0:4,71),,1]
dat$N_yais_end[1:6,c(0:4,71),,1]
rowSums(dat$N_yais_end[1:25,,,1])
# dat$SSB_yi[1:7,]
dat$SSB_yk[1:25,]
dat$R_yk[1:25,]
# dat$R_yi[1:7,]
# # 
# dat$catch_afk_TEMP[,8:9,]
# dat$catch_yaf_pred[1:5,,8]
dat$catch_yf_pred[1:10,]
# dat$Zreal_yai[1:3,c(0:4,71),] ## no fleets here!
dat$Freal_yf[1:20,]
# rep1$F1_yf[1:3,,]
# rep1$F2_yf[1:3,,]
# 
# rep1$Length_yais_beg[1:3,,,1]
# 

# rep1$fsh_slx_yafs[4:6,,5,1]

# compile("C:/Users/public/shire.cpp")
# dyn.load(dynlib("C:/Users/public/shire"))

