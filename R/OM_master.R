## OM_Master.R
## M S Kapur 
## Inspiration  from J Sullivan N Jacobsen Summer 2020
## kapurm@uw.edu

library(TMB)
library(dplyr)
library(reshape2)
library(ggplot2)
library(r4ss)
library(here)
library(ggsidekick)
dllUSE = c("shire_v2L",'shire_v2L_1')[2]
compile(here("TMB",paste0(dllUSE,".cpp")))
dyn.load(dynlib(here("TMB",dllUSE)))

source(here("R","functions",'load_files_OM.R'))
df <- load_data_OM(nspace = 6, move = TRUE) ## data that works with OM

# df$v1 = 0.99; df$Fmax = 3;
df$v1 = 0.7;  df$Fmax = 1.5;
# df$v1 = 0.65; df$Fmax = 1.15;
df$niter = 20
df$yRun =  30 #df$tEnd-1
# df$mat_age <- rep(1e-3,df$nage)
df$selShape_fish[5:7] <-  -1 ## turn OFF all slx for BC

omega_0ij_map <- matrix(NA, nrow = 6, ncol = 6) ## turn everything off
# omega_0ij_map[1,] <- df$parms$omega_0ij[1,] ## estimate to/from C1 only

mappy <- list(
  # logh_k = factor(rep(NA, 4)),
  # logR_0k = factor(rep(NA, 4)), ## sum wc = 12
  omega_0ij = factor(omega_0ij_map),
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
array(exp(obj$par[names(obj$par)=='log_fsh_slx_pars']), dim = c(9,2,2),
      dimnames = list(df$fltnames_fish))
## up to 30s
system.time(rep1 <- obj$report()) ## one off caclulation using start pars
head(round(rep1$catch_yf_pred/df$catch_yf_obs[,2:10],2),df$yRun)
colSums(rep1$N_0ais) ## should not be super small anywhere
rep1$SSB_0i ## should not be small or negative
rep1$SSB_yi[3,] ## should match SSB0 without fishing
round(rep1$SSB_yi[1:df$yRun,]) ## should not be small or negative
round(rep1$R_yi[1:df$yRun,])
rowSums(rep1$N_yais_end[1:df$yRun,,,1])

# rep1$N_yais_beg[1:7,c(0:4,71),,1]
# rep1$N_yais_mid[1:7,c(0:4,71),,1]
# rep1$N_yais_end[1:7,c(0:4,71),,1]
# likes <- rep1$ans_tot %>% matrix(., ncol = length(.)) %>% data.frame()
# names(likes) = c("SDR","CATCH","SURVEY","SURVCOMP","CATCHCOMP","PRIORS")
# likes
source(here('R','functions','boundPars.R')) ## bound selex etc
array(exp(upper[names(upper)=='log_fsh_slx_pars']), dim = c(9,2,2),
      dimnames = list(df$fltnames_fish))
## about 65s for 22 years
## 7 hours for 44 yrs
## 3hrs for all years with bounds (just catch like)
system.time(opt <-
              TMBhelper::fit_tmb(
                obj,
                lower = lower,
                upper = upper,
                dll = dllUSE,
                getHessian = FALSE,
                control = list(eval.max = 1e6,
                               iter.max = 1e6,
                               rel.tol = 1e-4)
              )$opt) ## estimate; can repreat for stability)
# for (k in 1:2)  opt <- nlminb(obj$env$last.par.best, obj$fn, obj$gr) 
best <- obj$env$last.par.best ## update object with the best parameters
array(round(exp(best[names(best)=='log_fsh_slx_pars'])), dim = c(9,2,2),
      dimnames = list(df$fltnames_fish))
## 81 s
dat <- obj$report(par = best)
head(round(dat$catch_yf_pred/df$catch_yf_obs[,2:10],2),df$yRun)
colSums(dat$N_0ais) ## should not be super small anywhere
dat$SSB_0i ## should not be small or negative
round(dat$R_yi[1:df$yRun,])
round(dat$SSB_yi[1:df$yRun,]) ## should not be small or negative
rowSums(dat$N_yais_beg[1:df$yRun,,,1])
rowSums(dat$N_yais_end[1:df$yRun,,,1])

steep <- exp(opt$par[1:4])
names(steep) <- paste0("h","_R",1:4)
logR_0 <- opt$par[5:8]
names(logR_0) <- paste0("logR_0","_R",1:4)

system.time(rep <- sdreport(obj, par = best)) ## re-run & return values at best pars

likes <- dat$ans_tot %>% matrix(., ncol = length(.)) %>% data.frame()
names(likes) = c("SDR","CATCH","SURVEY","SURVCOMP","CATCHCOMP","PRIORS")
likes
## save everything and plot
cppname = substr(dllUSE,7,nchar(dllUSE))
writeOM(dat=dat,obj = obj, opt = opt, rep=rep, cppname =cppname,
        runname = paste0("-",df$yRun,"y_",cppname,"_M=",df$mat_age[1],"_noF"))

opt2$par
opt2$objective
opt$time_for_MLE
opt2$Convergence_check
opt$AIC



# # 
dat$N_yais_beg[1:7,c(0:4,71),,1]
dat$N_yais_mid[1:15,c(0:4,71),,1]
dat$N_yais_end[1:15,c(0:4,71),,1]
rowSums(dat$N_yais_end[1:25,,,1])
# dat$SSB_yi[1:7,]
dat$SSB_yk[1:25,]
dat$R_yk[1:25,]
# dat$R_yi[1:7,]
# # 
dat$catch_afk_TEMP[,6,]
# dat$catch_yaf_pred[1:5,,8]
dat$catch_yf_pred[1:10,]
dat$Zreal_yai[5,c(0:4,71),] ## no fleets here!
dat$Freal_yf[1:20,]
dat$F1_yf[1:10,,]
# rep1$F2_yf[1:3,,]
# 
# rep1$Length_yais_beg[1:3,,,1]
# 

# rep1$fsh_slx_yafs[4:6,,5,1]

# compile("C:/Users/public/shire.cpp")
# dyn.load(dynlib("C:/Users/public/shire"))

