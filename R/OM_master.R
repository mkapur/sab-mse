## OM_Master.R
## M S Kapur 
## Inspiration & code guidance from J Sullivan, N Jacobsen Summer 2020
## kapurm@uw.edu
rm(list = ls())

# devtools::install_github('kaskr/adcomp', subdir = 'TMB')
library(TMB)
library(dplyr)
library(reshape2)
library(ggplot2)
library(r4ss)
library(here)
library(ggsidekick)
dllUSE = c('shire_v4')[1]
# compile(here("TMB",paste0(dllUSE,".cpp")))
dyn.load(dynlib(here("TMB",dllUSE)))

source(here("R","functions",'load_files_OM.R'))
df <- load_data_OM(nspace = 6, 
                   move = TRUE,
                   b_y_max = 0.109) ## data that works with OM
df$surv_yf_obs[df$surv_yf_obs >0] <-  df$surv_yf_obs[df$surv_yf_obs >0]*1000
df$yRun <- df$tEnd ## number of years to run model
df$parms$mort_k <- c(0.2,0.2,0.2,0.2)
df$Neqn <- buildNeqn(df)
df$parms$b_y <- rep(1,df$tEnd) ## 1 is no ramp (exp(-0.5*B) in recruits; b*lnRy in like))
# df$selShape_surv[4] <- -1 # constant slx for bc vast
## if by is low, the likelihood is weighted more strongly, and the model is given less
## flexibility in generating R_ys in the context of SDRs (aka do a better job of fitting
# data during this period)

array(exp(df$parms$log_srv_slx_pars), dim = c(df$nfleets_surv+df$nfleets_acomp-4,2,max(df$srv_blks_size),2),
      dimnames = dimnames(df$parms$log_srv_slx_pars))


mappy <-
  buildMap(toFix =  c("omega_0ij",
                      "logh_k",
                      "logSDR",
                      # "tildeR_yk",
                      "b_y",
                      "epsilon_tau",
                      "logpi_acomp",
                      "log_fsh_slx_pars",
                      "log_srv_slx_pars",
                    "mort_k"),
           fixFlt = c("all_fsh",
                      "all_srv"))
                    # c( paste0(c(as.character(unlist(df$fltnames_surv)),
                    # as.character(unlist(df$fltnames_acomp)[c(2,4,5)])))[-c(1,4)] )))
# mappy$logh_k <- factor(c(NA,NA,2,3)) ##  fix WC regs
# mappy$b_y <- factor(c(1,rep(NA,59))) ## enable estimation of year 1 b_y ## consider mirroring for these guys
# mappy$tildeR_yk <- factor(sort(rep(1:(length(mappy$tildeR_yk)/df$nstocks), each = df$nstocks))) ## make each area x year mirrored


array(mappy$log_fsh_slx_pars, dim = c(df$nfleets_fish,2,max(df$fsh_blks_size),2),
      dimnames = dimnames(df$parms$log_fsh_slx_pars))
array(mappy$log_srv_slx_pars, dim = c(df$nfleets_surv+df$nfleets_acomp-4,2,max(df$srv_blks_size),2),
      dimnames = dimnames(df$parms$log_srv_slx_pars))

system.time(obj <- MakeADFun(df,
                 parameters = df$parms,
                 dll = dllUSE,
                 # random = "tildeR_y",
                 map = mappy, ## fix everything for testing eigen fails
                 checkParameterOrder = TRUE)) 

system.time(rep1 <- obj$report()) ## one off caclulation using start pars
# dat = rep1;attach(dat)

rep1$N_yais_beg[,c(1,2,65:71),1,1]
# rep1$N_yais_mid[,c(1,2,65:71),1,1]
# rep1$N_yais_end[,c(1,2,65:71),1,1]
rep1$Length_yais_beg[4,c(1,2,65:71),1,1]
rep1$Length_yais_mid[4,c(1,2,65:71),1,1]
rep1$Length_yais_end[4,c(1,2,65:71),1,1]

rep1$surv_yf_pred

bounds <- boundPars(obj,
                    r0_lower = 0, 
                    boundSlx = c(NA,'fsh','srv')[2:3])


## inspect survey bounds in proper format
exp(bounds$srv_bnds_lwr)
exp(bounds$srv_bnds_upr)

# system.time(opt <- nlminb(
#   obj$par,
#   obj$fn,
#   obj$gr,
#   lower = bounds$lower,
#   upper = bounds$upper,
#   hessian = NULL,
#   control = list(eval.max = 1e6, iter.max = 1e6, rel.tol = 1e-4)
# )
# )

## tmbhelper is returning null OPTS
system.time(opt <-nlminb(obj$par, obj$fn, obj$gr,
              lower = bounds$lower,
              upper = bounds$upper))
# 
# system.time(opt <-
#               TMBhelper::fit_tmb(
#                 obj,
#                 lower = bounds$lower,
#                 upper = bounds$upper,
#                 # dll = dllUSE,
#                 getHessian = FALSE,
#                 getsd = FALSE,
#                 control = list(eval.max = 1e6,
#                                iter.max = 1e6,
#                                rel.tol = 1e-4)
#               )$opt) ## estimate; can repeat for stability)


# for (k in 1:2)  opt <- nlminb(obj$env$last.par.best, obj$fn, obj$gr) 
best <- obj$env$last.par.best ## update object with the best parameters
dat <- obj$report(par = best)
dat$surv_yf_pred/df$surv_yf_obs
# dat$catch_yf_pred_total/df$catch_yf_obs[,2:ncol(df$catch_yf_obs)]
## save everything and plot
cppname = substr(dllUSE,7,nchar(dllUSE))
writeOM(justPlots = FALSE,
  dat=dat,
  obj = obj, 
        opt = opt, 
        rep=rep, 
        cppname =cppname, 
        mappy = mappy,
        runname = paste0("-",df$yRun,"y_",
                         cppname,
                         "",
                         "_truncateestBC_VAST",
                         "_Bramp=1.0"))


# system.time(rep <- sdreport(obj, par = best)) ## re-run & return values at best pars
array(exp(best[names(best)== "log_fsh_slx_pars"]), dim = c(2,2,1,2))
array(exp(best[names(best)== "log_srv_slx_pars"]), dim = c(5,2,3,2))
steep <- exp(best[names(best) == 'logh_k']); names(steep) <- paste0("h","_R",1:4);steep
logR_0 <- best[names(best) == 'logR_0k'];names(logR_0) <- paste0("logR_0","_R",1:4);logR_0
epstau <- exp(best[names(best) == 'epsilon_tau']); names(epstau) <- paste0("epstau_",
                                                                      c('C1','C2','B2','B3','A3','A4'));epstau
by <- round(best[names(best) == 'b_y'],3); 
names(by) <- paste0("b_y",1:length((best[names(best) == 'b_y'])));by 

likes <- dat$ans_tot %>% matrix(., ncol = length(.)) %>% data.frame()
names(likes) = c("SDR","CATCH","SURVEY","SURVCOMP","CATCHCOMP","PRIORS")
likes
