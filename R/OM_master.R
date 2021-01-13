## OM_Master.R
## M S Kapur 
## Inspiration & code guidance from J Sullivan, N Jacobsen Summer 2020
## kapurm@uw.edu
rm(list = ls())

library(TMB)
library(dplyr)
library(reshape2)
library(ggplot2)
library(r4ss)
library(here)
library(ggsidekick)
dllUSE = c('shire_v4')[1]
compile(here("TMB",paste0(dllUSE,".cpp")))
dyn.load(dynlib(here("TMB",dllUSE)))

source(here("R","functions",'load_files_OM.R'))
df <- load_data_OM(nspace = 6, 
                   move = TRUE) ## data that works with OM
df$surv_yf_obs[df$surv_yf_obs >0] <-  df$surv_yf_obs[df$surv_yf_obs >0]*1000
df$yRun <- 5# df$tEnd ## number of years to run model
df$parms$mort_k <- c(0.2,0.2,0.2,0.2)
df$Neqn <- buildNeqn(df)

mappy <-
  buildMap(toFix =  c("omega_0ij",
                      # "epsilon_tau",
                      "log_fsh_slx_pars",
                      "log_srv_slx_pars",
                    "mort_k"),
           fixFlt = c("all_fsh", "AK_GOA_SURV", "BC_StRS", "BC_SS"))

array(mappy$log_fsh_slx_pars, dim = c(df$nfleets_fish,2,max(df$fsh_blks_size),2),
      dimnames = dimnames(df$parms$log_fsh_slx_pars))
array(mappy$log_srv_slx_pars, dim = c(df$nfleets_surv+df$nfleets_acomp-4,2,max(df$srv_blks_size),2),
      dimnames = dimnames(df$parms$log_srv_slx_pars))

system.time(obj <- MakeADFun(df,
                 parameters = df$parms,
                 dll = dllUSE,
                 map = mappy, ## fix everything for testing eigen fails
                 checkParameterOrder = TRUE)) 

system.time(rep1 <- obj$report()) ## one off caclulation using start pars
dat = rep1;attach(dat)

rep1$N_yais_beg[,c(1,2,65:71),1,1]
# rep1$N_yais_mid[,c(1,2,65:71),1,1]
# rep1$N_yais_end[,c(1,2,65:71),1,1]
rep1$Length_yais_beg[4,c(1,2,65:71),1,1]
rep1$Length_yais_mid[4,c(1,2,65:71),1,1]
rep1$Length_yais_end[4,c(1,2,65:71),1,1]

bounds <- boundPars(obj,
                    r0_lower = 0, 
                    boundSlx = c(NA,'fsh','srv')[3])

# with(bounds, array(
#   exp(lower[names(lower) == 'log_fsh_slx_pars']),
#   dim = c(7, 2, max(df$srv_blks), 2),
#   dimnames = list(df$fltnames_fish)
# ))
with(bounds, array(exp(lower[names(lower) == 'log_srv_slx_pars']),
                   dim = c(5, 2, max(df$srv_blks_size), 2)))
with(bounds, array(exp(upper[names(upper) == 'log_srv_slx_pars']),
                   dim = c(5, 2, max(df$srv_blks_size), 2)))

system.time(opt <-
              TMBhelper::fit_tmb(
                obj,
                lower = bounds$lower,
                upper = bounds$upper,
                dll = dllUSE,
                getHessian = FALSE,
                control = list(eval.max = 1e6,
                               iter.max = 1e6,
                               rel.tol = 1e-4)
              )$opt) ## estimate; can repeat for stability)
# for (k in 1:2)  opt <- nlminb(obj$env$last.par.best, obj$fn, obj$gr) 
best <- obj$env$last.par.best ## update object with the best parameters
dat <- obj$report(par = best)
dat$surv_yf_pred/df$surv_yf_obs
## save everything and plot
cppname = substr(dllUSE,7,nchar(dllUSE))
writeOM(justPlots = FALSE,
  dat=dat,
        obj = obj, 
        opt = opt, 
        rep=rep, 
        cppname =cppname, 
        mappy = mappy,
        runname = paste0("-",df$yRun,"y_",cppname,
                         "_fixAcompslx",
                         "_epstauon"))


# system.time(rep <- sdreport(obj, par = best)) ## re-run & return values at best pars

array(exp(best[names(best)== "log_srv_slx_pars"]), dim = c(5,2,1,2))
steep <- exp(best[names(best) == 'logh_k']); names(steep) <- paste0("h","_R",1:4);steep
logR_0 <- best[names(best) == 'logR_0k'];names(logR_0) <- paste0("logR_0","_R",1:4);logR_0
epstau <- exp(best[names(best) == 'epsilon_tau']); names(epstau) <- paste0("epstau_",
                                                                      c('C1','C2','B2','B3','A3','A4'));epstau

likes <- dat$ans_tot %>% matrix(., ncol = length(.)) %>% data.frame()
names(likes) = c("SDR","CATCH","SURVEY","SURVCOMP","CATCHCOMP","PRIORS")
likes
