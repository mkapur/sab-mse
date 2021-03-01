## OM_Master.R
## Code to condition & forecast 6-area Operating model for Transboundary Sablefish MSE
## M S Kapur kapurm@uw.edu
## Inspiration & code guidance from J Sullivan, N Jacobsen Summer 2020++

rm(list = ls())

# devtools::install_github('kaskr/adcomp', subdir = 'TMB')
library(TMB)
library(dplyr)
library(reshape2)
library(ggplot2)
library(r4ss)
library(here)
library(ggsidekick)
dllUSE = c('shire_v4_1')[1]
TMB::compile(here("TMB",paste0(dllUSE,".cpp")))#, flags = "-Wno-ignored-attributes")
dyn.load(dynlib(here("TMB",dllUSE)))

source(here("R","functions",'load_files_OM.R'))
yr_future <- 0
df <- load_data_OM(nspace = 6, 
                   move = TRUE,
                   yr_future  = yr_future,
                   b_y_max = 0.109) ## data that works with OM
# df$surv_yf_obs[df$surv_yf_obs >0] <-  df$surv_yf_obs[df$surv_yf_obs >0]*1000

df$parms$mort_k <- c(0.2,0.2,0.2,0.2)
df$Neqn <- buildNeqn(df)
df$parms$b_y <- rep(1,df$tEnd) ## 1 is no ramp (exp(-0.5*B) in recruits; b*lnRy in like))
df$F_yf_HCR <- array(0.2, dim = c(ifelse(yr_future == 0,1,yr_future), df$nfleets_fish))
# df$selShape_surv[4] <- -1 # constant slx for bc vast
## if by is low, the likelihood is weighted more strongly, and the model is given less
## flexibility in generating R_ys in the context of SDRs (aka do a better job of fitting
# data during this period)
mappy <-
  buildMap(toFix =  c("omega_0ij",
                      "logh_k",
                      "logSDR",
                      # "tildeR_y",
                      "b_y",
                      "epsilon_tau",
                      "logpi_acomp",
                      "log_fsh_slx_pars",
                      # "log_srv_slx_pars",
                      "mort_k"),
           fixFlt = c("all_fsh"))#,
                      # "all_srv"))
                      # colnames(df$srv_blks)[6] ))

# mappy$logh_k <- factor(c(NA,NA,2,3)) ##  fix WC regs
# mappy$b_y <- factor(c(1,rep(NA,59))) ## enable estimation of year 1 b_y ## consider mirroring for these guys
# mappy$tildeR_yk <- factor(sort(rep(1:(length(mappy$tildeR_yk)/df$nstocks), each = df$nstocks))) ## make each area x year mirrored

system.time(obj <- MakeADFun(df,
                 parameters = df$parms,
                 dll = dllUSE,
                 # random = "tildeR_y",
                 map = mappy, 
                 checkParameterOrder = TRUE)) 

system.time(rep1 <- obj$report());## one off caclulation using start pars

# sim <- obj$simulate(par=obj$par, complete=TRUE) ; 
# sim$comm_acomp_yafs_pred[77,20,1,1]

bounds <- boundPars(obj,
                    r0_lower = 0, 
                    boundSlx = c(NA,'fsh','srv')[2:3])

## inspect survey bounds in proper format
# exp(bounds$srv_bnds_lwr)
# exp(bounds$srv_bnds_upr)

## tmbhelper is returning null OPTS
system.time(opt <-nlminb(obj$par,
                         obj$fn,
                         obj$gr,
              lower = bounds$lower,
              upper = bounds$upper))
# 
system.time(opt <-
              TMBhelper::fit_tmb(
                obj,
                # lower = bounds$lower,
                # upper = bounds$upper,
                # dll = dllUSE,
                getHessian = FALSE,
                getsd = FALSE,
                control = list(eval.max = 1e6,
                               iter.max = 1e6,
                               rel.tol = 1e-4)
              )$opt) ## estimate; can repeat for stability)


# for (k in 1:2)  opt <- nlminb(obj$env$last.par.best, obj$fn, obj$gr) 
best <- obj$env$last.par.best ## update object with the best parameters
dat <- obj$report(par = best)
# dat$surv_yf_pred/df$surv_yf_obs
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
                         "_",
                         "_BC_DesignBased",
                         "_Bramp=1.0"))

## this repeats the entire estimation (conditioning) step for a number of replicates
## each replicate is simulated from your first conditioned OM. 
## teh technique is more applicable to checking stability. really we just need
## multiple reps (seeds) of OM simulations.
sim <- replicate(5, {
  set.seed(runif(1,1,1000)) ## randomize the seed for nrep
  obj$simulate(par=best, complete=TRUE) ## simulate from last obj, using best pars
  ## seed will vary recdevs.
  # simdata0 <- ## input, would be same as rep
  # simdata <- obj$simulate(par=obj$par, complete=TRUE) ## simulate from last obj,
  ## The default parameter values used for the simulation is obj$env$last.par
  # obj2 <- MakeADFun(simdata, df$parms, DLL=dllUSE, silent=TRUE) ## prep new mod with new data and og parms
  ## obj would need to be inclusive of fyears
  # nlminb(obj2$par, obj2$fn, obj2$gr)$par
})
head(dat$SSB_ym)
head(simdata$SSB_ym)
tail(simdata0$surv_yf_obs, yr_future)
simdata0$SSB_ym == simdata$SSB_ym

sim['surv_yf_obs',1] %>% data.frame() %>% tail()
sim['surv_yf_obs',4] %>% data.frame() %>% tail()
sim['surv_yf_obs',2] %>% data.frame() %>% tail()
sim['surv_yf_obs',3] %>% data.frame() %>% tail()
