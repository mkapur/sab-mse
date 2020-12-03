## OM_Master.R
## M S Kapur 
## Inspiration from J Sullivan, N Jacobsen Summer 2020
## kapurm@uw.edu

library(TMB)
library(dplyr)
library(reshape2)
library(ggplot2)
library(r4ss)
library(here)
library(ggsidekick)
dllUSE = c("shire_v3L",'shire_v4L')[2]
# compile(here("TMB",paste0(dllUSE,".cpp")))
dyn.load(dynlib(here("TMB",dllUSE)))

source(here("R","functions",'load_files_OM.R'))
df <- load_data_OM(nspace = 6, move = TRUE) ## data that works with OM
df$yRun =  35 # df$tEnd-1 ## number of years to run model
df$parms$mort_k <- c(0.2,0.2,0.2,0.2)
df$Neqn <- buildNeqn(df)
df$parms$logq_f <- rep(log(1e-5),length(df$parms$logq_f))
mappy <-
  buildMap(toFix =  c("omega_0ij",
                      "epsilon_tau", 
                      # "logh_k"))
                     "mort_k")) #,
           # fixFlt = c("BC_LL","BC_TRAP","BC_TWL")) )

## ~90s with full years
system.time(obj <- MakeADFun(df,
                 parameters = df$parms,
                 dll =dllUSE,
                 map = mappy, ## fix everything for testing eigen fails
                 checkParameterOrder = TRUE)) 

# system.time(rep1 <- obj$report()) ## one off caclulation using start pars

bounds <- boundPars(obj, 
                    r0_lower = 0, 
                    boundSlx = c('fsh','srv')[1])

# with(bounds, array(exp(lower[names(lower)=='log_fsh_slx_pars']), dim = c(7,2,1,2),
#                    dimnames = list(df$fltnames_fish)))
# with(bounds, array(exp(upper[names(upper)=='log_srv_slx_pars']), dim = c(7,2,1,2),
# dimnames = list(df$fltnames_fish)))

system.time(opt <-
              TMBhelper::fit_tmb(
                obj,
                # lower = bounds$lower,
                # upper = bounds$upper,
                dll = dllUSE,
                getHessian = FALSE,
                control = list(eval.max = 1e6,
                               iter.max = 1e6,
                               rel.tol = 1e-4)
              )$opt) ## estimate; can repreat for stability)
  # for (k in 1:2)  opt <- nlminb(obj$env$last.par.best, obj$fn, obj$gr) 
best <- obj$env$last.par.best ## update object with the best parameters
dat <- obj$report(par = best)


steep <- exp(opt$par[names(opt$par) == 'logh_k']); names(steep) <- paste0("h","_R",1:4);steep
logR_0 <- opt$par[names(opt$par) == 'logR_0k'];names(logR_0) <- paste0("logR_0","_R",1:4);logR_0
# epstau <- opt$par[names(opt$par) == 'epsilon_tau']; names(epstau) <- paste0("epstau_",inames)

likes <- dat$ans_tot %>% matrix(., ncol = length(.)) %>% data.frame()
names(likes) = c("SDR","CATCH","SURVEY","SURVCOMP","CATCHCOMP","PRIORS")
likes
## save everything and plot
cppname = substr(dllUSE,7,nchar(dllUSE))
writeOM(dat=dat,obj = obj, opt = opt, rep=rep, cppname =cppname,
        runname = paste0("-",df$yRun,"y_",cppname,"_M=",
                         paste(df$parms$mort_k,collapse="-"),
                         "_survOFF"))


system.time(rep <- sdreport(obj, par = best)) ## re-run & return values at best pars

