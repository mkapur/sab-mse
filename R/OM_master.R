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

mappy <- list(
  logh_k = factor(rep(NA, 4)),
  logR_0k = factor(rep(NA, 4)), ## sum wc = 12
  # omega_0ij = factor(matrix(NA, nrow = nrow(df$parms$omega_0ij), ncol = nrow(df$parms$omega_0ij))),
  logq_f = factor(rep(NA, 5)),
  b =  factor(rep(NA, 60)),  
  logpi_acomp = factor(rep(NA,df$nfleets_acomp)),
  logSDR = factor(NA)#,
  ## structure is fleet x alpha, beta x time block (1 for now)x sex 
  # log_fsh_slx_pars = factor(array(NA, dim = c(df$nfleets_fish,2,1,2))),
  # log_srv_slx_pars =  factor(array(NA, dim = c( df$nfleets_surv+(df$nfleets_acomp-5),2,1,2)))
)
p = proc.time()
obj <- MakeADFun(df,
                 parameters = df$parms,
                 map = mappy, ## fix everything
                 checkParameterOrder = TRUE,
                 DLL= "shireAEP") # Run the assessment, in TMB folder
proc.time()-p

# reps <- obj$report() ## return values with uncertainty
# opt <- TMBhelper::fit_tmb(obj) ## estimate
# compile("C:/Users/public/shire.cpp")
# dyn.load(dynlib("C:/Users/public/shire"))

