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
reps <- obj$report() ## return values with uncertainty
# for (k in 1:3) opt <- TMBhelper::fit_tmb(obj) ## estimate
proc.time()-p
# # 
# reps$N_yais_beg[1:7,c(0:4,71),,1]
# reps$N_yais_mid[1:7,c(0:4,71),,1]
# reps$N_yais_end[1:7,c(0:4,71),,1]
# reps$SSB_yi[1:7,]
# reps$SSB_yk[1:7,]
# reps$R_yk[1:7,]
# reps$R_yi[1:7,]
# # 
# reps$catch_afk_TEMP[,8:9,]
# reps$catch_yaf_pred[1:5,,8]
reps$catch_yf_pred[1:5,]
reps$Zreal_yai[1:3,c(0:4,71),] ## no fleets here!
reps$Freal_yf[1:3,]
# reps$F1_yf[1:3,,]
# reps$F2_yf[1:3,,]
# 
# reps$Length_yais_beg[1:3,,,1]
# 
likes <- reps$ans_tot %>% matrix(., ncol = length(.)) %>% data.frame()
names(likes) = c("SDR","PSEL","CATCH","SURVEY","SURVCOMP","CATCHCOMP","PRIORS")
likes
# reps$fsh_slx_yafs[4:6,,5,1]

# compile("C:/Users/public/shire.cpp")
# dyn.load(dynlib("C:/Users/public/shire"))

