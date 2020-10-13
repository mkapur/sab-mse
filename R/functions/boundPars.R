
lower <- obj$par-Inf
upper <- obj$par+Inf
lower[names(lower) == 'logh_k'] <- log(0.0001)
upper[names(upper) == 'logh_k'] <- log(0.99)
lower[names(lower) == 'logR_0k'] <- log(0.0001)
lower[names(lower) == 'logSDR'] <- log(0.0001)
## lower bound for p1 (a50 or mean)
lower[names(lower) == 'log_fsh_slx_pars'] <- log(0.0001)
## upper bound for p1 (a50 or mean)
upper[names(upper) == 'log_fsh_slx_pars'][c(1:9,19:28)] <- log(70)
## lower bound p2 = a95 (first four fleets)
lower[names(lower) == 'log_fsh_slx_pars'][c(c(1:4,19:22)+df$nfleets_fish)] <- log(30)
## upper bound p2 = a95 (first four fleets)
upper[names(upper) == 'log_fsh_slx_pars'][c(c(1:4,19:22)+df$nfleets_fish)] <- log(70)
## upper bound p2 =sd (fleets 5:9)
upper[names(upper) == 'log_fsh_slx_pars'][c(c(5:9,23:27)+df$nfleets_fish)] <- log(50)

## currently srv slx all logistic with a95, a50
nsurvsel = dim(df$parms$log_srv_slx_pars)[1]
## lower for everything
lower[names(lower) == 'log_srv_slx_pars'] <- log(0.0001)
## lower bound for p2 (a95)
lower[names(lower) == 'log_srv_slx_pars'][c(c(1:nsurvsel,17:(16+nsurvsel))+nsurvsel)] <- log(70)
## upper bound for p1 (a50 or mean)
upper[names(upper) == 'log_srv_slx_pars'][c(c(1:nsurvsel,17:(16+nsurvsel)))] <- log(70)
## upper bound for p2 (a95)
upper[names(upper) == 'log_srv_slx_pars'][c(c(1:nsurvsel,17:(16+nsurvsel))+nsurvsel)]  <- log(70)
lower[names(lower) == 'omega_0ij'] = 0
upper[names(upper) == 'omega_0ij'] = 1

## sanity check

## last five flts p2 should be 10; p2 for first 4 fleets should be > p1
array(exp(upper[names(upper) == 'log_fsh_slx_pars']), dim = dim(df$parms$log_fsh_slx_pars))
##  p2 should be > p1
array(exp(upper[names(upper) == 'log_srv_slx_pars']), dim = dim(df$parms$log_srv_slx_pars))
## all zero and/or
array(exp(lower[names(lower) == 'log_fsh_slx_pars']), dim = dim(df$parms$log_fsh_slx_pars))
array(exp(lower[names(lower) == 'log_srv_slx_pars']), dim = dim(df$parms$log_srv_slx_pars))
# upper[names(upper) == 'PSEL'] <- 9
# upper[names(upper) == 'logh'] <- log(0.999)
# upper[names(upper) == 'F0'] <- 2
