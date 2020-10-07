
p = proc.time()

proc.time()-p


reps$R_yk
reps$N_yais_end

save(obj, file = here("TMB",paste0('obj_',Sys.Date(),".rdata")))
save(opt, file = here("TMB",paste0('opt_',Sys.Date(),".rdata")))
save(reps, file = here("TMB",paste0('reps_',Sys.Date(),".rdata")))





plot.figures = FALSE # Set true for printing to file 
## OM MODEL INIT ----
set.seed(731)
# Initialize the model parameters. Make a version with movement and no seasons (simple)
# df.simple <- load_data_seasons(nseason = 1, 
#                                nspace = 2, 
#                                bfuture = 0.5, 
#                                movemaxinit = 0.5, 
#                                movefiftyinit =8) # Prepare data for operating model
# # Run the model using 'run.agebased.true.catch()' -- will fail if pars not == nspace
# sim.data.simple <- run.agebased.true.catch(df.simple)
# 
# ## sanity checks
# sim.data.simple$SSB %>% 
#   data.frame() %>%
#   mutate('totalSSB' = X1+X2, year = as.numeric(row.names(.))) %>%
#   melt(id = 'year') %>%
#   ggplot(., aes(x = year, y = value, color = variable)) + 
#   geom_line(size =2 ) +  theme_sleek()


## OM MODEL CONDITIONING ----
# code from runomem and run om condition


# time <- 1
# yrinit <- df$nyear
# ### Run the OM and the EM for x number of years in the MSE 
# ### Set targets for harvesting etc 
# # df$parms$initN <- df$parms$initN*0
# # df$parms$Rin <- df$parms$Rin*0
# # df$F0 <- 0*df$F0
# simyears <- 25 # Project 30 years into the future (2048 that year)
# year.future <- c(df$years,(df$years[length(df$years)]+1):(df$years[length(df$years)]+simyears))
# N0 <- NA
# sim.data <- runOM_datagen(df) ##
# sim.data <- run.agebased.true.catch(df)

# simdata0 <- sim.data # The other one is gonna get overwritten. 


# Plot stuff 
# parms <- getParameters_OM(trueparms = FALSE, df = df)

##  Create a data frame to send to runsabassessment 

# df.new <- create_TMB_data(sim.data, df, history = TRUE)

# parms.new <- parms
# F0 <- rowSums(sim.data$Fout)
# Rdev <- parms$Rin
# parms.new$F0 <- F0
# parms.new$Rin <- Rdev

## testing without movement
# df.new$omega_ai[,] <- 0.5 ## stationary dist
# df.new$X_ija[,] <- 






source(here("R","functions","plotChecks.R"))

# plots returned from plotChecks.R
(pNinit  | pNzero)/(pNage1  | pNage2)
pSRR
pLAA1  | pLAA2

plot(reps$N_0ai[,1]) ## init l at age
points(reps$N_0ai[,2],pch = 19) ## init l at age

plot(reps$N_yai_beg[1,,1]) ## init l at age
points(reps$N_yai_beg[5,,1],pch = 19) ## init l at age
points(reps$N_yai_beg[53,,1],pch = 19) ## init l at age

reps$survey_bio_f_obs
surv0 <- data.frame(cbind(reps$survey_bio_f_est, 1966:2018)) %>% filter(X6 >1978) %>% cbind(.,reps$survey_bio_f_obs)
names(surv0) = c(paste('EST_',1:5),'Year',paste('OBS_',1:5))
surv0 %>%
  melt(id = 'Year') %>%
  mutate(TYPE = substr(variable,1,3), FLEET = substr(variable,6,6)) %>%
  ggplot(., aes(x = Year, y = value, color = TYPE)) +
  geom_line() + theme_sleek() +
  facet_wrap(~FLEET)




points(reps$Length_yai_beg[3,,1],pch = 19) ## init l at age
points(reps$Length_yai_beg[6,,1],pch = 19) ## init l at age

a1 <- 
  reps$LengthAge_alyi_beg[,,1,2] %>%
  melt() %>%
  group_by(Var2) %>% 
  mutate(sumP = sum(value), pbin = value/sumP)
ggplot(a1,aes(x = Var1, y = Var2, fill = pbin)) +
  geom_tile() +
  labs(x = 'age', y = 'len', title = 'year 1', subtitle= 'subarea 2')


reps$F1_yf
reps$term0; reps$term1;reps$term2

## simulate and re-estimate from OM https://kaskr.github.io/adcomp/Simulation.html
sim <- replicate(5, {
  simdata <- obj$simulate(par=obj$par, complete=TRUE)
  obj2 <- MakeADFun(simdata, parms.new, DLL="runsabassessment", silent=TRUE) ## use original start pars
  reps2 <- obj2$report()
  nlminb(obj2$par, obj2$fn, obj2$gr)$par
})

library(lattice)
df <- data.frame(estimate=as.vector(sim), parameter=names(obj$par)[row(sim)])
densityplot( ~ estimate | parameter, data=df, layout=c(3,1))
dim(reps$N_beg)==dim(reps$N_beg2[[1]])
dim(reps$N_mid)==dim(reps$N_mid2[[1]])

reps$surv_pred ## there should be zeros if the fleet doesn't happen in either area, and doubles if in both, ncol == nfleets_suv
ncol(reps$surv_pred) == df.new$nfleets_surv
ncol(reps$Catch) == df.new$nfleets_fish



# ----
# SSB 
plot(df$years,rowSums(sim.data$SSB.weight))
lines(df$years,SSB.ss3*0.5)
lines(df$years,reps$SSB, col = 'red')
# Survey 
plot(df$years,sim.data$survey.true, type ='l')
points(df$years[df$survey_x == 2],df$survey[df$survey_x == 2,])
points(df$years[df$survey_x == 2],df.new$survey[df$survey_x == 2], col = 'red')
points(df$years[df$survey_x == 2],reps$Surveyobs[df$survey_x == 2], col = 'green')


lower <- obj$par-Inf
lower[names(lower) == 'F0'] <- 0.001
upper <- obj$par+Inf
upper[names(upper) == 'PSEL'] <- 9
upper[names(upper) == 'logphi_catch'] <- log(15)
upper[names(upper) == 'logh'] <- log(0.999)
upper[names(upper) == 'F0'] <- 2
upper[names(upper) == 'psel_fish'] <- 3
lower[names(lower) == 'psel_fish'] <- 0.0001

## optimize!
system.time(opt<-nlminb(obj$par,obj$fn,obj$gr,lower=lower,upper=upper,
                        control = list(iter.max = 1e6, 
                                       eval.max = 1e6))) 

rep <- sdreport(obj)
sdrep <- summary(rep)
rep.values<-rownames(sdrep)
nyear <- df$tEnd

R <- data.frame(name = sdrep[rep.values == 'R',1])

SSB <- data.frame(name = sdrep[rep.values == 'SSB',1])
SSB$SE <- sdrep[rep.values == 'SSB',2]
SSB$min <- SSB$name-2*SSB$SE
SSB$max <- SSB$name+2*SSB$SE
SSB$year <- df$years

Catch <- getUncertainty('Catch', df.new, sdrep)
SSB <- getUncertainty('SSB',df.new, sdrep)

plot(SSB$value)
lines(rowSums(sim.data$SSB))
lines(SSB.ss3*0.5, col = 'red')


# Compare the estimated parameters 
df.p <- as.data.frame(rep$par.fixed)
df.p$name <- names(rep$par.fixed)
df.p$idx <- 1:nrow(df.p)
df.p$model <- 'TMB'
names(df.p)[1] <- 'parameter'
df.p2 <- as.data.frame(unlist(parms.new))
df.p2$name <- names(rep$par.fixed)

df.p2$idx <- 1:nrow(df.p2)
names(df.p2)[1] <- 'parameter'
df.p2$model <- 'SS3'
df.plot <- rbind(df.p,df.p2)

# Fix the log values 
idx <- grep('log', df.plot$name)
df.plot$parameter[idx] <- exp(df.plot$parameter[idx])

ggplot(df.plot, aes(x=  idx, y = parameter, color = model))+geom_point()+
  facet_wrap(~name,scales = 'free')+theme_classic()


ss.exp <- rep(NA, df$nyear)
SSB.ss3 <- mod$derived_quants$Value[grep('SSB_1966', mod$derived_quants$Label):grep('SSB_2018', mod$derived_quants$Label)]


df.ss <- data.frame(year = rep(df$years,3), 
                    SSB = c(SSB$value,
                            reps$SSB,
                            SSB.ss3*0.5
                    ),
                    model = rep(c('TMB est','TMB SS3 parms','SS3'), each = df$tEnd)
)

cols <- PNWColors::pnw_palette('Starfish', n = length(unique(df.ss$model)))

p.ssb <- ggplot(df.ss[df.ss$model != 'Obs' & df.ss$model != 'SS3',], aes(x = year, y = SSB*1e-6, color = model))+
  geom_line()+
  geom_point(data = df.ss[df.ss$model %in% c('SS3','Obs'),])+
  theme_classic()+scale_y_continuous('Spawning biomass\n (million tonnes)')+
  scale_color_manual(values = cols)+
  theme(legend.position = c(0.8,0.8),
        legend.title = element_blank())+
  geom_ribbon(data = SSB, aes(ymin = min*1e-6, ymax = max*1e-6, x = year,
                              y= value*1e-6, color = NA, group = NA),
              fill = alpha('gray', alpha = 0.3), color = NA)

p.ssb