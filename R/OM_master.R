## OM_Master.R
## M S Kapur mod N Jacobsen Summer 2020
## kapurm@uw.edu
library(TMB)
library(dplyr)
library(reshape2)
library(ggplot2)
library(r4ss)
library(here)
library(ggsidekick)
source(here("R","functions",'load_files_OM.R'))
compile(here("TMB","runsabassessment.cpp"))
dyn.load(dynlib(here("TMB","runsabassessment")))


## OM MODEL INIT ----
# Initialize the model parameters. Make a version with movement and no seasons (simple)
df.simple <- load_data_seasons(nseason = 1, 
                               nspace = 2, 
                               bfuture = 0.5, 
                               movemaxinit = 0.5, 
                               movefiftyinit =8) # Prepare data for operating model
# Run the model using 'run.agebased.true.catch()' -- will fail if pars not == nspace
sim.data.simple <- run.agebased.true.catch(df.simple)

## sanity checks
sim.data.simple$SSB %>% 
  data.frame() %>%
  mutate('totalSSB' = X1+X2, year = as.numeric(row.names(.))) %>%
  melt(id = 'year') %>%
  ggplot(., aes(x = year, y = value, color = variable)) + 
  geom_line(size =2 ) +  theme_sleek()
  
## OM MODEL CONDITIONING ----
  # runomem and run om condition
set.seed(731)
plot.figures = FALSE # Set true for printing to file 
# Run the simulation model
assessment <- read.csv(here("input","data",'assessment_MLE.csv')) ## I believe this comes from SS3
assessment <- assessment[assessment$year > 1965 &assessment$year < 2018 ,]
Catch.obs <- read.csv(here("input","data",'hake_totcatch.csv'))
df <- load_data_seasons(nspace =2)
df$Catch <- Catch.obs$Fishery
time <- 1
yrinit <- df$nyear
### Run the OM and the EM for x number of years in the MSE 
### Set targets for harvesting etc 
# df$parms$initN <- df$parms$initN*0
# df$parms$Rin <- df$parms$Rin*0
# df$F0 <- 0*df$F0
simyears <- 25 # Project 30 years into the future (2048 that year)
year.future <- c(df$years,(df$years[length(df$years)]+1):(df$years[length(df$years)]+simyears))
N0 <- NA
sim.data <- run.agebased.true.catch(df)
simdata0 <- sim.data # The other one is gonna get overwritten. 

# Plott stuff 

parms <- getParameters_OM(trueparms = TRUE, df = df)

##  Create a data frame to send to runsabassessment 

df.new <- create_TMB_data(sim.data, df, history = TRUE)

parms.new <- parms
F0 <- rowSums(sim.data$Fout)
Rdev <- parms$Rin
parms.new$F0 <- F0
parms.new$Rin <- Rdev


obj <- MakeADFun(df.new,
                 parms.new,
                 DLL= "runsabassessment") # Run the assessment, in TMB folder
reps <- obj$report()
dim(reps$N_beg)==dim(reps$N_beg2[[1]])
dim(reps$N_mid)==dim(reps$N_mid2[[1]])
reps$Nzero3
reps$N_beg2[[1]][1:100,]
reps$SSBzero2
reps$SSB2 ## dim time x nspace
reps$CatchN


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

