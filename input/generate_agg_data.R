require(dplyr)
require(reshape2)
require(data.table)

# setwd("C:/Users/mkapur/Dropbox/UW/assessments/sab_2019/200.00_base_files_29May/100.00_base_files/hklpot_agg")
## read in CSVs and drop  'trawl'
disrate_cs <- read.csv(here('input','raw_data','catch',
                            "sablefish_OB_DisRatios_cs_2019_Coastwide_trawl_fixed_2019-06-27.csv")) %>% 
  filter(gear3 != 'Trawl')
disrate_ncs <- read.csv(here('input','raw_data','catch',
                             "sablefish_OB_DisRatios_ncs_2019_Coastwide_trawl_fixed_2019-06-27.csv")) %>%
  filter(gear3 != 'Trawl')
avgwt <- read.csv("Sablefish_WCGOP_Lengths_HookAndPot.csv")

## cs/ncs happend to both fleets spearately so weighting has to happen
# just this time it wont be hkl/pot separately

## combine catches
mod <- SS_output("./base_model")
 mod$catch %>% 
   filter(Fleet %in% c(1,2))  %>% 
   select(Yr, Seas, Fleet, Obs, se) %>% 
   group_by(Yr) %>% 
   summarise(Seas = 1, Fleet = 1, Obs = sum(Obs), se = 0.01 ) %>%
   write.csv(.,file = paste0("./generated/hklpot_catch_SS_",Sys.Date(),".csv"),row.names = FALSE) ## input line 52
 


## Get weighted average meanbw; line 653---
avgwt$Wghtd.AVG_W_KG <- avgwt$Wghtd.AVG_W*0.453592
avgwt$CV_2019 <- avgwt$Wghtd.AVG_W.SD/avgwt$Wghtd.AVG_W
meanbw <- matrix(NA, ncol = 8, nrow = length(unique(avgwt$Year))) %>% data.frame(.)
meanbw[,1] <- unique(avgwt$Year) 
meanbw[,2] <- 7
meanbw[,3] <- 1
meanbw[,4] <- 1
meanbw[,5] <- 2
meanbw[,6] <- round(avgwt$Wghtd.AVG_W_KG[avgwt$Gear == 'HookAndPot'],6)
meanbw[,7] <- round(avgwt$CV_2019[avgwt$Gear == 'HookAndPot'],6) ## is CV but called stderr in SS
meanbw[,8] <- paste("#HKL_POT")
names(meanbw) <- c("yr", "month", "fleet", "part", "type", "obs", "stderr", "note")

meanbw %>%
  write.table(.,file = paste0("./generated/meanbw_SS_",Sys.Date(),".csv"),row.names = FALSE, sep = " ", quote = F) # ## input line 653

## Get proportions from retained lbs, ----
## these will be multiplied by discard rate by sector -- only 2011 onward
discard <-  disrate_ncs  %>% 
  select(yr = ryear, "OBS" = Observed_Ratio, 'sd' = StdDev.Boot_Ratio  ) ## raw data
discard_late <- merge(disrate_cs %>% select(ryear,'CS_LBS' = Observed_RETAINED.LBS, "CS_RATIO" = Observed_Ratio ) ,
      disrate_ncs  %>% select(ryear,'NCS_LBS' = Median.Boot_RETAINED.LBS , "NCS_RATIO" = Observed_Ratio ) , by = 'ryear') %>%
  mutate(tot = CS_LBS+ NCS_LBS) %>%
  group_by(ryear) %>%
  summarise(cs_prop = CS_LBS/tot,
            ncs_prop = NCS_LBS/tot,
            cs_propxrate  = cs_prop*CS_RATIO,
            ncs_propxrate  = ncs_prop*NCS_RATIO,
            total_disrate = cs_propxrate+ncs_propxrate)

discard[discard$yr > 2010,"OBS"] <-discard_late$total_disrate 
discard %>% mutate(fleet = 1, month = 7, OBS = round(OBS,6), sd = round(sd,6),
                   note = paste("#HKL_POT") ) %>% select(yr, month, fleet, OBS,sd, note) #%>%
  # write.table(.,file = paste0("./generated/discard_rates_SS_",Sys.Date(),".csv"),row.names = FALSE, sep = " ", quote = F) ## input line 596
with(discard, plot(yr,OBS,type = 'p', ylim = c(0,0.5))) 
with(base$discard[base$discard$Fleet < 3,], points(Yr, Obs, add = T, col = 'red'))

## Get lengths from discard. Going to use the "weighted" columns
require(tidyr); require(data.table)
ldis <- read.csv("./sablefish_discard_lengths.csv")
ndis  <- read.csv("./sablefish_discard_nuniquetrips.csv") ## use N_unique_trips for nsamp
names(ldis)[1] <- 'Year'
names(ndis)[1] <- 'Gear'

ldis.temp <- ldis %>% filter(Gear == 'HookAndPot') %>% 
  select(Year, Prop.wghtd,Lenbin)  %>% 
  pivot_wider(id_cols = c(Year), names_from = c(Lenbin), values_from= Prop.wghtd) %>%
  mutate(X18 = 0) %>% ## ad empty 18 bin
  select(Year, X18, everything())

names(ldis.temp)[2:ncol(ldis.temp)] <- paste0('X',seq(18,90,2))
ldis.temp2 <- cbind(ldis.temp, ldis.temp[,-1])
ldis.temp2[,-1] <- round(ldis.temp2[,-1],6)
names(ldis.temp2[39:ncol(ldis.temp2)]) <- paste0(names(ldis.temp2)[2:38],".1")

Nsamp <- ndis %>% filter(Gear == 'HookAndPot') %>% select(N_unique_Trips)
cbind(ldis.temp2,Nsamp) %>% data.frame() %>%
  ## sex combined -- just duplicated as in examples
  mutate(month = 7, fleet = 1, sex = 0, part = 1, Nsamp = N_unique_Trips ) %>%
  select(-N_unique_Trips) %>%
  select(Year, month, fleet, sex, part, Nsamp,  X18, everything()) %>%

  write.table(.,
              file = paste0("./generated/lcomp_discard_SS_",Sys.Date(),".csv"),
              row.names = FALSE, sep = " ", quote = F) 

## Comps data reformatting; do all at once ----
## lengths ----
dat = SS_readdat_3.30(file = "C:/Users/mkapur/Dropbox/UW/assessments/sab_2019/200.00_base_files_29May/100.00_base_files/hklpot_agg/base_model/data.ss")
origdat <- dat$lencomp %>% 
  mutate(FltSvy = as.numeric(FltSvy)) %>%
  filter(FltSvy != 1  | FltSvy != 3)  ## we're replacing these
  # mutate(FltSvy = as.numeric(FltSvy) -1) ## move down
# origdat[,1:6] <- lapply(origdat[,1:6],as.numeric)
# origdat[,1:6] <- lapply(origdat[,1:6],round,3)
## new disc comps 
with(subset(lc, FltSvy == 1 & Part == 2), summary(as.numeric(L30)))
with(subset(dat$lencomp, FltSvy == 1 & Part == 2), summary(as.numeric(f30)))
with(subset(dat$lencomp, FltSvy == 2 & Part == 2), summary(as.numeric(f30)))

## old 1 and 2
# lnames = c(paste0("L", seq(18,90,2)),paste0("L", seq(18,90,2),".1")) ## male then female
# lNnames = paste0(lnames,"N")

plot(dat$lencomp$f18[dat$lencomp$FltSvy == 1][1:37] ~ lc$L18)
lc <- read.csv("C:/Users/mkapur/Dropbox/UW/assessments/sab_2019/200.00_base_files_29May/100.00_base_files/hklpot_agg/lcomps_fixtwl_18.csv") %>% filter(FltSvy == 1)
# lc$FltSvy[lc$FltSvy == 2] <- 3
# lc[,lNnames] <- lc[,lnames]*lc[,"Nsamp"] ## back out actual numbers
# 
# ## total times sample sizes -- save these proportions
# lc2 <- lc %>% 
#   select(Yr,lNnames) %>% 
#   group_by(Yr) %>% 
#   summarise_all(funs(sum))
# 
# newNsamp <- lc2 %>% select(-Yr) %>% rowSums(.)
# lc3 <- matrix(NA, nrow = length(unique(lc2$Yr)), ncol = length(lNnames)) %>% data.frame(.)
# names(lc3) <- lNnames
# for(i in 1:length(lNnames)){
#   for(y in 1:length(unique(lc2$Yr))){
#     lc3[y,lNnames[i]] <- round(lc2[y,lNnames[i]] /newNsamp[y],5) ## get new proportions
#   } ## end filling all years (Rows)
# } ## end filling all lengths (cols)
# lc4 <- lc3 %>% mutate(Yr = unique(lc2$Yr), Seas = 7, FltSvy = 1, Gender = 0, Part = 2, Nsamp = round(newNsamp,0)) %>%
#   select(Yr, Seas, FltSvy, Gender, Part, Nsamp, everything()) 
# 
# names(dat$lencomp) <- names(read.csv("./lcomps_fixtwl_18.csv"))
write.table(lc,file = paste0("./generated/lcomp_fixtwl_SS_",Sys.Date(),".csv"),row.names = FALSE, sep = " ", quote = F) ## input line 596

templen = rbind(origdat,read.csv("./lcomps_fixtwl_18.csv"))
unique(templen$FltSvy)
unique(dat$lencomp$FltSvy)
subset(dat$lencomp,Yr == 1986)
subset(templen,Yr == 1986)

## ages ----
# These actually go until 50m kelli sent till 60, truncating for now
# anames = c(paste0("A", seq(0,50,1)),paste0("A", seq(0,50,1),".1")) ## male then female

# aNnames <- paste0(anames,"N")
dat <- SS_readdat_3.30(file = "./base_model/data.ss")
origdat <- dat$agecomp %>% 
  mutate(FltSvy = as.numeric(FltSvy)) %>%
  filter(FltSvy > 3) %>% 
  mutate(FltSvy = as.numeric(FltSvy) -1) ## move down everything above twl
  
origdat[,1:6] <- lapply(origdat[,1:6],as.numeric)
origdat[,1:6] <- lapply(origdat[,1:6],round,3)

ac0 <- read.csv("./acomps_fixtwl_60.csv")
ac0[,60] <- rowSums(ac0[,60:70]) ## sum fem into plus group
ac0[,121] <- rowSums(ac0[,121:131]) ## sum mal into plus group
ac <- ac0[,c(1:60,71:121)] ## drop extra cols
names(ac) <- names(origdat)
newacomps <- rbind(origdat,ac)
newacomps <- newacomps[order(newacomps$FltSvy,newacomps$Yr),]## full dataframe to paste

## sanity check summation happened correctly
as.numeric(newacomps[1,'Nsamp']) == sum(as.numeric(dat$agecomp$Nsamp[ dat$agecomp$Yr == "1986" &  dat$agecomp$FltSvy == "1"]),
                            as.numeric(dat$agecomp$Nsamp[ dat$agecomp$Yr == "1986" &  dat$agecomp$FltSvy == "2"]))
write.table(newacomps,
            file = paste0("./generated/acomp_fixtwl_SS_",Sys.Date(),".csv"),
            row.names = FALSE, sep = " ", quote = F) ## input line 596


## checking something for kelli
## load the dat sent on 05 Jul at 3:56pm
dat <- SS_readdat("C:/Users/mkapur/Dropbox/UW/assessments/sab_2019/data_05Jul356pm.ss")
mhdat <- SS_readdat("C:/Users/mkapur/Dropbox/UW/assessments/sab_2019/200.00_base_files_29May/100.00_base_files/hklpot_agg/201.00_base_files_1July_ToSTAR_FixedGear/data.ss")
## see what the age comp values are
dat$agecomp %>% filter(FltSvy == 1 & Yr == 1991)
mhdat$agecomp %>% filter(FltSvy == 1 & Yr == 1991)

## they should match the summed stuff from what she emailed me
ac %>% filter(FltSvy == 1 & Yr == 1991)



subset(dat$agecomp, Yr == 1986)
subset(newacomps, Yr == 1986)

 # %>%  write.table(.,file = paste0("./generated/acomp_fixtwl_SS_",Sys.Date(),".csv"),row.names = FALSE, sep = " ", quote = F) ## input line 596
base <- SS_output("C:/Users/mkapur/Dropbox/UW/assessments/sab_2019/200.00_base_files_29May/100.00_base_files/hklpot_agg/base_model/", forecast = F)
tmp <- SS_output("C:/Users/mkapur/Dropbox/UW/assessments/sab_2019/200.00_base_files_29May/100.00_base_files/hklpot_agg/100.79/", forecast = F)
tmp$likelihoods_by_fleet
tmp$likelihoods_used
tmp0$likelihoods_by_fleet
tmp$age_comp_fit_table
tmp$agedbase[is.na(tmp$agedbase$effN),]
# tmpdat <- SS_readdat("C:/Users/mkapur/Dropbox/UW/assessments/sab_2019/200.00_base_files_29May/100.00_base_files/hklpot_agg/100.77/data.ss")
# basedat <- SS_readdat("C:/Users/mkapur/Dropbox/UW/assessments/sab_2019/200.00_base_files_29May/100.00_base_files/hklpot_agg/base_model/data.ss")
# 
# SS_readctl_3.30("C:/Users/mkapur/Dropbox/UW/assessments/sab_2019/200.00_base_files_29May/100.00_base_files/hklpot_agg/100.71/control.ss")
# SS_readctl_3.30("C:/Users/mkapur/Dropbox/UW/assessments/sab_2019/200.00_base_files_29May/100.00_base_files/hklpot_agg/base_model/control.ss")


