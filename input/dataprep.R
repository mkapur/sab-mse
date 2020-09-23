## One-timer script to process data in raw_data folder. 
## If applicable, notes are included on manipulations to data that couldn't be done in R or were done elsewhere.
## all processed files are stored in input_data for use in load_data_seasons.R
## on the fly visualizations of input data are also produced here.
## fall 2020 kapurm@uw.edu

require(tidyverse)
require(reshape2)
require(r4ss)
require(here)
require(ggplot2)
require(ggsidekick)
require(PNWColors)

## 2019 Sab WC SS ----
wc <- SS_output(here("input","raw_data","2019 WC Stock Assessment"))

## Fltnames ----
fltnames <- read.table(here("input","input_data","OM_fleetnames.txt"), header = TRUE)
fltnames_fish <- fltnames$NAME[fltnames$COMM]
fltnames_surv <- fltnames$NAME[fltnames$SURV]
fltnames_acomp <- fltnames$NAME[fltnames$ACOMP]
fltnames_lcomp <- fltnames$NAME[fltnames$LCOMP]

nfleets_fish <- length(fltnames$NAME[fltnames$COMM])
nfleets_surv <- length(fltnames$NAME[fltnames$SURV])
nfleets_acomp <- length(fltnames$NAME[fltnames$ACOMP])
nfleets_lcomp <- length(fltnames$NAME[fltnames$LCOMP])
## Color palettes ----

## catch & discards ----

## might need to convert some units
# BC in tonnes (1 tonne = 1 metric ton = 1000kg)
## WC SS in mt (same as BC)
## AK in tons (1000 tons = 1 metric ton)
## format: year x fleet

#* catches ----
omcatch0 <- data.frame(Year = 1960:2018)

## read files
akcatch <- read.csv(here("input","raw_data","catch","AK_comm_catch.csv")) ## this is separated at 145W, corresponding to A1 and A2
names(akcatch)[2:ncol(akcatch)] = rev(paste(fltnames$NAME[fltnames$M == 'AK' & fltnames$COMM]))

bccatch <- read.csv(here("input","raw_data","catch","BC_om_landingCatch.csv"))[,c(1,4,3,5)] ## year, trap trawl and LL
names(bccatch)[2:4] = paste(fltnames$NAME[fltnames$M == 'BC' & fltnames$COMM]  )


wccatch <- wc$catch %>% select(Yr, Obs, Fleet_Name) %>%
  pivot_wider(., names_from = Fleet_Name, values_from = Obs) %>%
  select(-NONE) %>% filter(Yr > 1959)
names(wccatch)[2:3] = paste(fltnames$NAME[fltnames$M == 'WC' & fltnames$COMM]  )

names(bccatch)[1] <- names(wccatch)[1] <- names(akcatch)[1] <- 'Year'

omcatch <- merge(omcatch0, akcatch,  by = "Year", all = TRUE) %>% 
  merge(., bccatch, by = "Year", all = TRUE) %>% 
  merge(., wccatch, by = "Year", all = TRUE) 

write.csv(omcatch %>% select(Year,fltnames_fish) , file = here("input","input_data","OM_catch.csv"), row.names = FALSE)

omcatch <- read.csv(here("input","input_data","OM_catch.csv") )
## plot values
omcatch %>%
  melt(id = "Year") %>%
  ggplot(., aes(x = Year, y = value, color = variable)) +
  theme_sleek() + theme(legend.position = c(0.8,0.8)) +
  scale_color_manual(values = fishfltPal) +
  scale_x_continuous(breaks = seq(1960,2020,10)) +
  geom_line(lwd = 1) +
  labs(x = 'Year', y = 'Obs Catch (mt)', color = 'Comm. Fleet')

ggsave(last_plot(),
       file = here('input','input_data','input_figs','om_catches.png'),
       height = 6, width = 6, unit = 'in', dpi = 420)

#* discards ----
akdis <- read.csv(here("input","raw_data","catch","AK_FIX_discard.csv"),fileEncoding="UTF-8-BOM")[,1:3]
names(akdis)[2:ncol(akdis)] <- paste(fltnames$NAME[fltnames$M == 'AK' & fltnames$DISCARD])

bcdis <- read.csv(here("input","raw_data","catch","BC_om_releasesCatch.csv"))[,c(1,4:5)]
names(bcdis)[2:3] = paste(fltnames$NAME[fltnames$M == 'BC' & fltnames$DISCARD])

wcdis <- wc$discard %>% mutate(Obs = Obs*1000) %>% select(Yr, Obs, Fleet_Name) %>%
  pivot_wider(., names_from = Fleet_Name, values_from = Obs)  %>% 
  filter(Yr > 1959)
names(wcdis)[2:3] <- paste(fltnames$NAME[fltnames$M == 'WC' & fltnames$DISCARD])

names(bcdis)[1] <- names(wcdis)[1] <- names(akdis)[1] <- 'Year'

omdis0 <- data.frame(Year = 1960:2018)
omdis <- merge(omdis0, akdis,  by = "Year", all = TRUE) %>% 
  merge(., bcdis, by = "Year", all = TRUE) %>% 
  merge(., wcdis, by = "Year", all = TRUE) 

omdis[omdis == -1.0] <- NA
save(omdis, file = here("input","input_data","OM_discard.csv"))


## plot values
omdis %>%
  melt(id = "Year") %>%
  ggplot(., aes(x = Year, y = value, color = variable)) +
  theme_sleek() + theme(legend.position = c(0.8,0.8)) +
  scale_color_manual(values = fishfltPal[which(colnames(fishfltPal) %in% names(omdis)[2:ncol(omdis)])]) +
  scale_x_continuous(breaks = seq(1970,2020,10)) +
  geom_line(lwd = 1) +
  labs(x = 'Year', y = 'Obs discard (mt)', color = 'Comm. Fleet',
       subtitle="WC 2019 SS Values multiplied by 1000") 

ggsave(last_plot(),
       file = here('input','input_data','input_figs','om_discards.png'),
       height = 6, width = 6, unit = 'in', dpi = 420)


## length comps  ----
# the dims are YEAR X [AGE or LENGTH] x FLEET
# will RBIND everything with a mutation for fleet for SS use

## ak lcomps
## I manually cut & paste the values into CSV for e/w m/f for each fishery 

## this loop fills in missing length bins (every two for AK)

# for(i in

#){
#   temp <- read.csv(i) %>% 
#     bind_cols(., data.frame(matrix(0, ncol = length(seq(1,40,2)), nrow = 30))) 
#   names(temp) <- c('Year', c(seq(41,99,2),seq(1,40,2)))
#   temp <- temp[,c(1,32:ncol(temp),2:31)] ## reorder
#   write.csv(temp, 
#             file = here('input','raw_data',"comps", paste0(basename(tools::file_path_sans_ext(i)), "_filled.csv")),
#             row.names = FALSE)
#   rm(temp)
# }

#* ak len comps ----
aklc.files <-   list.files(here('input','raw_data',"comps"), pattern = "lcomps", full.names = TRUE)

akfemidx <- c(1,3,9,11,5,7)
akmalidx <- c(1:12)[!c(1:12) %in% akfemidx]
ak_lencomps_female <- ak_lencomps_male <- array(NA, dim = c(length(1960:2018), 
                                                            31,
                                                            length(fltnames$NAME[fltnames$LCOMP & fltnames$M == 'AK']) ))
dimnames(ak_lencomps_female)[[3]] <- dimnames(ak_lencomps_male)[[3]] <- fltnames$NAME[fltnames$LCOMP & fltnames$M == 'AK']

for(i in 1:length(akfemidx)){
  ak_lencomps_female[,,i] <- as.matrix(merge(data.frame('Year' = 1960:2018),
                                             read.csv(aklc.files[akfemidx[i]],fileEncoding="UTF-8-BOM"),
                                             by = 'Year', all.x = TRUE))
}

for(i in 1:length(akmalidx)){
  ak_lencomps_male[,,i] <- as.matrix(merge(data.frame('Year' = 1960:2018),
                                           read.csv(aklc.files[akmalidx[i]],fileEncoding="UTF-8-BOM"),
                                           by = 'Year', all.x = TRUE))
}

save(ak_lencomps_female, file = here('input','raw_data',"comps","ak_lencomps_female_BINDME.rdata"))
save(ak_lencomps_male, file  = here('input','raw_data',"comps","ak_lencomps_male_BINDME.rdata"))




#* bc len comps ----
bcdat <- read.csv(here('input','raw_data',"comps",'BC_LWMSO_1970-present.csv')) ## slow

## calculate raw comp values
bcgears <- c("TRAP","LONGLINE","BOTTOM TRAWL") ## how they appear in data from Brendan
bc_lencomps_female <- bc_lencomps_male <- array(NA, dim = c(length(1978:2018), length(1:88), length(bcgears)),
                                                dimnames = list(paste(1978:2018), paste(1:88),  bcgears))

bc_lencomps_yearN_female <-bc_lencomps_yearN_male <- array(NA, dim = c(length(1978:2018), 2 , length(bcgears)),
                                                           dimnames = list(paste(1978:2018), c('YEAR','n'),  bcgears))

for(i in 1:length(bcgears)){
  bc_lencomps0 <-  bcdat %>%
    filter(GEAR_DESC == bcgears[i] & 
             TRIP_TYPE == 'OBSERVED COMMERCIAL' &
             !is.na(SPECIMEN_AGE) & SPECIMEN_SEX_CODE != 3) %>% 
    select(TRIP_TYPE, GEAR_DESC, NS_AREA, YEAR, SPECIMEN_SEX_CODE, SPECIMEN_AGE) %>%
    group_by(YEAR,SPECIMEN_SEX_CODE,SPECIMEN_AGE) %>%
    summarise(total = n() )%>%   mutate(freq = total/sum(total))  %>%
    tidyr::complete(., SPECIMEN_AGE = 1:88, 
                    SPECIMEN_SEX_CODE = factor(1:2), 
                    YEAR = factor(1978:2018),
                    fill = list(TRIP_TYPE = 'OBSERVED COMMERCIAL',
                                GEAR_DESC ='BOTTOM TRAWL', 
                                total = 0)) 
  
  
  bc_lencomps0$SPECIMEN_AGE <- as.numeric(bc_lencomps0$SPECIMEN_AGE)
  bc_lencomps0$YEAR <- as.numeric(bc_lencomps0$YEAR)
  bc_lencomps0$SPECIMEN_SEX_CODE <- as.factor(bc_lencomps0$SPECIMEN_SEX_CODE)
  
  bc_lencomps0$freq[bc_lencomps0$total==0] <- 0 ## if no real samples, overwrite nan to zero
  
  
  bc_lencomps1 <- bc_lencomps0 %>% distinct() ## rm duplicates created via complete()
  
  ## sanity check - sum freq for every year-sex combo should be 1 or 0 (no samples at all)
  bc_lencomps1 %>% group_by(YEAR, SPECIMEN_SEX_CODE) %>% summarise(sumf = sum(freq))
  bc_lencomps1 %>% filter(YEAR == 1978 & SPECIMEN_SEX_CODE == 1 ) %>% summarise(sum(freq))
  bc_lencomps0 %>% filter(YEAR == 1978 & SPECIMEN_SEX_CODE == 1 ) 
  
  ## should be in the low 100s
  bc_lencomps_yearN_female[,,i] <- bc_lencomps1 %>%     filter(SPECIMEN_SEX_CODE == 1) %>% group_by(YEAR) %>% summarise(n = sum(total)) %>% as.matrix()
  bc_lencomps_yearN_male[,,i] <- bc_lencomps1 %>%     filter(SPECIMEN_SEX_CODE == 2) %>%group_by(YEAR) %>% summarise(n = sum(total)) %>% as.matrix
  
  
  bc_lencomps_female[,,i] <- bc_lencomps1 %>%
    ungroup() %>%
    filter(SPECIMEN_SEX_CODE == 1) %>%
    select(-SPECIMEN_SEX_CODE, -total) %>% #View()
    tidyr::pivot_wider(., values_from = freq, names_from = SPECIMEN_AGE) %>% 
    select(-YEAR) %>%
    as.matrix()
  
  bc_lencomps_male[,,i] <- bc_lencomps1 %>%
    ungroup() %>%
    filter(SPECIMEN_SEX_CODE == 2) %>%
    select(-SPECIMEN_SEX_CODE, -total) %>% #View()
    tidyr::pivot_wider(., values_from = freq, names_from = SPECIMEN_AGE) %>%
    select(-YEAR) %>%
    as.matrix()
}
dimnames(bc_lencomps_male)[[3]] <- dimnames(bc_lencomps_female)[[3]] <- fltnames$NAME[fltnames$LCOMP & fltnames$M == 'BC'][c(2,1,3)]
save(bc_lencomps_female, file = here('input','raw_data',"comps","bc_lencomps_female_BINDME.rdata"))
save(bc_lencomps_male, file = here('input','raw_data',"comps","bc_lencomps_male_BINDME.rdata"))

#* plot bc lcomps ----
# plist = list(); idx = 1
# for(i in 1:length(bcgears)){
#   # par(mfrow = c(3, ceiling(sum(bc_lencomps_yearN_female[,2,i] != 0)/3)))
#   for(y in 1:length(1978:2018)){
#     if(all(bc_lencomps_female[y,,i] == 0) & all(bc_lencomps_male[y,,i] == 0) ) next()
#     plist[[idx]] <- bc_lencomps_female[y,,i] %>%
#       melt() %>%
#       ggplot(., aes(x = 1:88, y = value)) +
#       geom_line(color = sexPal[1], lwd = 1.1, alpha = 0.5) +
#       geom_line(data = melt(bc_lencomps_male[y,,i]),
#                 color = sexPal[2], lwd = 1.1, alpha = 0.5) +
#       ggsidekick::theme_sleek() +
#       labs(x = 'Length (cm)', y = '', title = paste0(1978+y," ", dimnames(bc_lencomps_female)[[3]][i]))
#     idx = idx+1
#   }
# }
# ggsave(Rmisc::multiplot(plotlist = plist,  layout = matrix(1:24, ncol = 4, nrow = 6,byrow = TRUE)),
#        file = here('input','input_data','input_figs','bc_lcomps_female.png'),
#        height = 18, width = 12, unit = 'in', dpi = 420)

#* wc len comps ----
## these are only from NWCBO
wcdat <- SS_readdat(here('input',"raw_data","2019 WC Stock Assessment","data.ss"))
wcLC <- wcdat$lencomp %>% select(-Seas, -Gender, -Part, -Nsamp)

wc_lencomps_female <- bind_cols( Year = wcLC[,1],data.frame(matrix(0, ncol = length(seq(0,16,2)), nrow = 16)),wcLC[,3:39])
wc_lencomps_male <- bind_cols(Year = wcLC[,1],data.frame(matrix(0, ncol = length(seq(0,16,2)), nrow = 16)), wcLC[,40:76])
names(wc_lencomps_female) <- names(wc_lencomps_male) <- c('Year',seq(0,90,2))


save(wc_lencomps_female, file = here('input','raw_data',"comps","wc_lencomps_female_BINDME.rdata"))
save(wc_lencomps_male, file  = here('input','raw_data',"comps","wc_lencomps_male_BINDME.rdata"))

## now bind all and save
OM_lencomps_female <- OM_lencomps_male <- list()

list.files(here('input','raw_data','comps'), pattern = "BINDME", full.names = TRUE) %>%
  lapply(load,.GlobalEnv)

OM_lencomps_male <- list(ak_lencomps_male, bc_lencomps_male, wc_lencomps_male )
save(OM_lencomps_female, file = here('input','input_data',"OM_lencomps_female.rdata"))
save(OM_lencomps_male, file = here('input','input_data',"OM_lencomps_male.rdata"))


## age comps  ----

#* ak agecomps ----
## fixed gear and GOA survey. note they are NOT sex specific.
ak_agecomps  <- array(NA, dim = c(length(1960:2018),
                                  length(0:70),
                                  length(fltnames$NAME[fltnames$ACOMP & fltnames$M == 'AK']) ))
dimnames(ak_agecomps) <-  list(c(1960:2018),
                               c(0:70),
                               c(paste(fltnames$NAME[fltnames$M == 'AK' & fltnames$ACOMP])))

## fix e, fixw, goasurv
ak_agecomps[,c(1,3:31),1] <- as.matrix(merge(data.frame('Year' = 1960:2018),
                                    read.csv(here("input","raw_data","comps","AK_fishery_fixedgear_E_agecomp.csv")),
                                    by = 'Year',all.x = TRUE)%>% select(-Year))
ak_agecomps[,c(1,3:31),2] <- as.matrix(merge(data.frame('Year' = 1960:2018),
                                    read.csv(here("input","raw_data","comps","AK_fishery_fixedgear_W_agecomp.csv")),
                                    by = 'Year',all.x = TRUE) %>% select(-Year))

# ak_agecomps[,,3] <- as.matrix(merge(data.frame('Year' = 1960:2018),
#                                     read.csv(here("input","raw_data","comps","AK_fishery_fixedgear_E_agecomp.csv")),
#                                     by = 'Year',all.x = TRUE))
#Age.pop is the estimate of numbers of fish after running through the age length key. 
# -9 is a catch all for unassigned fish to round out the estimated abundance in that year. So your age comps
# Age, Sex, Age.Pop/sum(Age.Pop) for that sex and year combination. 
## for formatting purposes these are presnetly females only.
actemp <- merge(read.csv(here("input","raw_data","comps","GOA Age Composition Totals.csv")) %>%
                  select(Survey, Year, Sex, Age..years., Age.Pop),read.csv(here("input","raw_data","comps","GOA Age Composition Totals.csv")) %>%
                  select(Survey, Year, Sex, Age..years., Age.Pop) %>%
                  group_by(Survey, Year, Sex) %>%
                  dplyr::summarise(sumap = sum(Age.Pop)),
                by = c('Year','Sex'), all = TRUE) %>%
  mutate(value = Age.Pop/sumap) %>%
  select(Year, Sex, Age..years., value) %>%
  filter(Sex == 'Female' & Age..years. != -9) %>%
  tidyr::pivot_wider(., names_from= Age..years., values_from = value) %>%
  select(-Sex) %>%
  mutate(`19` = NA) %>%
  dplyr::relocate(`1`,.after = Year) %>%
  dplyr::relocate(`19`,.before= `20`) %>%
  dplyr::relocate(`18`,.before = `19`)

ak_agecomps[,2:21,3] <- as.matrix(merge(data.frame('Year' = 1960:2018),actemp,
                                        by = 'Year', all = TRUE) %>% select(-Year))
                                        

  

# read.csv(here("input","raw_data","comps","GOA Age Composition Totals.csv")) %>% View()
# select(Survey, Year,Sex, Age..years.) %>%
#   group_by(Survey, Year, Sex, Age..years.) %>%
#   dplyr::summarise(n = n()) %>%
#   merge(.,goa_agecomps_Nsamp, by = c('Survey','Year','Sex') ) %>%
#   mutate(freq = n/yearN) %>% 
#   tail()
# ungroup() %>%
# group_by(Survey, Year, Sex) %>%
# dplyr::summarise(yearN = n()) 

#* bc agecomps ----
## bc has comps from two surveys but they aren't otherwise used
## deleted the duplicate years on SS.

# read.csv(here("input","raw_data","comps",
#               "BC_om_FemaleSSAgeProp.csv"))[,2:35]%>%
#   filter_all(., any_vars(. != -1)) %>% 
#   mutate('Year' = 1990:2009) %>%
#   select(Year, everything()) %>%
#   write.csv(here("input","raw_data","comps",
#                  "BC_om_FemaleSSAgeProp.csv"), row.names = FALSE) 
# read.csv(here("input","raw_data","comps",
#               "BC_om_MaleSSAgeProp.csv"))[,2:35]%>%
#   filter_all(., any_vars(. != -1)) %>% 
#   mutate('Year' = 1990:2009) %>%
#   select(Year, everything()) %>%
#   write.csv(here("input","raw_data","comps",
#                  "BC_om_MaleSSAgeProp.csv"), row.names = FALSE) 


bc_agecomps_female <- bc_agecomps_male <- array(NA, dim= c(length(1960:2018), 
                                                           35,
                                                           3)) ## SS, commercial, strs

bc_agecomps_female[,,1] <- as.matrix(merge(data.frame('Year' = 1960:2018),
                                           read.csv(here("input","raw_data","comps",
                                                         "BC_om_FemaleCommercialTrapAgeProp.csv")),
                                           by = 'Year', all.x = TRUE))
bc_agecomps_female[,,2] <- as.matrix(merge(data.frame('Year' = 1960:2018),
                                           read.csv(here("input","raw_data","comps",
                                                         "BC_om_FemaleStRSAgeProp.csv")),
                                           by = 'Year', all.x = TRUE))
bc_agecomps_female[,,3] <-  as.matrix(merge(data.frame('Year' = 1960:2018),
                                            read.csv(here("input","raw_data","comps",
                                                          "BC_om_FemaleSSAgeProp.csv")),
                                            by = 'Year', all.x = TRUE))

bc_agecomps_male[,,1] <- as.matrix(merge(data.frame('Year' = 1960:2018),
                                         read.csv(here("input","raw_data","comps",
                                                       "BC_om_MaleCommercialTrapAgeProp.csv")),
                                         by = 'Year', all.x = TRUE))
bc_agecomps_male[,,2] <- as.matrix(merge(data.frame('Year' = 1960:2018),
                                         read.csv(here("input","raw_data","comps",
                                                       "BC_om_MaleStRSAgeProp.csv")),
                                         by = 'Year', all.x = TRUE))
bc_agecomps_male[,,3] <-  as.matrix(merge(data.frame('Year' = 1960:2018),
                                          read.csv(here("input","raw_data","comps",
                                                        "BC_om_MaleSSAgeProp.csv")),
                                          by = 'Year', all.x = TRUE))

dimnames(bc_agecomps_female)[[3]] <- 
  dimnames(bc_agecomps_male)[[3]] <- paste(fltnames$NAME[fltnames$M == 'BC' & fltnames$ACOMP])

save(bc_agecomps_female, file  = here('input','raw_data',"comps","bc_agecomps_female_BINDME.rdata"))
save(bc_agecomps_male, file  = here('input','raw_data',"comps","bc_agecomps_male_BINDME.rdata"))


#* wc agecomps ----
wcdat <- SS_readdat(here('input',"raw_data","2019 WC Stock Assessment","data.ss"))

wc_agecomps_female <- wc_agecomps_male <- array(NA, dim= c(length(1960:2018), 
                                                           length(0:50)+1,
                                                           2),
                                                           #length( paste(fltnames$NAME[fltnames$M == 'WC' & 
                                                            #                             fltnames$ACOMP]))),
                                                dimnames =list( 1960:2018,
                                                                c('Year',0:50), 
                                                                paste(fltnames$NAME[fltnames$M == 'WC' & 
                                                                                      fltnames$ACOMP])[2:3])) 
## only compare fisheries and NWCBO
wcAC_fem <- wcdat$agecomp[,1:60] %>% 
  select(-Seas, -Gender, -Part, -Nsamp, -Ageerr,-Lbin_lo,-Lbin_hi) %>% 
  filter(FltSvy %in% c(1,3)) %>%
  mutate(Year = Yr) %>%
  select(-Yr)

## note that these are CAALs for FlT8
wcAC_mal <- wcdat$agecomp[,c(1:9,61:ncol(wcdat$agecomp))] %>% 
  select(-Seas, -Gender, -Part, -Nsamp, -Ageerr,-Lbin_lo,-Lbin_hi) %>%
  filter(FltSvy %in% c(1,3)) %>%
  mutate(Year = Yr) %>%
  select(-Yr) 

for(flt in 1:length(c(1,3))){
  wc_agecomps_female[,,flt] <- as.matrix(merge(data.frame('Year' = 1960:2018),
                                               wcAC_fem %>% 
                                                 filter(FltSvy == c(1,3,8)[flt]) %>%
                                                 select(-FltSvy),
                                               by = 'Year', all = TRUE))
  wc_agecomps_male[,,flt] <- as.matrix(merge(data.frame('Year' = 1960:2018),
                                             wcAC_mal %>% 
                                               filter(FltSvy == c(1,3,8)[flt]) %>%
                                               select(-FltSvy),
                                             by = 'Year', all = TRUE))
  
}
## deal with CAALs flt 8; Nsamp drama, maybe later

# wc_agecomps_male[,,3] 
# wcAC_fem8 <- wcdat$agecomp[,1:60] %>%  
#   filter(FltSvy ==8 & Gender == 1 ) %>%  
#   mutate(Year = Yr) %>%
#   select(-Seas, -Yr, -Gender, -Part, -Nsamp, -Ageerr, -Lbin_lo, -Lbin_hi,-FltSvy) %>%
#   group_by(Year) %>%
#   summarise_all(funs(sum)) #%>%
# 
# 
# wcAC_mal8 <- wcdat$agecomp[,c(1:9,61:ncol(wcdat$agecomp))] %>%  
#   filter(FltSvy ==8 & Gender == 2 ) %>%  
#   mutate(Year = Yr) %>%
#   select(-Seas, -Yr, -Gender, -Part, -Nsamp, -Ageerr, -Lbin_lo, -Lbin_hi,-FltSvy) %>%
#   group_by(Year) %>%
#   summarise_all(funs(sum)) #%>%
# # group_by(Year) %>%
# # mutate_all(.,  ~(scales::rescale(.) %>% as.vector))
# for(y in 1:nrow(wcAC_fem8)){
#   for(a in 2:ncol(wcAC_fem8)){
#     wcAC_fem8[y,a] <- wcAC_fem8[y,a]/rowSums(wcAC_fem8[y,2:ncol(wcAC_fem8)])
#     wcAC_mal8[y,a] <- wcAC_mal8[y,a]/rowSums(wcAC_mal8[y,2:ncol(wcAC_mal8)])
#   }
# }
# wc_agecomps_female[,,3] <- as.matrix(merge(data.frame('Year' = 1960:2018),
#                                            wcAC_fem8,
#                                            by = 'Year',all = TRUE))
# 
# wc_agecomps_male[,,3] <- as.matrix(merge(data.frame('Year' = 1960:2018),
#                                          wcAC_mal8,
#                                          by = 'Year',all = TRUE))

# save(wc_agecomps_female, file = here('input','raw_data',"comps","wc_agecomps_female_BINDME.rdata"))
# save(wc_agecomps_male, file  = here('input','raw_data',"comps","wc_agecomps_male_BINDME.rdata"))

OM_agecomps_female <- list(ak_agecomps, bc_agecomps_female, wc_agecomps_female )
OM_agecomps_male <- list(ak_agecomps, bc_agecomps_male, wc_agecomps_male )
save(OM_agecomps_female, file = here('input','input_data',"OM_agecomps_female.rdata"))
save(OM_agecomps_male, file = here('input','input_data',"OM_agecomps_male.rdata"))


## mortality ----
## keeping this simple, just one value for each R
M_k <- c(0.7,rep(0.2,3))
save(M_k, file = here('input','input_data','M_k.rdata'))
## growth ----
## load growth params - this used to happen in load_data_seasons
growPar <- read.csv(here("input","raw_data","demography","Table3_2020-05-06phase2.csv"))
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}
growPar$Sex <- substrRight(as.character(growPar$Sex), 1)## overwrite for loops

# array of year x stock x sex
Linf_yk <- L1_yk <- kappa_yk <- sigmaG_yk <- array(NA, 
                                                   dim = c(length(1960:2018),
                                                           length(unique(growPar$Region)),2))

for(s in 1:2){
  for(r in 1:4){
    temp <- subset(growPar, Sex == c("F","M")[s] & Region == paste0("R",c(1,2,4,5)[r]))
    if("pool" %in% temp$Period){
      Linf_yk[,r,s] <- temp$Linf
      L1_yk[,r,s] <- temp$L1
      kappa_yk[,r,s] <- temp$k
      sigmaG_yk[,r,s] <- temp$Sigma
      next()
    }
    Linf_yk[1:(2009-1960),r,s] <- temp$Linf[temp$Period == 'early']
    Linf_yk[(2009-1960):nrow(Linf_yk),r,s] <- temp$Linf[temp$Period == 'late']
    
    L1_yk[1:(2009-1960),r,s] <- temp$L1[temp$Period == 'early']
    L1_yk[(2009-1960):nrow(L1_yk),r,s] <- temp$L1[temp$Period == 'late']
    
    kappa_yk[1:(2009-1960),r,s] <- temp$k[temp$Period == 'early']
    kappa_yk[(2009-1960):nrow(Linf_yk),r,s] <- temp$k[temp$Period == 'late']
    
    sigmaG_yk[1:(2009-1960),r,s] <- temp$Sigma[temp$Period == 'early']
    sigmaG_yk[(2009-1960):nrow(Linf_yk),r,s] <- temp$Sigma[temp$Period == 'late']
  }
}
growthPars <- list("Linf_yk"=Linf_yk, "L1_yk"=L1_yk,
             "kappa_yk"=kappa_yk, "sigmaG_yk"=sigmaG_yk)
save(growthPars, 
     file = here('input','input_data',"OM_growthPars.rdata"))


## weight at length
wtatlen_kab <- matrix(NA, nrow = 4, ncol = 2) ## stock x a,b for aL^b
rownames(wtatlen_kab) <- paste0("R",4:1)

wtatlen_kab[4,1] <- wc$Growth_Parameters$WtLen1
wtatlen_kab[4,2] <- wc$Growth_Parameters$WtLen2

## Hybrid of WC and BC, from Cox 2011
wtatlen_kab[3,1] <- sum(8.58e-6,wc$Growth_Parameters$WtLen1)/2
wtatlen_kab[3,2] <- sum(3.05,wc$Growth_Parameters$WtLen2)/2

## mean BC and inflated BC used for AK-BC
wtatlen_kab[2,1] <- sum(8.58e-6,8.58e-6 *1.15)/2
wtatlen_kab[2,2] <-  sum(3.05,3.05*1.15)/2
# wtatlen_kab[2,2] <-  4

## AK uses wt at age. Let's just increase the BC values a bit
wtatlen_kab[1,1] <- 8.58e-6
wtatlen_kab[1,2] <- 3.05*1.15

save(wtatlen_kab, 
     file = here('input','input_data',"OM_wtatlen_kab.rdata"))



# Linf_yk[1:(2009-1966),1:nstocks] <-  growPar$Linf[growPar$Period == 'early'][1:nstocks]
# Linf_yk[(2009-1966):nrow(Linf_yk),1:nstocks] <-  growPar$Linf[growPar$Period == 'late'][1:nstocks]

## maturity ----
## Ben sent me wiggly curves and p25/50/75 values for FEMALES for the defined regions
## note that we use log15 for the 75%ile
logistic <- function(age, a50, a75){
  selage <- 1/(1+exp(-(log(15))*(age-a50)/(a75-a50)))
  return(selage)
}

matPar <- data.frame('region' = paste0('R',1:4),
                     'a50' = c(8.3,4.5,4.6,6.8),
                     'a75' = c(8.7,5.1,4.97,7.7))

mat_ak <- matrix(NA, nrow = length(0:70), ncol = 4)

for(a in 1:nrow(mat_ak)){
  for(k in 1:nrow(matPar)){
    mat_ak[a,k] <- logistic(age = a, a50 = matPar$a50[k], a75 = matPar$a75[k])
  }
}
save(mat_ak, file =  here('input','input_data',"OM_maturity_ak.rdata"))

mat_ak %>%
  data.frame() %>%
  mutate(age = 0:70) %>%
  melt(id = 'age') %>%
  ggplot(., aes(x = age, y = value, color = variable)) +
  geom_line(lwd = 1.1) +
  scale_color_manual(values = rev(demPal), labels =  paste0('R',1:4)) +
  theme_sleek() +theme(legend.position = c(0.7,0.7)) +
  labs(x = 'Age', y = 'Proportion Mature', color = 'Stock')

ggsave(last_plot(),
       file = here('input','input_data','input_figs','om_maturity.png'),
       height = 3, width = 3.5, unit = 'in', dpi = 420)

## movement ----
## see raw_data/demography/movement/movement-estimates-toAge.R
## the mod needs to be run first to come up with the conversion factors --


## selex [fishery] ----
## array year x age x fleet for each sex
## these are really ballparks and will need to be tweaked to condition OM

OM_surv_selex_yafs <- array(NA, dim = c(length(1960:2018),
                                        length(0:70),
                                        length(paste(fltnames$NAME[fltnames$SURV])),
                                        2))
OM_fish_selex_yafs <- array(NA, dim = c(length(1960:2018),
                                        length(0:70),
                                        length(paste(fltnames$NAME[fltnames$COMM])),
                                        2))
dimnames(OM_surv_selex_yafs) =  list( 1960:2018,   c(0:70),  paste(fltnames$NAME[fltnames$SURV]))
dimnames(OM_fish_selex_yafs) =  list( 1960:2018,   c(0:70),  paste(fltnames$NAME[fltnames$COMM]))


#* ak aselex ----
## AK values from an email with Kari on 25 Aug 2020 from tem.par
logistic2 <- function(age, a50, delta){
  selage <- 1/(1+exp(-delta*(age-a50)))
  return(selage)
}
## Longline/fixed gear fishery for post-IFQ years (1995 onward), log scale values
## right now assuming TRAWL has the same
a50_fem_ak <- 1.07528048238
del_fem_ak <- 0.764917256584
a50_mal_ak <- 1.14893337353
del_mal_ak <- 0.945517826919

for(a in 1:nage){
  OM_fish_selex_yafs[,a,1:4,1] <-  logistic2(age = a, a50 = a50_fem_ak, delta = del_fem_ak)
  OM_fish_selex_yafs[,a,1:4,2] <- logistic2(age = a, a50 = a50_mal_ak, delta = del_mal_ak)
}


#* bc lselex ----
## BC only has LEN selex

## code from brendan
# Code below should give you everything you need to estimate selectivity at length for each gear/survey type.
# Trap (g = 1) and longline (g = 2) are dome-shaped using a normal distribution with mean = alpha and SD = beta (selType = 2)
# Trawl (g = 3) is dome shaped with gamma distribution with shape= alpha and scale = beta (selType = 3)
# Surveys (Std: g = 4; StRs: g = 5) are asymptotic with a logistic 
# parameterization, with L50 = alpha - beta, and SD = beta (selType = 1)

# Initial values. LL is pos 3, Trap is pos 2, Trawl is pos 4
alpha_g1 <- c(62.8329, 63.6959, 33.8898, 54.1045, 64.2127)
beta_g1 <- c(7.04483, 3.09715, 1.41494, 4.55724, 12.9197)
selType <- c(2,2,3,1,1)
nG <- length(alpha_g1)
# Calculate selectivity by gear and length
len <- 32:75
sel_lg <- array(0, dim = c(length(len),5))
for(g in 1:nG){
  if(selType[g] == 1)
  {
    sel_lg[,g] <- 1 / ( 1 + exp( - log(19) * (len - alpha_g1[g] + beta_g1[g]) / beta_g1[g] ) )
  }
  if(selType[g] == 2)
  {
    sel_lg[,g] <- exp(-(0.5 * (len - alpha_g1[g])/beta_g1[g])^2 )
  }
  if(selType[g] == 3)
  {
    sel_lg[,g] <- len ^(alpha_g1[g] - 1) * exp(-len/beta_g1[g])
    sel_lg[,g] <- sel_lg[,g] / max(sel_lg[,g])
  }
}
selMat <- cbind(len, sel_lg)
colnames(selMat) <- c("Length","Trap","LL","Trawl","Std","StRS")

## find age at length using saming init_LAA key as for movement [dims alis]
# load(here("input","raw_data","demography","movement", "init_LAA.rda")) ## from prelim runs
# B1 = B2 = matrix(NA, nrow = 82, ncol = 3)
# B1[,1] <- B2[,1] <- 0:81
# for(s in 1:2){
#   for(l in 0:81){
#     B1[l,s+1] <-which.max(init_LAA[,l,3,s])
#     B2[l,s+1] <-which.max(init_LAA[,l,4,s])
#   }
# }
# BCALK = data.frame(cbind(B1,B2[,2:3]))
# names(BCALK) = c('len','ageb1F','ageb1M','ageb2F','ageb2M')

## LL, TRAP, TRAWL
for(y in 1:dim(OM_fish_selex_yafs)[[1]]){
  # OM_fish_selex_yafs[y,,5,1:2] <- c(rep(0, length(0:31)),t(selMat[1:39,'Trap'])) ## LL
  # OM_fish_selex_yafs[y,,6,1:2] <- c(rep(0, length(0:31)),t(selMat[1:39,'Trap'])) ## TRAP
  # OM_fish_selex_yafs[y,,7,1:2] <- c(rep(0, length(0:31)),t(selMat[1:39,'Trawl'])) ## TRAWL
  
  
  OM_fish_selex_yafs[y,,5,1:2] <- c(rep(0, length(0:31)),t(selMat[1:39,'Trap'])) ## LL
  OM_fish_selex_yafs[y,,6,1:2] <- c(rep(0, length(0:31)),t(selMat[1:39,'Trap'])) ## TRAP
  OM_fish_selex_yafs[y,,7,1:2] <- c(rep(0, length(0:31)),t(selMat[1:39,'Trawl'])) ## TRAWL
}



logistic3 <- function(age, a50, a95){
  selage <- 1/(1+exp(-log(19)*(age-a50)/(a95-a50)))
  return(selage)
}
## for ease we should use the BC growth curve (from their OM) and convert to exp age
## since this can get tweaked down the road I will eyeball from the B1 
# l50_bc <- 52.976
## assuming trap, trawl and fix have the same selex
# a50_fem_bc <- 17; a95_fem_bc <- 20
# a50_mal_bc <- 15; a95_mal_bc <- 17
# 
# for(a in 1:nage){
#   OM_fish_selex_yafs[,a,5:7,1] <-   logistic3(age = a, a50 = a50_fem_bc, a95 = a95_fem_bc)
#   OM_fish_selex_yafs[,a,5:7,2] <- logistic3(age = a, a50 = a50_mal_bc, a95 = a95_mal_bc)
# }
#* wc aselex ----
## 8 and 9 are WC fix and TWL, which correspond to fleets 1 and 3 in SS
## these entries are timeblocked. For ease I wrote this CSV and dragged the values.
# wcas0 <- wc$ageselex %>%
#   filter(Yr > 1880 & Factor == 'Asel' & Fleet %in% c(1,3)) %>%
#   filter(Sex == 1) %>%
#   mutate(Year = Yr) %>%
#   select(-Factor, -Seas, -Morph,-Label, -Yr) %>%
#   group_by(Fleet) %>%
#   complete(., Year = c(1890,1960:2018)) %>% write.csv(.,here('input','raw_data','selex','wcas0F.csv'))
# wcas0 <- wc$ageselex %>%
#   filter(Yr > 1880 & Factor == 'Asel' & Fleet %in% c(1,3)) %>%
#   filter(Sex == 2) %>%
#   mutate(Year = Yr) %>%
#   select(-Factor, -Seas, -Morph,-Label, -Yr) %>%
#   group_by(Fleet) %>%
#   complete(.,Year = c(1890,1960:2018)) %>% write.csv(.,here('input','raw_data','selex','wcas0M.csv'))


# for(s in 1:2){
#   for(flt in c(1,3)){
#     wcas0[wcas0$Year < 1996 & wcas0$Sex == s & wcas0$Fleet == flt,4:ncol(wcas0)] <-
#       wcas0[wcas0$Year == 1996 & wcas0$Sex == s & wcas0$Fleet == flt,4:ncol(wcas0)]
#     
#     wcas0[wcas0$Year %in% 1997:2002 & wcas0$Sex == s & wcas0$Fleet == flt,4:ncol(wcas0)] <-
#       wcas0[wcas0$Year%in% 1997:2002  & wcas0$Sex == s & wcas0$Fleet == flt,4:ncol(wcas0)]
#     wcas0[wcas0$Year %in% 2003:2017 & wcas0$Sex == s & wcas0$Fleet == flt,4:ncol(wcas0)] <-
#       wcas0[wcas0$Year%in% 2003:2017  & wcas0$Sex == s & wcas0$Fleet == flt,4:ncol(wcas0)]
#     
#     wcas0[wcas0$Year > 2017 & wcas0$Sex == s & wcas0$Fleet == flt,4:ncol(wcas0)] <-
#       wcas0[wcas0$Year > 2017  & wcas0$Sex == s & wcas0$Fleet == flt,4:ncol(wcas0)]
#   }
# }

wcas0 <- read.csv(here('input','raw_data','selex','wcas0.csv'))

for(flt in 1:2){
OM_fish_selex_yafs[,,c(8,9)[flt],1] <-  as.matrix( merge(data.frame('Year' = 1960:2018),
                                                        wcas0 %>%
                   filter(Sex == 1 & Fleet == c(1,3)[flt]) %>%
                   select(-Sex) %>%
                   select(-Fleet),
                 by= 'Year', 
                 all.x = TRUE)
                 %>% select(-Year))

OM_fish_selex_yafs[,,c(8,9)[flt],2] <- as.matrix( merge(data.frame('Year' = 1960:2018),
                                                        wcas0 %>%
                                                          filter(Sex == 2 & Fleet == c(1,3)[flt]) %>%
                                                          select(-Sex) %>%
                                                          select(-Fleet),
                                                 by= 'Year', all.x = TRUE) %>% select(-Year)) 


}

## use NWSLP as placeholder because NWCBO is CAAL
for(y in 1:nyear){
  for(flt in 1:nfleets_surv){
    OM_surv_selex_yafs[y,,flt,1] <-  as.matrix(wc$ageselex %>% filter(Fleet == 8 &
                                                               Factor == 'Asel' & Yr == 2018,
                                                             Sex == 1) %>%
      select(-Factor, -Seas, -Morph,-Label, -Yr,-Fleet,-Sex))
    
    
    OM_surv_selex_yafs[y,,flt,2] <-as.matrix(wc$ageselex %>% filter(Fleet == 8 &
                                                                      Factor == 'Asel' & Yr == 2018,
                                                                    Sex == 2) %>%
                                               select(-Factor, -Seas, -Morph,-Label, -Yr,-Fleet,-Sex))
  } ## end flt
} ## end yr

## selex [survey] ---


# 
# for(flt in 1:3){
#   
#   wc_fem_asel[,,flt] <-   as.matrix( merge(data.frame('Year' = 1960:2018),
#                                            wc$ageselex %>% 
#                                              filter(Yr > 1959 & Yr < 2019 & Factor == 'Asel') %>%
#                                              mutate(Year = Yr) %>%
#                                              select(-Factor, -Seas, -Morph,-Label, -Yr) %>%
#                                              filter(Sex == 1 & Fleet == c(1,3,8)[flt]) %>%
#                                              select(-Sex) %>%
#                                              select(-Fleet),
#                                            by= 'Year', all.x = TRUE) )
#   
#   wc_mal_asel[,,flt] <-    as.matrix( merge(data.frame('Year' = 1960:2018),
#                                             wc$ageselex %>% 
#                                               filter(Yr > 1959 & Yr < 2019 & Factor == 'Asel') %>%
#                                               mutate(Year = Yr) %>%
#                                               select(-Factor, -Seas, -Morph,-Label, -Yr) %>%
#                                               filter(Sex == 2 & Fleet == c(1,3,8)[flt]) %>%
#                                               select(-Sex) %>%
#                                               select(-Fleet),
#                                             by= 'Year', all.x = TRUE) )
#   
# }    


save(OM_fish_selex_yafs, file = here('input','input_data',"OM_fish_selex_yafs.rdata"))
save(OM_surv_selex_yafs, file = here('input','input_data',"OM_surv_selex_yafs.rdata"))

#* plot input selex ----
png(here('input','input_data','input_figs','fishery_selex.png'),
    height = 8, width = 6, unit = 'in', res = 420)
par(mfrow = c(3,3) )
for(flt in 1:nfleets_fish){
  for(s in 1:2){
    tmp <- OM_fish_selex_yafs[59,,flt,s]
    if(s == 1) plot(tmp, 
                    col = sexPal[1], 
                    type = 'l', lwd = 2, 
                    xlab = ifelse(fltnames$SELTYPE[fltnames$COMM][flt] == 'AGE',
                                  'Age','Length'), 
                    ylab = 'Selectivity',
                    lty = 1,
                    ylim = c(0,1), main = fltnames_fish[flt], xlim = c(0,75),
                    col.main  = c(rep(mgmtPal[1],4), rep(mgmtPal[2],3),rep(mgmtPal[3],2))[flt])
    box(which = 'plot', lty = 'solid', 
        col = c(rep(mgmtPal[1],4), rep(mgmtPal[2],3),rep(mgmtPal[3],2))[flt], 
        lwd = 2)
    if(s == 2) lines(tmp, col = sexPal[2], type = 'l', lty = 2, lwd = 2)
  }
}
dev.off()

png(here('input','input_data','input_figs','survey_selex.png'),
    height = 8, width = 6, unit = 'in', res = 420)
par(mfrow = c(2,3) )
for(flt in 1:nfleets_surv){
  for(s in 1:2){
    tmp <- OM_surv_selex_yafs[59,,flt,s]
    if(s == 1) plot(tmp, 
                    col = sexPal[1], 
                    type = 'l', lwd = 2, 
                    xlab = ifelse(fltnames$SELTYPE[fltnames$SURV][flt] == 'AGE',
                                  'Age','Length'), 
                    ylab = 'Selectivity',
                    lty = 1,
                    ylim = c(0,1), main = fltnames_surv[flt], xlim = c(0,75),
                    col.main  = c(rep(mgmtPal[1],4), rep(mgmtPal[2],3),rep(mgmtPal[3],2))[flt])
    box(which = 'plot', lty = 'solid', 
        col = c(rep(mgmtPal[1],4), rep(mgmtPal[2],3),rep(mgmtPal[3],2))[flt], 
        lwd = 2)
    if(s == 2) lines(tmp, col = sexPal[2], type = 'l', lty = 2, lwd = 2)
  }
}
dev.off()
## survey ----
## spatial matrix -- for matching on region, stock, sub_area, etc
spmat <- data.frame(subarea = c('A1',"A2","B2","B1","C2","C1"),
                    stock = c("R4","R3","R3","R2","R2","R1"),
                    mgmt = c("AI","AK", rep("BC",2), rep("CC",2)))
bcnom <- read.csv(here("input","raw_data","survey","BC_early_index.csv")) %>%
  mutate(SE = 0.317, lci = nominal.Trap.CPUE-1.96*SE, uci =nominal.Trap.CPUE+1.96*SE, Fleet = "BC_early") %>%
  select(YEAR, nominal.Trap.CPUE, SE, Fleet)
names(bcnom) <- c('Year','value', 'sigma', 'fleet')
#* survey error ----
## reformat this and save            
survsig <- read.csv(here("input","raw_data","survey","Indices_SS3_2020-01-23v3.csv"))  %>% ## VAST stdization
  distinct(Fleet, Year, Estimate_metric_tons, .keep_all = TRUE) %>% ## remove any dupes
  filter(Fleet != "AllAreas" & Fleet != "Eastern_Bering_Sea") %>%
  # merge(.,spmat, by.x = "Fleet", by.y = "mgmt", all.y = FALSE) %>%
  mutate(value = Estimate_metric_tons,
         sigma = SD_log, 
         fleet = Fleet, 
         mgmt = fleet, 
         type = 'survey') %>%
  select(Year, value, sigma, fleet) %>%
  bind_rows(bcnom) %>%
  select(-value) %>%
  pivot_wider(., id_cols = Year, names_from = fleet, values_from = sigma) %>%
  merge(., data.frame('Year' = 1960:2018), all = TRUE) %>%
  select(-Year) 

  # arrange(.,Year,fleet) #%>%

names(survsig) <- paste(fltnames$NAME[fltnames$SURV][c(5,4,1,2,3)])
write.csv(survsig %>% select(fltnames_surv),here("input","input_data","OM_indices_sigma.csv"),row.names = FALSE)

## make columns as fleets, include extra  bc surv
names(bcnom) <- c('Year','value', 'sigma', 'Fleet')
#* survey biomass ----
vast0 <- read.csv(here("input","raw_data","survey","Indices_SS3_2020-01-23v3.csv"))  %>% ## VAST stdization
  filter(Fleet != "AllAreas" & Fleet != "Eastern_Bering_Sea") %>%
  merge(.,spmat, by.x = "Fleet", by.y = "mgmt", all.y = FALSE) %>%
  distinct(Fleet, Year, Estimate_metric_tons, .keep_all = TRUE) %>% ## remove any dupes
  mutate(value = Estimate_metric_tons,
         sigma = SD_log, 
         Fleet = Fleet, 
         mgmt = Fleet, 
         type = 'survey') %>%
  select(Year, Fleet, type, value, sigma) %>%
  arrange(.,Year,Fleet) 
surv_vals <- vast0 %>%
  select(Year, Fleet, value) %>%
  rbind(., bcnom[,c(1,4,2)]) %>%
  arrange(.,Year,Fleet) %>%
  tidyr::pivot_wider(names_from= Fleet, values_from = value) %>%
  filter(Year > 1964 & Year < 2019) %>%
merge(., data.frame('Year' = 1960:2018), all = TRUE) 
surv_vals[surv_vals == -1] <- NA
names(surv_vals)[2:6] <- paste(fltnames$NAME[fltnames$SURV][c(3,2,1,4,5)]) 

write.csv(surv_vals %>% select(fltnames_surv), here("input","input_data","OM_indices.csv"),row.names = FALSE) ## save in special order

surv_vals %>%
  # select(-BC_EARLY) %>%
  mutate(BC_EARLY = BC_EARLY*1000) %>%
  melt(id = "Year") %>%
  ggplot(., aes(x = Year, y = value, color = variable)) +
  theme_sleek() + theme(legend.position = c(0.8,0.8)) +
  scale_color_manual(values = survfltPal)+
  scale_x_continuous(breaks = seq(1970,2020,10)) +
  geom_line(lwd = 1) +
  labs(x = 'Year', y = 'Index of Relative Abundance', color = 'Survey Fleet') +
  labs(subtitle = "BC_EARLY has been multiplied by 1000 for comparison")
ggsave(last_plot(),
       file = here('input','input_data','input_figs','OM_indices.png'),
       height = 6, width = 6, unit = 'in', dpi = 420)

## Aging error ----
## from MH on google drive; just use the "first" for each mgmt region
## array of age and sd x ages
ageerr_SD <- list.files(here("input","raw_data","survey"), pattern = "*SS_format*", full.names = TRUE) %>%
  lapply(read.csv) %>%
  bind_rows() %>% 
  filter(X == 'SD')

ageerr_ExpAge <- list.files(here("input","raw_data","survey"), pattern = "*SS_format*", full.names = TRUE) %>%
  lapply(read.csv) %>%
  bind_rows() %>% 
  filter(X == 'Expected_age')

rownames(ageerr_SD) <- rownames(ageerr_ExpAge) <- c('AK','BC','WC')

save(ageerr_SD, file = here("input","input_data","ageerr_SD.rdata"))
save(ageerr_ExpAge, file = here("input","input_data","ageerr_ExpAge.rdata"))

# Load the age comps 
# age_survey.tmp <- read.csv(here("input","data",'age_survey_ss.csv'))
# age_survey.tmp2 <- array(NA, dim = c(nrow(age_survey.tmp),ncol(age_survey.tmp), nfleets_surv)) ## placeholder; the last term should be nfleets-acomp
# age_survey.tmp2[,,1] <- age_survey.tmp2[,,2] <- as.matrix(age_survey.tmp)
# age_catch.tmp <- read.csv(here("input","data",'age_catch_ss.csv'))
# ac.data <- read.csv(here("input","data",'ac_data.csv'))




