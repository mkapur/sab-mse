## Dataprep
## M S KAPUR
## wrangling of raw synthesis/catch/misc stuff received from management regions
## this does NOT get it in TMB ready format, but writes the CSVs that serve as the baseline
## and therefore only need to get updated if we've received new data from a management region


## spatial matrix -- for matching on region, stock, sub_area, etc
spmat <- data.frame(subarea = c('A1',"A2","B3","B2","B1","C2","C1"),
                    stock = c("R5","R4","R4","R3","R2","R2","R1"),
                    mgmt = c("AI","AK", rep("BC",3), rep("CC",2)))



# Landings & Discards ----
ccbase <- SS_output("./input/raw/CC_2019_100.00/", forecast = FALSE, covar = FALSE) ## from STAR - possibly not the most updated


## METHOD 1 [me]
cc_catdis <- ccbase$catch %>%
  select(Fleet_Name, Yr, Obs, ret_bio) %>%
  mutate(discard = Obs-ret_bio) %>%
  select(-ret_bio) %>%
  reshape2::melt(id = c("Yr", "Fleet_Name")) %>%
  mutate(Year = Yr,
         fleet = Fleet_Name,
         type = ifelse(as.character(variable) == 'Obs','Catch','Discard'),
         stock = NA,
         subarea = NA,
         mgmt = "CC") %>%
  select(Year, fleet, type, subarea, stock, mgmt,  everything(), -Yr, -Fleet_Name, -variable) 
cc_catdis %>%
  filter(type == 'Catch') %>%
  write.csv(., file = paste0("input/temp/cc_catch_",Sys.Date(),".csv"), row.names = FALSE)

## METHOD 2 [SAM]
cn <- ccbase$catch %>%
  select(Fleet_Name, Yr, Obs, ret_bio) %>%
  mutate(discard = Obs-ret_bio) %>%
  select(Fleet_Name, Yr, Obs) %>%  
  # tidyr::pivot_wider(names_from = Fleet_Name,  values_from = Obs) %>%
  # group_split(Fleet_Name) %>% data.frame() %>% 
  #   tibble::column_to_rownames('Yr') %>%
  # select(Fleet_Name, Obs) %>%
  group_split(Fleet_Name) 


rownames(cn) <- cn$Yr
save(], file = paste0("input/temp/cc_catch_",Sys.Date(),".Rdata"))

cc_catdis0
cc_catdis %>%
  filter(type == 'Catch') %>%
  write.csv(., file = paste0("input/temp/cc_catch_",Sys.Date(),".csv"), row.names = FALSE)



cc_catdis %>%
  filter(type == 'Discard') %>%
  write.csv(., file = paste0("input/temp/cc_discard_",Sys.Date(),".csv"), row.names = FALSE)


## ALASKA ##
akcat <- read.csv("./input/raw/AK_catch_ByGear_and_ByArea.csv")
akdis <- read.csv("./input/raw/AK_discards.csv")

ak_catdis <- data.frame(Year = akdis$Year, Gear = akdis$Gear, akdis[grep(paste("BSAI*","GOA*",sep= "|"), names(akdis))]) %>%
  select(-BSAI_PctDiscard, -GOA_PctDiscard) %>%
  reshape2::melt(id = c("Year","Gear")) %>%
  mutate(variable = as.character(variable)) %>%
  mutate(type = stringr::str_split_fixed(variable, "_", 2)[,2],
         stock_temp =  stringr::str_split_fixed(variable, "_", 2)[,1]) %>%
  
  mutate(fleet = Gear,
         stock = ifelse(stock_temp == 'BSAI', "R5","R4")) %>%
  merge(.,spmat, by.x = "stock", by.y = "stock", all.y = FALSE) %>%
  select(Year,fleet, type, subarea, stock, mgmt,  everything(), -Gear, -stock_temp, -variable) 

ak_catdis %>%
  filter(type == 'Catch') %>%
  write.csv(., file = paste0("input/temp/ak_catch_",Sys.Date(),".csv"), row.names = FALSE)

ak_catdis %>%
  filter(type == 'Discard') %>%
  write.csv(., file = paste0("input/temp/ak_discard_",Sys.Date(),".csv"), row.names = FALSE)


## combine landings/discards and write to CLEAN
list.files("./input/temp/", full.names = TRUE)[grep("_catch_", list.files("./input/temp/", full.names = TRUE))] %>%
  lapply(.,read.csv) %>% bind_rows() %>%
  write.csv(., "./input/cleaned/clean_landings.csv",row.names = FALSE)

list.files("./input/temp/", full.names = TRUE)[grep("_discard_", list.files("./input/temp/", full.names = TRUE))] %>%
  lapply(.,read.csv) %>% bind_rows() %>%
  write.csv(., "./input/cleaned/clean_discard.csv",row.names = FALSE)



## Surveys - from VAST ----
read.csv("./input/raw/Indices_SS3_2020-01-23v3.csv")  %>% ## VAST stdization
  distinct(Fleet, Year, Estimate_metric_tons, .keep_all = TRUE) %>% ## remove any dupes
  filter(Fleet != "AllAreas") %>%
  merge(.,spmat, by.x = "Fleet", by.y = "mgmt", all.y = FALSE) %>%
  mutate(value = Estimate_metric_tons,
         sigma = SD_log, 
         fleet = Fleet, 
         mgmt = fleet, 
         type = 'survey') %>%
  select(Year, fleet, type, subarea, stock, mgmt, value, sigma) %>%
  arrange(.,Year,fleet) %>%
  write.csv(., "./input/cleaned/clean_survey.csv",row.names = FALSE)

## make columns as fleets, include extra  bc surv
bcnom <- read.csv('./input/raw/om_indexSeries.csv') %>%
  mutate(SE = 0.317, lci = nominal.Trap.CPUE-1.96*SE, uci =nominal.Trap.CPUE+1.96*SE, Fleet = "BC_E") %>%
  select(YEAR, nominal.Trap.CPUE, SE, Fleet)
names(bcnom) <- c('Year','value', 'sigma', 'Fleet')

vast0 <- read.csv("./input/raw/Indices_SS3_2020-01-23v3.csv")  %>% ## VAST stdization
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
  filter(Year > 1978 & Year < 2019)
surv_vals[is.na(surv_vals)] <- -1
# row.names(surv_vals) <- surv_vals$Year
write.csv(surv_vals[,c(3,4,2,5,6)], "./input/cleaned/clean_survey.csv",row.names = FALSE) ## save in special order

modyr

## survey - SAM format
surveys <- read.csv("./input/raw/Indices_SS3_2020-01-23v3.csv")  %>% ## VAST stdization
  distinct(Fleet, Year, Estimate_metric_tons, .keep_all = TRUE) %>% ## remove any dupes
  filter(Fleet != "AllAreas") %>%
  merge(.,spmat, by.x = "Fleet", by.y = "mgmt", all.y = FALSE) %>%
  mutate(value = Estimate_metric_tons,
         sigma = SD_log, 
         fleet = Fleet, 
         mgmt = fleet, 
         type = 'survey') %>%
  select(Year, fleet, type, subarea, stock, mgmt, value, sigma) %>%
  arrange(.,Year,fleet) %>%
  filter( fleet == 'CC') %>%
  group_split(fleet) 
save(surveys, file = paste0("input/temp/cc_surv_",Sys.Date(),".Rdata"))


## Comps ----

bcdat <- read.csv(here('input','raw','BC_LWMSO_1970-present.csv')) ## slow

## calculate raw comp values
bcgears <- c("TRAP","LONGLINE","BOTTOM TRAWL")
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

plist = list(); idx = 1

for(i in 1:length(bcgears)){
  # par(mfrow = c(3, ceiling(sum(bc_lencomps_yearN_female[,2,i] != 0)/3)))
  for(y in 1:length(1978:2018)){
    if(all(bc_lencomps_female[y,,i] == 0) & all(bc_lencomps_male[y,,i] == 0) ) next()
    plist[[idx]] <- bc_lencomps_female[y,,i] %>%
      melt() %>%
    ggplot(., aes(x = 1:88, y = value)) +
      geom_line(color = 'red', lwd = 1.1, alpha = 0.5) +
      geom_line(data = melt(bc_lencomps_male[y,,i]), 
                color = 'blue', lwd = 1.1, alpha = 0.5) +
      ggsidekick::theme_sleek() +
      labs(x = 'Length CM', title = paste0(1978+y," ",bcgears[i]))
    idx = idx+1
  }
}


Rmisc::multiplot(plotlist = plist, cols = 4)
save(bc_lencomps_female, here('input',"cleaned","bc_Lcomps_FEMALE.rdata"))
save(bc_lencomps_male, paste(here('input',"cleaned","bc_Lcomps_MALE.rdata")))

ccdat <- SS_readdat(file = "./input/raw/CC_2019_100.00/data.ss")
ccdat$agecomp$fleet <-   ccdat$fleetnames[ccdat$agecomp$FltSvy ] ## rename fleets
ccdat$lencomp$fleet <-   ccdat$fleetnames[ccdat$lencomp$FltSvy ] ## rename fleets

## METHOD 1
ccdat$agecomp %>%
  mutate(Year = Yr, subarea = NA, stock = NA, mgmt = "CC", type = 'age_comp' ) %>%
  select(-Seas, - FltSvy, -Gender, -Part, -Ageerr, -Lbin_lo, -Yr, -Lbin_hi) %>%
  select(Year, fleet, type, subarea,stock, mgmt, Nsamp, everything()) %>%
  write.csv(., "./input/temp/CC_agecomps.csv",row.names = FALSE)

ccdat$lencomp %>%
  mutate(Year = Yr, subarea = NA, stock = NA, mgmt = "CC", type = 'age_comp' ) %>%
  select(-Seas, - FltSvy, -Gender, -Part, -Yr) %>%
  select(Year, fleet, type, subarea,stock, mgmt, Nsamp, everything()) %>%
  write.csv(., "./input/temp/CC_lencomps.csv",row.names = FALSE)


## METHOD 2
ccdat$agecomp %>%
  mutate(Year = Yr, subarea = NA, stock = NA, mgmt = "CC", type = 'age_comp' ) %>%
  select(-Seas, - FltSvy, -Gender, -Part, -Ageerr, -Lbin_lo, -Yr, -Lbin_hi) %>%
  select(Year, fleet, Nsamp, everything()) %>%
  group_split(fleet)


write.csv(., "./input/temp/CC_agecomps.csv",row.names = FALSE)


## combine comps and write to CLEAN
list.files("./input/temp/", full.names = TRUE)[grep("_lencomps", list.files("./input/temp/", full.names = TRUE))] %>%
  lapply(.,read.csv) %>% bind_rows() %>%
  write.csv(., "./input/cleaned/clean_lencomps.csv",row.names = FALSE)

list.files("./input/temp/", full.names = TRUE)[grep("_agecomps", list.files("./input/temp/", full.names = TRUE))] %>%
  lapply(.,read.csv) %>% bind_rows() %>%
  write.csv(., "./input/cleaned/clean_agecomps.csv",row.names = FALSE)