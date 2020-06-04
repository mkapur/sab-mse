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


## Comps -- would be great to get these @ subarea if possible ----
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
