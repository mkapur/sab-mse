## Dataprep
## M S KAPUR
## wrangling of raw synthesis/catch/misc stuff received from management regions
## this does NOT get it in TMB ready format, but writes the CSVs that serve as the baseline
## and therefore only need to get updated if we've received new data from a management region

# Landings & Discards ----

## ALASKA ##
akcat <- read.csv("./input/raw/AK_catch_ByGear_and_ByArea.csv")
akdis <- read.csv("./input/raw/AK_discards.csv")

ak_catdis <- data.frame(Year = akdis$Year, Gear = akdis$Gear, akdis[grep(paste("BSAI*","GOA*",sep= "|"), names(akdis))]) %>%
  select(-BSAI_PctDiscard, -GOA_PctDiscard) %>%
  reshape2::melt(id = c("Year","Gear")) %>%
  mutate(variable = as.character(variable)) %>%
  mutate(type = stringr::str_split_fixed(variable, "_", 2)[,2],
        stock_temp =  stringr::str_split_fixed(variable, "_", 2)[,1]) %>%

  mutate(stock = ifelse(stock_temp == 'BSAI', "R5","R4"),
         subarea = ifelse(stock_temp == 'BSAI', "A1","A2"),
         mgmt = "AK") %>%
  select(Year, type, subarea, stock, mgmt,  everything(), -stock_temp, -variable) 

ak_catdis %>%
  filter(type == 'Catch') %>%
  write.csv(., file = paste0("input/cleaned/ak_catch_",Sys.Date(),".csv"), row.names = FALSE)

ak_catdis %>%
  filter(type == 'Discard') %>%
  write.csv(., file = paste0("input/cleaned/ak_discard_",Sys.Date(),".csv"), row.names = FALSE)
