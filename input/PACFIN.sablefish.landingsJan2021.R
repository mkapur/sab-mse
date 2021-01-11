##############################################################################################################
# Purpose: Summarize Sablefish Commercial Landings 81- present from PacFIN
# Written: Melissa Haltuch, 2019
# Update M Kapur 2021
##############################################################################################################
#Note that OR catches from 1981-1986 should not be used, use the historical catch reconstrution instead

#Read in and setup PacFIN Catches 81 - present
#Catches can be double counted here so be careful
#e.g. "ALL" in ARID is all catch by gear and state
#e.g. catches in 1A = CP as this is essentially the same area
#e.g. In the PCID column "AWA" stands for All WA ports.
require(here)
require(dplyr)
rm(list=ls())

#read in pacfin catches

newDat = load(here('input','raw_data','catch',"PacFIN.SABL.Catch.03.Dec.2020.RData")) # takes a sec
newDat <- PacFIN.SABL.Catch.03.Dec.2020$PacFIN.INPFC.Summary.Catch
newDat %>% select("COUNCIL","YEAR","PERIOD","SPID","ARID","GRID","GRGROUP","PCID","CATCH.KG")
names(newDat)
DataYr <- 2020 #2018

load( here('input','raw_data','catch',"PacFIN.SABL.Catch.INPFC.28.Mar.2019.dmp"))
oldDat = PacFIN.SABL.Catch.INPFC.28.Mar.2019
# names(newDat)
#"COUNCIL"  "YEAR"     "PERIOD"   "SPID"     "ARID"     "GRID"     "GRGROUP"  "PCID"     "CATCH.KG"
#PERIOD IS MONTH HERE!
newDat$Catch.MT <- newDat$CATCH.KG/1000	#convert to MT
Nyears <- DataYr-1981+1

#subset data to keep only areas of interest
newDat$flag <- 0 				
newDat$flag[(newDat$ARID %in% c("CP","MT","EK","CL","VN","UI"))] <- 1 	
newDat <- newDat[newDat$flag %in% 1,]      # remove unused records
#"CP" - conception
#"MT" - monetery
#"EK" - eureka
#"CL" - columbia
#"VN" - vancouver
#"UI" - unknown pcouncil

#aggregate by gear groups
# unique(newDat$GRGROUP)
#[1] "POT" "TWL" "NET" "HKL" "TWS" "TLS" "MSC" "DRG"
#HKL gear fleet = TLS, and HKL
#POT gear fleet
#TWL fleet = add non-trawl nets, (NET)to trawl (TWL), shrimp trawl (TWS),  "MSC",  "DRG"
#Adding a hake sablefish bycatch fleet using the DAHL SECTOR = 3 code to select the appropraite vessels

#catches by port of landing
Catch <- newDat[,c("YEAR","PCID","Catch.MT", "PERIOD", "GRGROUP","DAHL_SECTOR")]

#assign catches to state using PCID
Catch$State <- "NA"
Catch$State[Catch$PCID %in% c("AWA")] <- "WA"
Catch$State[Catch$PCID %in% c("AOR")] <- "OR"
Catch$State[Catch$PCID %in% c("ACA")] <- "CA"

#sum catches for each STATE
Catch.sum.state <- aggregate(Catch$Catch.MT,list(year=Catch$YEAR, state=Catch$State, gear=Catch$GRGROUP) ,sum)
names(Catch.sum.state) <- c("Year", "State", "Gear", "Catch.MT")
write.csv(Catch.sum.state, file = here('input','raw_data','catch',"Catch.sum.state.allgears.csv"))

#subset state hkl gear table
Catch.state.hkl <- subset(Catch, Catch$GRGROUP %in% c("HKL","TLS"))
Catch.sum.state.hkl <- aggregate(Catch.state.hkl$Catch.MT,list(year=Catch.state.hkl$YEAR, state=Catch.state.hkl$State) ,sum)
names(Catch.sum.state.hkl) <- c("Year", "State", "Catch.MT")
write.csv(Catch.sum.state.hkl, file = here('input','raw_data','catch',"Catch.sum.state.hklgears.csv"))

#subset state pot gear table
Catch.state.pot <- subset(Catch, Catch$GRGROUP %in% c("POT"))
Catch.sum.state.pot <- aggregate(Catch.state.pot$Catch.MT,list(year=Catch.state.pot$YEAR, state=Catch.state.pot$State) ,sum)
names(Catch.sum.state.pot) <- c("Year", "State", "Catch.MT")
write.csv(Catch.sum.state.pot, file = here('input','raw_data','catch',"Catch.sum.state.potgears.csv"))



# separate DAHL_SECTOR = 3 AND DAHL_SECTOR != 3
Catch.w3 <- Catch[ which(Catch$DAHL_SECTOR == "03"),]
Catch.no3 <- Catch[ which(Catch$DAHL_SECTOR != "03"),]

#subset state trawl gear table
Catch.state.trawl <- subset(Catch.no3, Catch.no3$GRGROUP  %in% c("TWL", "TWS", "NET", "MSC", "DRG"))
Catch.sum.state.trawl <- aggregate(Catch.state.trawl$Catch.MT,list(year=Catch.state.trawl$YEAR, state=Catch.state.trawl$State) ,sum)
names(Catch.sum.state.trawl) <- c("Year", "State", "Catch.MT")
write.csv(Catch.sum.state.trawl, file = here('input','raw_data','catch',"Catch.sum.state.trawlgears.csv"))

Catch.state.trawl.only <- subset(Catch.no3, Catch.no3$GRGROUP  %in% c("TWL"))
Catch.sum.state.trawl.only <- aggregate(Catch.state.trawl.only$Catch.MT,list(year=Catch.state.trawl.only$YEAR, state=Catch.state.trawl.only$State) ,sum)
names(Catch.sum.state.trawl.only) <- c("Year", "State", "Catch.MT")
write.csv(Catch.sum.state.trawl.only, file = here('input','raw_data','catch',"Catch.sum.state.trawlgears.only.csv"))

#subset state trawl gears with dalh_sector = 3, these are hake trips that catch small sablefish
Catch.state.trawl.hake <- subset(Catch.w3, Catch.w3$GRGROUP  %in% c("TWL", "TWS", "NET", "MSC", "DRG"))
Catch.sum.state.trawl.hake <- aggregate(Catch.state.trawl.hake$Catch.MT,list(year=Catch.state.trawl.hake$YEAR, state=Catch.state.trawl.hake$State) ,sum)
names(Catch.sum.state.trawl.hake) <- c("Year", "State", "Catch.MT")
write.csv(Catch.sum.state.trawl.hake, file = here('input','raw_data','catch',"Catch.sum.state.trawlgears.hake.csv"))



## reweight via discards
disrate_cs <- read.csv(here('input','raw_data','catch',
                            "sablefish_OB_DisRatios_cs_2019_Coastwide_trawl_fixed_2019-06-27.csv")) %>% 
  filter(gear3 != 'Trawl')
disrate_ncs <- read.csv(here('input','raw_data','catch',
                             "sablefish_OB_DisRatios_ncs_2019_Coastwide_trawl_fixed_2019-06-27.csv")) %>%
  filter(gear3 != 'Trawl')


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


## sum HKL + pot as before, now reweighted
rbind(Catch.state.pot,Catch.state.hkl) %>%
  group_by(YEAR) %>%
  summarise(Obs=sum(Catch.MT)) %>% tail()
  # filter(YEAR > 2005) %>%
  # merge(., filter(wc$catch, 
  #                 # Yr >2005, 
  #                 Fleet_Name == 'FIX') %>% 
  #         select(Yr, Obs), by.x = 'YEAR', by.y = 'Yr')  %>% 
  write.csv(., 
            file = here('input','raw_data','catch',"Catch.sum.state.HKL_POT.csv"))


