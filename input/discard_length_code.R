## Get lengths from discard. Going to use the "weighted" columns
require(tidyr); require(data.table)
ldis <- read.csv("./sablefish_discard_lengths.csv") ## first tab
ndis  <- read.csv("./sablefish_discard_nuniquetrips.csv") ## use N_unique_trips for nsamp
names(ldis)[1] <- 'Year'
names(ndis)[1] <- 'Gear'

ldis.temp <- ldis %>% filter(Gear == 'HookAndPot') %>% 
  select(Year, Prop.wghtd,Lenbin)  %>% 
  ## similar to making a pivot table
  pivot_wider(id_cols = c(Year), names_from = c(Lenbin), values_from= Prop.wghtd) %>%
  mutate(X18 = 0) %>% ## add empty 18 bin
  select(Year, X18, everything())


## duplicate the columms for second sex
names(ldis.temp)[2:ncol(ldis.temp)] <- paste0('X',seq(18,90,2)) ## rename cols
ldis.temp2 <- cbind(ldis.temp, ldis.temp[,-1]) ## bind them again
ldis.temp2[,-1] <- round(ldis.temp2[,-1],6) ## round values
names(ldis.temp2[39:ncol(ldis.temp2)]) <- paste0(names(ldis.temp2)[2:38],".1")

Nsamp <- ndis %>% filter(Gear == 'HookAndPot') %>% select(N_unique_Trips)
cbind(ldis.temp2,Nsamp) %>% data.frame() %>%
  ## sex combined -- just duplicated as in examples
  mutate(month = 7, fleet = 1, sex = 0, part = 1, Nsamp = N_unique_Trips ) %>%
  select(-N_unique_Trips) %>%
  select(Year, month, fleet, sex, part, Nsamp,  X18, everything()) 