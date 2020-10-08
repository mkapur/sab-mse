require(ggplot2)
require(patchwork)
library(gridExtra)
require(dplyr)

years <- 1960:2019
nyear <- length(years)
tEnd <- length(years)
age <- 0:70 
nage <- length(age)


reps$Ninit_ais[,,1] %>%
  data.frame() %>%
  mutate('Age' = age) %>%
  reshape2::melt(id = c('Age')) %>%
  ggplot(., aes(x = Age, y = value, color = variable )) +
  scale_color_manual(values = subareaPal) +
  geom_line(lwd = 2) + 
  labs(x = 'Age in Initial Years',y = 'Initial Numbers', color = 'subarea') +
  ggsidekick::theme_sleek()+
  facet_wrap(~variable,scales = 'free_y')

reps$SSB_yk %>%
  data.frame() %>%
  mutate('Yr' = years) %>%
  reshape2::melt(id = c('Yr')) %>%
  ggplot(., aes(x = Yr, y = value, color = variable )) +
  scale_color_manual(values = demPal) +
  geom_line(lwd = 2) + labs(x = 'Modeled Year',y = 'SSB', color = 'stock') +
  ggsidekick::theme_sleek() + theme( legend.position = c(0.8,0.8))



## ssb_ym with compare ----
spmat <- data.frame(subarea = c('A1',"A2","B2","B1","C2","C1"),
                    stock = c("R4","R3","R3","R2","R2","R1"),
                    mgmt = c("AK","AK", rep("BC",2), rep("CC",2)))

assSB <- read.csv(here('input','downloads','AssessmentDat_thru2018.csv'),fileEncoding="UTF-8-BOM") %>%
  mutate(REG = substr(Index,1,2), assSSBMT = Value) %>% filter(Type == 'SpawnBiomass')

names(reps$SSB_yi) <- spmat$subarea
SSB_ym0 <-  reps$SSB_yi %>%
  data.frame() %>%
  mutate('Yr' =years) %>%
  reshape2::melt(id = c('Yr')) %>%
  merge(., spmat, by.x = "variable", by.y = "subarea") %>% 
  group_by(Yr, mgmt) %>%
  summarise(totSSBkg = sum(value), 
            totSSBmt = totSSBkg/1e6) %>% select(Yr, mgmt, totSSBmt)
levels(SSB_ym0$mgmt) <- c('AI','AK','BC','CC','WC')

SSB_ym0$mgmt[SSB_ym0$mgmt == 'AI'] <- 'AK'  
SSB_ym0$mgmt[SSB_ym0$mgmt == 'CC'] <- 'WC'  
SSB_ym <- SSB_ym0 %>% group_by(Yr, mgmt) %>%
  summarise(omSSBMT = sum(totSSBmt))

merge(assSB, SSB_ym, by.x = c('Year','REG'), by.y = c('Yr','mgmt')) %>%
  select(Year, REG, CV,assSSBMT, omSSBMT) %>%
  filter(Year > 1980) %>%
  ggplot(., aes(x = Year, y = assSSBMT, color = REG)) +
  geom_line(aes(y = omSSBMT),lwd = 1.1) +
  geom_errorbar(aes(ymin = assSSBMT-CV*assSSBMT,ymax= assSSBMT+CV*assSSBMT, color = REG)) +
  scale_color_manual(values = mgmtPal)+
  geom_point()+  ggsidekick::theme_sleek() +
  labs(x = 'Modeled Year',y = 'SSB', color = 'Mgmt Region') +
  facet_wrap(~REG,scales = 'free_y')


## catch pred by fleet ----
catch_yf_predt <- data.frame(reps$catch_yf_pred)
names(catch_yf_predt) <- df$fltnames_fish
catch_yf_predt <- catch_yf_predt %>%
  mutate(Year = years) %>%
  group_by(Year) %>%
  mutate("AK_FIX (aggregate)" = sum(AK_FIX_W, AK_FIX_E),
         "AK_TWL (aggregate)"= sum(AK_TWL_W, AK_TWL_E)) %>%
  select(-AK_TWL_W,-AK_TWL_E,-AK_FIX_W,-AK_FIX_E) %>%
  melt(id = 'Year') %>%
  mutate(Type = 'PRED') %>%
  mutate(REG = substr(variable,0,2)) #%>%
# filter(REG == 'BC') #%>% View()

catch_yf_obst <- df$catch_yf_obs %>%  data.frame() %>%select(-Year) 

catch_yf_obst <- catch_yf_obst %>%
  mutate(Year = years) %>%
  group_by(Year) %>%
  mutate("AK_FIX (aggregate)" = sum(AK_FIX_W, AK_FIX_E),
         "AK_TWL (aggregate)"= sum(AK_TWL_W, AK_TWL_E)) %>%
  select(-AK_TWL_W,-AK_TWL_E,-AK_FIX_W,-AK_FIX_E) %>%
  melt(id = 'Year') %>%
  ## convert CV to SD via CV = mean/sd
  mutate(Type = 'OBS', 
         lci = value - 1.96*(0.1*value),
         uci = value + 1.96*(0.1*value)) %>%
  mutate(REG = substr(variable,0,2)) #%>% head()

ggplot(data = catch_yf_obst, 
       aes(x = Year, y = value, color = variable)) +
  geom_line(data = catch_yf_predt, lwd = 0.75) +
  scale_color_manual(values = fishfltPal) +
  scale_x_continuous(limits = c(1980,2020)) +
  geom_point(pch = 1, fill = NA, col = 'black') +
  geom_errorbar(aes(ymin = lci, ymax = uci), col = 'black') +
  theme_sleek() + 
  theme(legend.position = 'none')+
  labs(y = 'catch', color = 'Fishing Fleet')+
  facet_wrap(~variable, scales = "free_y")
