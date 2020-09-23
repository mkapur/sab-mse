require(ggplot2)
require(patchwork)
library(gridExtra)
require(dplyr)

## N init and Nzero ----
pNinit <- Ninit_ais[,,1] %>%
  data.frame() %>%
  mutate('Age' = age) %>%
  reshape2::melt(id = c('Age')) %>%
  ggplot(., aes(x = Age, y = value, color = variable )) +
  scale_color_manual(values = subareaPal) +
  geom_line(lwd = 2) + 
  labs(x = 'Age in Initial Years',y = 'Initial Numbers', color = 'subarea') +
  ggsidekick::theme_sleek()

ggsave(pNinit,
       file = here('figs',
                   "NZeroNinit_Move_BothSex_normalizedOmega.png"),
       width = 6, height = 4, unit = 'in',
       dpi = 420)

pNzeroF <- N_0ais[,,1] %>% data.frame() %>%
  mutate('Age' = 1:71) %>%
  reshape2::melt(id = c('Age')) %>%
  ggplot(., aes(x = Age, y = value, color = variable )) +
  geom_line(lwd = 1) + 
  scale_color_manual(values = subareaPal) +
  labs(x = 'Age in Year 0',y = 'Unfished Numbers', color = 'subarea') +
  ggsidekick::theme_sleek()+theme(legend.position = 'none')

pNzeroM <- sim.data$N_0ais[,,2] %>% data.frame() %>%
  mutate('Age' = 1:71) %>%
  reshape2::melt(id = c('Age')) %>%
  ggplot(., aes(x = Age, y = value, color = variable )) +
  geom_line(lwd = 1) + 
  scale_color_manual(values = subareaPal) +
  labs(x = 'Age in Year 0',y = 'Unfished Numbers', color = 'subarea') +
  ggsidekick::theme_sleek() 

ggsave(pNzeroF  | pNzeroM,
       file = here('figs',
                   "NZero_Move_2Sex_normalizedOmega.png"),
       width = 6, height = 4, unit = 'in',
       dpi = 420)

## N at age by area by year (females) ----
png(file = here('figs','N_age_years1-15.png'),
    width = 10, height = 8, unit = 'in', res = 420)
par(mfrow = c(2,3))
for(i in 1:6){
  plot(N_yais_beg[1,,i,1],
       type = 'l',
       lwd = 2, 

       col = 'black',
       main = inames[i], 
       col.main =subareaPal[i], 
       ylim = c(0,20000),
       xlim = c(0,70),
       xlab = "Age", 
       ylab = 'Numbers')
  box(which = 'plot', lty = 'solid', 
      col = subareaPal[i],
      lwd = 2)
  for (y in 2:15) {
    lines(N_yais_beg[y, , 1, 1],
          # lty = c(2:6)[y],
          col = gray.colors(25, start = 0.1, end = 0.9)[y],
          lwd = 2)
  }

}
dev.off()

## total nums in area by year
png(file = here('figs','N_iy.png'),
    width = 10, height = 8, unit = 'in', res = 420)
par(mfrow = c(2,3))
for(i in 1:6){
  plot(rowSums(N_yais_beg[,,i,]),
       type = 'l',
       lwd = 2, 
       col =subareaPal[i],
       main = inames[i], 
       col.main =subareaPal[i], 
       ylim = c(0,150000),
       xlim = c(0,nyear),
       xlab = "Model Year", 
       ylab = 'Numbers (M+F)')
  # lines(rowSums(N_yais_beg[,,i,2]),
  #      type = 'l',
  #      lwd = 3, 
  #      col =subareaPal[i])
}
dev.off()

## SSB y ----
pSSByi <- SSB_yi %>%
  data.frame() %>%
  mutate('Yr' = 1:nrow(SSB_yi)) %>%
  reshape2::melt(id = c('Yr')) %>%
  ggplot(., aes(x = Yr, y = value, color = variable )) +
  scale_color_manual(values = subareaPal) +
  geom_line(lwd = 2) + 
  labs(x = 'Modeled Year',y = 'SSB', color = 'subarea') +
  ggsidekick::theme_sleek() + theme( legend.position = c(0.8,0.8))

pSSByk <- SSB_yk[1:nyear,] %>%
  data.frame() %>%
  mutate('Yr' = 1:nyear) %>%
  reshape2::melt(id = c('Yr')) %>%
  ggplot(., aes(x = Yr, y = value, color = variable )) +
  scale_color_manual(values = demPal) +
  geom_line(lwd = 2) + labs(x = 'Modeled Year',y = 'SSB', color = 'stock') +
  ggsidekick::theme_sleek() + theme( legend.position = c(0.8,0.8))

## R y ----
pRyi <- R_yi %>% data.frame() %>%
  mutate('Yr' = 1:nrow(.)) %>%
  reshape2::melt(id = c('Yr')) %>%
  ggplot(., aes(x = Yr, y = value, color = variable )) +
  scale_color_manual(values = subareaPal) +
  geom_line(lwd = 2) + 
  labs(x = 'Model Year',y = 'Recruits', color = 'subarea') +
  ggsidekick::theme_sleek()+ theme( legend.position = c(0.8,0.8))

pRyk <- R_yk %>% data.frame() %>%
  mutate('Yr' = 1:nrow(.)) %>%
  reshape2::melt(id = c('Yr')) %>%
  ggplot(., aes(x = Yr, y = value, color = variable )) +
  scale_color_manual(values = demPal) +
  geom_line(lwd = 2) + 
  labs(x = 'Model Year',y = 'Recruits', color = 'stock') +
  ggsidekick::theme_sleek()+ theme( legend.position = c(0.8,0.8))


## SRR ----
getRYK <- function(h,r0,ssb,ssb0){
  RYK = (4*h*r0*ssb)/
    (ssb0*(1-h)+
       ssb*(5*h-1))
  return(RYK)
}
# plot(0:1000, getRYK(h = 0.5, r0 =1500, ssb = 0:1000, ssb0 = 500),
#      ylim = c(0,3000), xlim = c(0,1000))

SRR <- data.frame('SSByk' = 0:10000, 
                  'RYK' = getR)

pSRR <- R_yk[1:15,] %>% 
  data.frame() %>%
  mutate('Yr' = 1:nrow(.))  %>%
  reshape2::melt(.,id = c('Yr')) %>%
  mutate(RYK = value) %>%
  select(-value) %>%
  bind_cols(.,
            SSB_yk[1:15,] %>% 
              data.frame() %>%
              mutate('Yr' = 1:nrow(.))  %>%
              reshape2::melt(.,id = c('Yr')) %>%
              mutate(SSByk = value) %>%
              select(-value,-variable)) %>%
  ggplot(., aes(x = SSByk, y = RYK, color = variable, group = Yr...1)) +
  scale_color_manual(values = demPal) +
  geom_point() +
  labs(x = 'SSB',y = 'Recruits #', color = 'stock') +
  scale_y_continuous(limits = c(0,25000)) +
  # scale_x_continuous(limits = c(0,400000)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~ variable)

ggsave(pSRR,
       file = here('figs',
                   "SRR.png"),
       width = 6, height = 6, unit = 'in',
       dpi = 420)

## deterministic/expected length at ages
pLAA1F <- Length_yais_beg[,,1,1] %>% 
  data.frame() %>%
  mutate('Yr' = 1:nrow(.)) %>%
  reshape2::melt(id = c('Yr')) %>%
  ggplot(., aes(x = as.numeric(substr(variable,2,3))-1, y = value)) +
  geom_point(color = 'orange') +
  labs(x = 'Age',y = 'Length') +
  scale_y_continuous(limits = c(0,75))+
  ggsidekick::theme_sleek()

pLAA2F <- Length_yais_beg[,,2,1] %>% data.frame() %>%
  mutate('Yr' = 1:nrow(.)) %>%
  reshape2::melt(id = c('Yr')) %>%
  ggplot(., aes(x = as.numeric(substr(variable,2,3))-1, y = value, color = Yr )) +
  geom_point(color = 'orange') +
  labs(x = 'Age',y = 'Length', color = 'orange') +
  scale_y_continuous(limits = c(0,75))+
  ggsidekick::theme_sleek()

pLAA1M <- Length_yais_beg[,,1,2] %>% data.frame() %>%
  mutate('Yr' = 1:nrow(.)) %>%
  reshape2::melt(id = c('Yr')) %>%
  ggplot(., aes(x = as.numeric(substr(variable,2,3))-1, y = value)) +
  geom_point(color = 'blue') +
  labs(x = 'Age',y = 'Length', color = 'Year') +
  scale_y_continuous(limits = c(0,75))+
  ggsidekick::theme_sleek()

pLAA2M <- Length_yais_beg[,,2,2] %>% data.frame() %>%
  mutate('Yr' = 1:nrow(.)) %>%
  reshape2::melt(id = c('Yr')) %>%
  ggplot(., aes(x = as.numeric(substr(variable,2,3))-1, y = value, color = Yr )) +
  geom_point(color = 'blue') +
  labs(x = 'Age',y = 'Length', color = 'Year') +
  scale_y_continuous(limits = c(0,75))+
  ggsidekick::theme_sleek()

(pLAA1F  | pLAA1M) /
  (pLAA2F  | pLAA2M)

## LAA dist ----
pA1 <- list() ## for area 1
for(y in 1:5){ #45:dim(LengthAge_alyi_beg)[3]){ ## loop years
  a1 <- 
    LengthAge_alyis_beg[,,y,1,1] %>%
    melt() %>%
    group_by(Var2) %>% 
    mutate(sumP = sum(value), ## total prob within A-L bin
           pbin = value/sumP) ## relative prob within A-L bin
  
  ggplot(a1,aes( x = Var2, y = pbin, color = factor(Var1))) +
    # geom_histogram(stat ='identity', position = 'stack') +
    geom_line() + 
    geom_area(aes(fill=factor(Var1))) +
    # geom_density( )+
    labs(x = 'len', y = 'prob(Bin)', fill = 'age', 
         title = paste("year ",y), 
         subtitle= 'subarea 1') +
    theme_sleek()
  
  ggplot(a1,aes( x = Var1, y = pbin, color = factor(Var2))) +
    # geom_histogram(stat ='identity', position = 'stack') +
    geom_line() + 
    # geom_area(aes(fill=factor(Var2))) +
    # geom_density( )+
    labs(x = 'age', y = 'prob(Bin)', fill = 'len', 
         title = paste("year ",y), 
         subtitle= 'subarea 1') +
    theme_sleek()
  
  pA1[[y]] <- ggplot(a1,aes(x = Var1, y = Var2, fill = pbin)) +
    geom_tile() +
    labs(x = 'age', y = 'len', 
         title = paste("year ",y), 
         subtitle= 'subarea 1') +
    theme_sleek()
  rm(a1)
  # p[[i]] <- qplot(1:10,10:1,main=i)
}
# png(here("figs","LAA_Dist_A1.png"), 
#     height = 10, width = 10, unit = 'in', res = 420)
do.call(gridExtra::grid.arrange,pA1)
dev.off()


pA2 <- list() ## for area 1
for(y in 1:10){ #45:dim(LengthAge_alyi_beg)[3]){ ## loop years
  a1 <- 
    LengthAge_alyi_beg[,,y,2] %>%
    melt() %>%
    group_by(Var2) %>% 
    mutate(sumP = sum(value), pbin = value/sumP)
  pA2[[y]] <- ggplot(a1,aes(x = Var1, y = Var2, fill = pbin)) +
    geom_tile() +
    labs(x = 'age', y = 'len', 
         title = paste("year ",y), 
         subtitle= 'subarea 2') +
    theme_sleek()
  rm(a1)
  # p[[i]] <- qplot(1:10,10:1,main=i)
}
png(here("figs","LAA_Dist_A2.png"), 
    height = 10, width = 10, unit = 'in', res = 420)
do.call(grid.arrange,pA2)
dev.off()

## CATCH pred by fleet ----


catch_yf_predt <- data.frame(catch_yf_pred)

catch_yf_predt <- catch_yf_predt %>%
  mutate(Year = year) %>%
  melt(id = 'Year') %>%
  mutate(Type = 'PRED') %>%
  mutate(REG = substr(variable,0,2)) #%>%
  filter(REG == 'BC') #%>% View()

catch_yf_obst <- catch_yf_obs[1:length(year),] %>%  data.frame() %>%select(-Year) 

catch_yf_obst <- catch_yf_obst %>%
  mutate(Year = year) %>%
  melt(id = 'Year') %>%
  ## convert CV to SD via CV = mean/sd
  mutate(Type = 'OBS', 
         lci = value - 1.96*(0.1*value),
         uci = value + 1.96*(0.1*value)) %>%
  mutate(REG = substr(variable,0,2)) #%>%
  filter(REG == 'BC') #%>% View()

ggplot(data = catch_yf_obst, 
       aes(x = Year, y = value, color = variable)) +
  geom_line(data = catch_yf_predt, lwd = 0.75) +
  scale_color_manual(values = fishfltPal) +
  geom_point(pch = 1, fill = NA, col = 'black') +
  geom_errorbar(aes(ymin = lci, ymax = uci), col = 'black') +
  theme_sleek() + 
  theme(legend.position = 'none')+
  labs(y = 'catch', color = 'Fishing Fleet')+
  facet_wrap(~variable, scales = "free_y")

ggsave(last_plot(),
       file = here('figs',
                   paste0('catch_fits_noLenSel_',
                          'v1=',v1,'Fmax=',Fmax,Sys.Date(),'.png')),
       width = 10, height = 6, unit = 'in',
       dpi = 420)
## catch pred by m ----
catch_yf_predm <- catch_yf_predt %>% 
  group_by(Year, REG) %>%
  summarise(totC = sum(value)) %>%  mutate(Type = 'PRED') 
catch_yf_obsm <- catch_yf_obst %>% 
  group_by(Year, REG) %>%
  summarise(totC = sum(value)) %>%
  mutate(Type = 'OBS', 
         lci = totC - 1.96*(0.1*totC),
         uci = totC + 1.96*(0.1*totC)) 

ggplot(data = catch_yf_obsm, 
       aes(x = Year, y = totC, color = REG)) +
  geom_line(data = catch_yf_predm, lwd = 1.1) +
  scale_color_manual(values = mgmtPal) +
  geom_point(pch = 1, fill = NA, col = 'black') +
  geom_errorbar(aes(ymin = lci, ymax = uci), col = 'black') +
  theme_sleek() + 
  theme(legend.position = 'none')+
  labs(y = 'catch', color = 'Fishing Fleet')+
  facet_wrap(~REG, scales = "free_y")

ggsave(last_plot(),
       file = here('figs',
                   paste0('catchm_fits_noLenSel_',
                          'v1=',v1,'Fmax=',Fmax,Sys.Date(),'.png')),
       width = 10, height = 6, unit = 'in',
       dpi = 420)
## survey preds ----


survey_yf_predt <- data.frame(survey_yf_pred)

survey_yf_predt <- survey_yf_predt %>%
  mutate(Year = year) %>%
  melt(id = 'Year') %>%
  mutate(Type = 'PRED') %>%
  mutate(REG = substr(variable,0,2)) #%>%
filter(REG == 'BC') #%>% View()

survey_yf_obst <- surv_yf_obs %>%  data.frame() 

survey_yf_obst <- survey_yf_obst %>%
  mutate(Year = year) %>%
  melt(id = 'Year') %>%
  ## convert CV to SD via CV = mean/sd
  mutate(Type = 'OBS', 
         lci = value - 1.96*(0.1*value),
         uci = value + 1.96*(0.1*value)) %>%
  mutate(REG = substr(variable,0,2)) #%>%
filter(REG == 'BC') #%>% View()

ggplot(data = survey_yf_obst, 
       aes(x = Year, y = value, color = variable)) +
  geom_line(data = survey_yf_predt, lwd = 0.75) +
  scale_color_manual(values = fishfltPal) +
  geom_point(pch = 1, fill = NA, col = 'black') +
  geom_errorbar(aes(ymin = lci, ymax = uci), col = 'black') +
  theme_sleek() + 
  theme(legend.position = 'none')+
  labs(y = 'survey', color = 'Fishing Fleet')+
  facet_wrap(~variable, scales = "free_y")

ggsave(last_plot(),
       file = here('figs',
                   paste0('survey_fits_selMult_',
                          'v1=',v1,'Fmax=',Fmax,Sys.Date(),'.png')),
       width = 10, height = 6, unit = 'in',
       dpi = 420)

# survey_yf_predt$value[survey_yf_predt$variable %in% c('BC_LL','BC_TRAP','BC_TWL')] <-
#   survey_yf_predt$value[survey_yf_predt$variable %in% c('BC_LL','BC_TRAP','BC_TWL')]/2E2
# 
# survey_yf_predt$value[survey_yf_predt$variable %in% c('AK_FIX_W','AK_FIX_E','AK_TWL_W')] <-
#   survey_yf_predt$value[survey_yf_predt$variable %in% c('AK_FIX_W','AK_FIX_E','AK_TWL_W')]*2.2


# 
# simdata$SSB_yk %>% 
#   data.frame() %>%
#   mutate('Yr' = 1:53) %>%
#   reshape2::melt(id = c('Yr')) %>%  
#   ggplot(., aes(x = Yr, y = value, color = variable )) +
#   geom_line(lwd = 2) + labs(x = 'Modeled Year',y = 'SSB', color = 'stock') +
#   ggsidekick::theme_sleek()

# }