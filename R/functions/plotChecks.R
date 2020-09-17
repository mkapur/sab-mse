require(ggplot2)
require(patchwork)
library(gridExtra)


pNinit <- Ninit_ais[,,1] %>%
  data.frame() %>%
  mutate('Age' = age) %>%
  reshape2::melt(id = c('Age')) %>%
  ggplot(., aes(x = Age, y = value, color = variable )) +
  scale_color_manual(values = subareaPal) +
  geom_line(lwd = 2) + 
  labs(x = 'Age in Initial Years',y = 'Initial Numbers', color = 'subarea') +
  ggsidekick::theme_sleek()


pNzeroF <- N_0ais[,,1] %>% data.frame() %>%
  mutate('Age' = 1:71) %>%
  reshape2::melt(id = c('Age')) %>%
  ggplot(., aes(x = Age, y = value, color = variable )) +
  geom_line(lwd = 1) + 
  scale_color_manual(values = subareaPal) +
  labs(x = 'Age in Year 0',y = 'Unfished Numbers', color = 'subarea') +
  ggsidekick::theme_sleek()+theme(legend.position = 'none')
pNzeroM <- N_0ais[,,2] %>% data.frame() %>%
  mutate('Age' = 1:71) %>%
  reshape2::melt(id = c('Age')) %>%
  ggplot(., aes(x = Age, y = value, color = variable )) +
  geom_line(lwd = 1) + 
  scale_color_manual(values = subareaPal) +
  labs(x = 'Age in Year 0',y = 'Unfished Numbers', color = 'subarea') +
  ggsidekick::theme_sleek() 

ggsave(pNzeroF  | pNzeroM,
       file = here('figs',
                   "NZero_Move_2Sex.png"),
       width = 6, height = 4, unit = 'in',
       dpi = 420)

## N at age by area by year (females)
png(file = here('figs','N_age_years1-5.png'),
    width = 10, height = 8, unit = 'in', res = 420)
par(mfrow = c(2,3))
for(i in 1:6){
  plot(N_yais_beg[1,,i,1],
       type = 'l',
       lwd = 2, 
       col = subareaPal[i],
       main = inames[i], 
       col.main =subareaPal[i], 
       ylim = c(0,500),
       xlim = c(0,70),
       xlab = "Age", 
       ylab = 'Numbers')
  for(y in 2:2){
    lines(N_yais_beg[y,,1,1],  col = subareaPal[i], lwd =2)
  }
}
dev.off()

pNage2 <- N_yai_beg[,,2] %>%
  data.frame() %>%
  mutate('Yr' = 1:nrow(.)) %>%
  reshape2::melt(id = c('Yr')) %>%
  mutate(age = as.numeric(substr(variable,2,2))-1) %>%
  ggplot(., aes(x = age, y = value, group = factor(Yr), color = factor(Yr))) +
  geom_point() +
  geom_line()  +
  # geom_boxplot() +
  scale_color_grey()+
  labs(x = 'Age in Year',y = 'Numbers', title = "AREA2") +
  ggsidekick::theme_sleek()+ theme(legend.position = 'none')

pSSByi <- SSB_yi[1:52,] %>%
  data.frame() %>%
  mutate('Yr' = 1:52) %>%
  reshape2::melt(id = c('Yr')) %>%
  ggplot(., aes(x = Yr, y = value, color = variable )) +
  geom_line(lwd = 2) + labs(x = 'Modeled Year',y = 'SSB', color = 'subarea') +
  ggsidekick::theme_sleek()

pSSByk <- SSB_yk[1:nyear,] %>%
  data.frame() %>%
  mutate('Yr' = 1:nyear) %>%
  reshape2::melt(id = c('Yr')) %>%
  ggplot(., aes(x = Yr, y = value, color = variable )) +
  geom_line(lwd = 2) + labs(x = 'Modeled Year',y = 'SSB', color = 'stock') +
  ggsidekick::theme_sleek()

pRyi <- R_yi %>% data.frame() %>%
  mutate('Yr' = 1:nrow(.)) %>%
  reshape2::melt(id = c('Yr')) %>%
  ggplot(., aes(x = Yr, y = value, color = variable )) +
  geom_line(lwd = 2) + 
  labs(x = 'Model Year',y = 'Recruits', color = 'subarea') +
  ggsidekick::theme_sleek()

pRyk <- R_yk %>% data.frame() %>%
  mutate('Yr' = 1:nrow(.)) %>%
  reshape2::melt(id = c('Yr')) %>%
  ggplot(., aes(x = Yr, y = value, color = variable )) +
  geom_line(lwd = 2) + 
  labs(x = 'Model Year',y = 'Recruits', color = 'stock') +
  ggsidekick::theme_sleek()


pSRR <- cbind(R_yk, SSB_yk) %>% 
  data.frame() %>%
  mutate('Yr' = 1:nrow(.), 
         RECS1 = X1, RECS2 = X2, SSBS1 = X3, SSBS2 = X4) %>%
  select(-X1,-X2,-X3,-X4) %>%
  reshape2::melt(id = c('Yr')) %>%
  mutate(variable2 = substr(variable,1,3), stock = substr(variable,4,5)) %>%
  select(-variable) %>%
  tidyr::pivot_wider(names_from = variable2) %>%
  ggplot(., aes(x = SSB, y = REC, color = Yr )) +
  geom_point() +
  labs(x = 'SSB',y = 'Recruits #', color = 'Y') +
  ggsidekick::theme_sleek() +
  facet_wrap(~ stock, scales = 'free')


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

## distributions
pA1 <- list() ## for area 1
for(y in 1:10){ #45:dim(LengthAge_alyi_beg)[3]){ ## loop years
  a1 <- 
    LengthAge_alyis_beg[,,y,1,1] %>%
    melt() %>%
    group_by(Var2) %>% 
    mutate(sumP = sum(value), pbin = value/sumP)
  

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
png(here("figs","LAA_Dist_A1.png"), 
    height = 10, width = 10, unit = 'in', res = 420)
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

## catches post F tuning

## new option - by area weighted Fs
catch_yf_predt <- data.frame(matrix(NA, nrow = tEnd, ncol = nfleets_fish))
for(flt in 1:nfleets_fish){
  for(y in 1:(tEnd-1)){
  catch_yf_predt[y,flt] <- sum(catch_yfi_pred[y,flt,])
  }
}

## old option -- fleet -specific
catch_yf_predt <- data.frame(catch_yf_pred)[1:(tEnd-1),]

## comparing options
# plot(Freal_yf[,2] ~ rowSums(F_area_yfi[,2,]), 
#      xlab = 'New Method sum of F x fleet x Area', ylab = 'old method F_fleet')
# abline(0,1,col = 'red',add = TRUE)
# 
# plot(catch_yf_pred[,2] ~   catch_yf_predt[,2], 
#      xlab = 'New Method sum of F x fleet x Area', ylab = 'old method F_fleet')
# abline(0,1,col = 'red',add = TRUE)



names(catch_yf_predt) <- paste('Fleet',1:3)

catch_yf_predt <- catch_yf_predt %>%
  mutate(Year = 1:nrow(catch_yf_predt)) %>%
  melt(id = 'Year') %>%
  mutate(Type = 'PRED')

catch_yf_obst<- catch_yf_obs[1:52,]%>%  data.frame() %>%select(-X1) 
names(catch_yf_obst) <- paste('Fleet',1:3)

catch_yf_obst <- catch_yf_obst %>%
  mutate(Year = 1:52) %>%
  melt(id = 'Year') %>%
  ## convert CV to SD via CV = mean/sd
  mutate(Type = 'OBS', 
         lci = value - 1.96*(0.1*value),
         uci = value + 1.96*(0.1*value)) 

ggplot(data = catch_yf_obst, aes(x = Year, y = value, color = variable)) +
  geom_line(data = catch_yf_predt, lwd = 0.75) +
  scale_color_manual(values = c('blue','blue1','blue2')) +
  geom_point(pch = 1, fill = NA, col = 'black') +
  geom_errorbar(aes(ymin = lci, ymax = uci), col = 'black') +
  theme_sleek() + theme(legend.position = 'none')+
  labs(y = 'Catch (lbs)', color = 'Fishing Fleet') +
  facet_wrap(~variable, scales = "free_y")

ggsave(last_plot(),
       file = here('figs',
                   paste0('catch_fits_',
                          'v1=',v1,Sys.Date(),'.png')),
       width = 10, height = 6, unit = 'in',
       dpi = 420)




# 
# simdata$SSB_yk %>% 
#   data.frame() %>%
#   mutate('Yr' = 1:53) %>%
#   reshape2::melt(id = c('Yr')) %>%  
#   ggplot(., aes(x = Yr, y = value, color = variable )) +
#   geom_line(lwd = 2) + labs(x = 'Modeled Year',y = 'SSB', color = 'stock') +
#   ggsidekick::theme_sleek()

# }