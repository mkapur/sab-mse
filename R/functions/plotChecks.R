require(ggplot2)
require(patchwork)
library(gridExtra)

## SANITY CHECK PLOTS
nspace <- dim(Ninit_Aai)[3]


# data.frame(matrix(Ninits, ncol=nspace, byrow=FALSE))




pNinit <- Ninit_ai %>%
  data.frame() %>%
  mutate('Yr' = 1:21) %>%
  reshape2::melt(id = c('Yr')) %>%
  ggplot(., aes(x = Yr, y = value, color = variable )) +
  geom_line(lwd = 2) + 
  labs(x = 'Age in Initial Years',y = 'Initial Numbers', color = 'subarea') +
  ggsidekick::theme_sleek()




pNzero <- Nzero %>% data.frame() %>%
  mutate('Yr' = 1:21) %>%
  reshape2::melt(id = c('Yr')) %>%
  ggplot(., aes(x = Yr, y = value, color = variable )) +
  geom_line(lwd = 2) + 
  labs(x = 'Age in Year 0',y = 'Unfished Numbers', color = 'subarea') +
  ggsidekick::theme_sleek()

ggsave(pNzero,
       file = here('figs',
                   "NZero_noOmega.png"),
       width = 6, height = 4, unit = 'in',
       dpi = 420)

pNage1 <- N_yai_beg[,,1] %>%
data.frame() %>%
  mutate('Yr' = 1:nrow(.)) %>%
  reshape2::melt(id = c('Yr')) %>%
  mutate(age = as.numeric(substr(variable,2,2))-1) %>%
  ggplot(., aes(x = age, y = value,  group = Yr )) +
  # geom_point() +
  # geom_boxplot() +
  geom_line() +
  labs(x = 'Age in Year',y = 'Numbers',  title = "AREA1") +
  ggsidekick::theme_sleek() + theme(legend.position = 'none') 


plot(N_yai_beg[1,,1],type = 'l') #, ylim = c(0,20e5))
for(i in 1:10){
  lines(N_yai_beg[i,,1], col = grey((i+0.2)/10), lwd =2)
}

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

pSSByk <- SSB_yk[1:52,] %>%
  data.frame() %>%
  mutate('Yr' = 1:52) %>%
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


pSRR <- cbind(R_yk, SSB_yk) %>% data.frame() %>%
  mutate('Yr' = 1:nrow(.), RECS1 = X1, RECS2 = X2, SSBS1 = X3, SSBS2 = X4) %>%
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



pLAA1 <- Length_yai_beg[,,1] %>% data.frame() %>%
  mutate('Yr' = 1:nrow(.)) %>%
  reshape2::melt(id = c('Yr')) %>%
  ggplot(., aes(x = as.numeric(substr(variable,2,3))-1, y = value, color = Yr )) +
  geom_point() +
  labs(x = 'Age',y = 'Length', color = 'Year') +
  ggsidekick::theme_sleek()

pLAA2 <- Length_yai_beg[,,2] %>% data.frame() %>%
  mutate('Yr' = 1:nrow(.)) %>%
  reshape2::melt(id = c('Yr')) %>%
  ggplot(., aes(x = as.numeric(substr(variable,2,3))-1, y = value, color = Yr )) +
  geom_point() +
  labs(x = 'Age',y = 'Length', color = 'Year') +
  ggsidekick::theme_sleek()


## distributions

pA1 <- list() ## for area 1
for(y in 1:10){ #45:dim(LengthAge_alyi_beg)[3]){ ## loop years
  a1 <- 
    LengthAge_alyi_beg[,,y,1] %>%
    melt() %>%
    group_by(Var2) %>% 
    mutate(sumP = sum(value), pbin = value/sumP)
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
do.call(grid.arrange,pA1)
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
catch_yf_predt <- data.frame(catch_yf_pred)[1:(tEnd-1),] 
names(catch_yf_predt) <- paste('Fleet',1:3)

catch_yf_predt <- catch_yf_predt %>%
  mutate(Year = 1:52) %>%
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
                   'v1=',v1,'.png')),
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