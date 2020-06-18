require(ggplot2)
require(patchwork)
## SANITY CHECK PLOTS
nspace <- dim(reps$Ninit_Aai)[3]


# data.frame(matrix(Ninits, ncol=nspace, byrow=FALSE))




pNinit <- reps$Ninit_ai %>%
  data.frame() %>%
  mutate('Yr' = 1:21) %>%
  reshape2::melt(id = c('Yr')) %>%
  ggplot(., aes(x = Yr, y = value, color = variable )) +
  geom_line(lwd = 2) + 
  labs(x = 'Age in Initial Years',y = 'Initial Numbers', color = 'subarea') +
  ggsidekick::theme_sleek()

Nzero = reps$N_0ai
pNzero <- Nzero %>% data.frame() %>%
  mutate('Yr' = 1:21) %>%
  reshape2::melt(id = c('Yr')) %>%
  ggplot(., aes(x = Yr, y = value, color = variable )) +
  geom_line(lwd = 2) + 
  labs(x = 'Age in Year 0',y = 'Unfished Numbers', color = 'subarea') +
  ggsidekick::theme_sleek()


pNage1 <- reps$N_yai_beg[,,1] %>%
data.frame() %>%
  mutate('Yr' = 1:nrow(.)) %>%
  reshape2::melt(id = c('Yr')) %>%
  mutate(age = as.numeric(substr(variable,2,2))-1) %>%
  ggplot(., aes(x = age, y = value,  group = Yr )) +
  # geom_point() +
  geom_boxplot() +
  # geom_line() +
  labs(x = 'Age in Year',y = 'Numbers',  title = "AREA1") +
  ggsidekick::theme_sleek() + theme(legend.position = 'none') 


pNage2 <- reps$N_yai_beg[,,2] %>%
  data.frame() %>%
  mutate('Yr' = 1:nrow(.)) %>%
  reshape2::melt(id = c('Yr')) %>%
  mutate(age = as.numeric(substr(variable,2,2))-1) %>%
  ggplot(., aes(x = age, y = value, group = Yr )) +
  # geom_point() + 
  # geom_line()  +
  geom_boxplot() +
  labs(x = 'Age in Year',y = 'Numbers', title = "AREA2") +
  ggsidekick::theme_sleek()+ theme(legend.position = 'none')

pSSByi <- reps$SSB_yi %>%
  data.frame() %>%
  mutate('Yr' = 1:53) %>%
  reshape2::melt(id = c('Yr')) %>%
  ggplot(., aes(x = Yr, y = value, color = variable )) +
  geom_line(lwd = 2) + labs(x = 'Modeled Year',y = 'SSB', color = 'subarea') +
  ggsidekick::theme_sleek()

pSSByk <- reps$SSB_yk %>%
  data.frame() %>%
  mutate('Yr' = 1:53) %>%
  reshape2::melt(id = c('Yr')) %>%
  ggplot(., aes(x = Yr, y = value, color = variable )) +
  geom_line(lwd = 2) + labs(x = 'Modeled Year',y = 'SSB', color = 'stock') +
  ggsidekick::theme_sleek()

pRyi <- reps$R_yi %>% data.frame() %>%
  mutate('Yr' = 1:nrow(.)) %>%
  reshape2::melt(id = c('Yr')) %>%
  ggplot(., aes(x = Yr, y = value, color = variable )) +
  geom_line(lwd = 2) + 
  labs(x = 'Model Year',y = 'Recruits', color = 'subarea') +
  ggsidekick::theme_sleek()

pRyk <- reps$R_yk %>% data.frame() %>%
  mutate('Yr' = 1:nrow(.)) %>%
  reshape2::melt(id = c('Yr')) %>%
  ggplot(., aes(x = Yr, y = value, color = variable )) +
  geom_line(lwd = 2) + 
  labs(x = 'Model Year',y = 'Recruits', color = 'stock') +
  ggsidekick::theme_sleek()


pSRR <- cbind(reps$R_yk, reps$SSB_yk) %>% data.frame() %>%
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



pLAA1 <- reps$Length_yai_beg[,,1] %>% data.frame() %>%
  mutate('Yr' = 1:nrow(.)) %>%
  reshape2::melt(id = c('Yr')) %>%
  ggplot(., aes(x = as.numeric(substr(variable,2,3))-1, y = value, color = Yr )) +
  geom_point() +
  labs(x = 'Age',y = 'Length', color = 'Year') +
  ggsidekick::theme_sleek()

pLAA2 <- reps$Length_yai_beg[,,2] %>% data.frame() %>%
  mutate('Yr' = 1:nrow(.)) %>%
  reshape2::melt(id = c('Yr')) %>%
  ggplot(., aes(x = as.numeric(substr(variable,2,3))-1, y = value, color = Yr )) +
  geom_point() +
  labs(x = 'Age',y = 'Length', color = 'Year') +
  ggsidekick::theme_sleek()


## distributions
library(gridExtra)
pA1 <- list() ## for area 1
for(i in 1:10){ #45:dim(reps$LengthAge_alyi_beg)[3]){ ## loop years
  a1 <- 
    reps$LengthAge_alyi_beg[,,i,2] %>%
    melt() %>%
    group_by(Var2) %>% 
    mutate(sumP = sum(value), pbin = value/sumP)
  pA1[[i]] <- ggplot(a1,aes(x = Var1, y = Var2, fill = pbin)) +
    geom_tile() +
    labs(x = 'age', y = 'len', 
         title = paste("year ",i), 
         subtitle= 'subarea 1') +
    theme_sleek()
  rm(a1)
  # p[[i]] <- qplot(1:10,10:1,main=i)
}

pA2 <- list() ## for area 1
for(i in 1:10){ #45:dim(reps$LengthAge_alyi_beg)[3]){ ## loop years
  a1 <- 
    reps$LengthAge_alyi_beg[,,i,2] %>%
    melt() %>%
    group_by(Var2) %>% 
    mutate(sumP = sum(value), pbin = value/sumP)
  pA2[[i]] <- ggplot(a1,aes(x = Var1, y = Var2, fill = pbin)) +
    geom_tile() +
    labs(x = 'age', y = 'len', 
         title = paste("year ",i), 
         subtitle= 'subarea 2') +
    theme_sleek()
  rm(a1)
  # p[[i]] <- qplot(1:10,10:1,main=i)
}
do.call(grid.arrange,pA1)
do.call(grid.arrange,pA2)


# 
# simdata$SSB_yk %>% 
#   data.frame() %>%
#   mutate('Yr' = 1:53) %>%
#   reshape2::melt(id = c('Yr')) %>%  
#   ggplot(., aes(x = Yr, y = value, color = variable )) +
#   geom_line(lwd = 2) + labs(x = 'Modeled Year',y = 'SSB', color = 'stock') +
#   ggsidekick::theme_sleek()

# }