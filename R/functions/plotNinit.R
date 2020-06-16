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
  ggplot(., aes(x = age, y = value, color = factor(Yr), group = Yr )) +
  geom_point() + geom_line() +
  labs(x = 'Age in Year',y = 'Numbers',  title = "AREA1") +
  ggsidekick::theme_sleek() + theme(legend.position = 'none') 

pNage2 <- reps$N_yai_beg[,,2] %>%
  data.frame() %>%
  mutate('Yr' = 1:nrow(.)) %>%
  reshape2::melt(id = c('Yr')) %>%
  mutate(age = as.numeric(substr(variable,2,2))-1) %>%
  ggplot(., aes(x = age, y = value, color = factor(Yr), group = Yr )) +
  geom_point() + geom_line()  +
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

(pNinit  | pNzero)/(pNage1  | pNage2)
(pRyi | pRyk)/(pSSByi  | pSSByk)


# 

# 
# simdata$SSB_yk %>% 
#   data.frame() %>%
#   mutate('Yr' = 1:53) %>%
#   reshape2::melt(id = c('Yr')) %>%  
#   ggplot(., aes(x = Yr, y = value, color = variable )) +
#   geom_line(lwd = 2) + labs(x = 'Modeled Year',y = 'SSB', color = 'stock') +
#   ggsidekick::theme_sleek()

# }