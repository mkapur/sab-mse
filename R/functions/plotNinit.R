require(ggplot2)
# plotNinit <- function(reps){
  nspace <- dim(reps$Ninit_Aai)[3]
  Ninits <- reps$Ninit_ai
  Ninits <- reps2$Ninit_ai
  N_yai_beg <- reps$N_yai_beg
  data.frame(matrix(Ninits, ncol=nspace, byrow=FALSE))
  
  
  Ninits %>%
    data.frame() %>%
    mutate('Yr' = 1:21) %>%
    reshape2::melt(id = c('Yr')) %>%
    ggplot(., aes(x = Yr, y = value, color = variable )) +
    geom_line(lwd = 2) + labs(x = 'Initializing Year',y = 'Numbers', color = 'subarea') +
    ggsidekick::theme_sleek()

  Nzero = reps$N_0ai
  Nzero %>% data.frame() %>%
    mutate('Yr' = 1:21) %>%
    reshape2::melt(id = c('Yr')) %>%
    ggplot(., aes(x = Yr, y = value, color = variable )) +
    geom_line(lwd = 2) + labs(x = 'Initializing Year',y = 'Unfished Numbers', color = 'subarea') +
    ggsidekick::theme_sleek()
  
  reps$SSB_yi %>% 
    data.frame() %>%
    mutate('Yr' = 1:53) %>%
    reshape2::melt(id = c('Yr')) %>%  
    ggplot(., aes(x = Yr, y = value, color = variable )) +
    geom_line(lwd = 2) + labs(x = 'Modeled Year',y = 'SSB', color = 'stock') +
    ggsidekick::theme_sleek()
  
  reps$SSB_yk %>% 
    data.frame() %>%
    mutate('Yr' = 1:53) %>%
    reshape2::melt(id = c('Yr')) %>%  
    ggplot(., aes(x = Yr, y = value, color = variable )) +
    geom_line(lwd = 2) + labs(x = 'Modeled Year',y = 'SSB', color = 'stock') +
    ggsidekick::theme_sleek()
  simdata$SSB_yk %>% 
    data.frame() %>%
    mutate('Yr' = 1:53) %>%
    reshape2::melt(id = c('Yr')) %>%  
    ggplot(., aes(x = Yr, y = value, color = variable )) +
    geom_line(lwd = 2) + labs(x = 'Modeled Year',y = 'SSB', color = 'stock') +
    ggsidekick::theme_sleek()
  
# }