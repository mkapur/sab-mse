## one-off make HCR schematic plots
require(here)
require(dplyr)
require(ggplot2)
require(ggsidekick)
mgmtPal <- c( "#015b58",  "#2c6184", "#984e73")

## hcr function ----
hcr <- function(biomass, depl, B0 =1000,rule = 1){
  if(rule == 1){
    ## PFMC 40-10 rule: line between B10 and B40, biomass*Fmsy thereafter
    FMSY = 0.5 ## where to anchor the upper limit reference point catch (for plotting)
    ulim = 0.4; llim =0.1
    FFt <- ifelse(  depl<= ulim & depl >= llim,
             FMSY/(ulim -llim)*depl-FMSY*llim/(ulim -llim),
             ifelse( biomass/(B0*ulim) > 1, FMSY, 0))
  } else if(rule == 2){
    ## bc mse 60/40
    FMSY = 0.5 
    ulim = 0.6; llim = 0.4
    FFt <- ifelse( depl/ulim <= 1 &   depl/ulim >= llim,
             FMSY*( depl/ulim-llim)/(1-llim),
             ifelse( depl /ulim > 1, FMSY, 0))
  } else if(rule == 3){
    ## npfmc tier 3
    ## ak upper is b40 or 40% unfished spawning biomass (see p 43 of recent assessment)
    FMSY = 0.5*0.95 ## actually F35
    ulim = 0.4; llim = 0.05 
    FFt <- ifelse(  biomass/(B0*ulim) <= 1 & biomass/(B0*ulim) >= llim,
             FMSY*(biomass/(B0*ulim)-llim)/(1-llim), 
             ifelse( biomass/(B0*ulim) > 1, FMSY, 0))
  } else if(rule == 4){
    ## constant catch
    FFt <- FMSY * 0.5
  } ## end rule 4
  return(FFt)
}

B0 = 1000 ## constant
bs = seq(0,1000,1)

## build wc df ----
wc <- data.frame('REG' ='PFMC (California Current)','B' = NA, 'DEPL' = NA,'FF_sq' = NA,'FF_NPFMC' = NA, 'CATCH' = NA)

# B40 = B0*0.4
for(b in seq_along(bs)){
  wc[b,'REG'] <- 'PFMC (California Current)'
  wc[b,'B'] <- bs[b]
  wc[b,'DEPL'] <-  wc[b,'B']/B0
  ## F rate is relative to depletion level not straight up B
  wc[b,'FF_sq'] <-  hcr(biomass =wc[b,'B'], depl = wc[b,'DEPL'], rule = 1)
  # cat(wc[b,'B'],wc[b,'DEPL'],wc[b,'FF_sq'] ,"\n")
  wc[b,'FF_NPFMC'] <-  hcr(biomass =wc[b,'B'], depl = wc[b,'DEPL'], rule = 3)
  
  # wc[b,'FF'] <-  ifelse(  wc[b,'DEPL'] <= ulim & wc[b,'DEPL'] >= llim,
  #                         FMSY/0.3*wc[b,'DEPL']-FMSY*0.1/0.3,
  #                         # FMSY*(  wc[b,'DEPL']/ulim-llim)/(1-llim),
  #                         ifelse(   wc[b,'B']/B40 > 1, FMSY, 0))
  wc[b,'CATCH'] = wc[b,'FF_sq']*  wc[b,'B'] 
}

## build bc df ----
bc <- data.frame('REG' ='BC Harvest Rule','B' = NA,'DEPL' = NA,'FF_sq' = NA, 'FF_NPFMC' = NA,'CATCH' = NA)

for(b in seq_along(bs)){
  bc[b,'REG'] <- 'BC Harvest Rule'
  bc[b,'B'] <-  bs[b]
  bc[b,'DEPL'] <-  bc[b,'B']/B0
  bc[b,'FF_sq'] <- hcr( b =bc[b,'B'], depl = bc[b,'DEPL'], rule = 2)
  bc[b,'FF_NPFMC'] <- hcr( b =bc[b,'B'], depl = bc[b,'DEPL'], rule = 3)
  # bc[b,'FF'] <-  ifelse(  bc[b,'DEPL']/ulim <= 1 &   bc[b,'DEPL']/ulim >= llim,
  #                          FMSY*(  bc[b,'DEPL']/ulim-llim)/(1-llim),
  #                         ifelse(  bc[b,'DEPL'] /ulim > 1, FMSY, 0))
  bc[b,'CATCH'] <-  bc[b,'FF_sq']* bc[b,'B']
}
## build ak df ----
ak <- data.frame('REG' ='NPFMC (Alaska)','B' = NA, 'DEPL' = NA,'FF_sq' = NA,'FF_NPFMC' = NA, 'CATCH' = NA)
for(b in seq_along(bs)){
  ak[b,'REG'] <- 'NPFMC (Alaska)'
  ak[b,'B'] <-  bs[b]
  ak[b,'DEPL'] <-  ak[b,'B']/B0
  ## NOTE THAT LLIM MEANS 5% OF B40, NOT B5%
  ak[b,'FF_sq'] <-  hcr( b =ak[b,'B'], depl = ak[b,'DEPL'], rule = 3)
  ak[b,'FF_NPFMC'] <-   ak[b,'FF_sq']
  # ak[b,'FF'] <- ifelse(  ak[b,'B']/B40 <= 1 & ak[b,'B']/B40 >= llim,
  #                         FMSY*(ak[b,'B']/B40-llim)/(1-llim), 
  #                        ifelse( ak[b,'B']/B40 > 1, FMSY, 0))
  ak[b,'CATCH'] <-  ak[b,'FF_sq']* ak[b,'B']
}

## plot status quo ----
supp.labs <-  c('AK (status quo)','BC (status quo)', 'Cal Current(status quo)')
names(supp.labs) <- c('NPFMC (Alaska)','BC Harvest Rule', 'PFMC (California Current)')

rbind(wc,ak,bc) %>%
  mutate(REG = factor(REG, levels = c('NPFMC (Alaska)','BC Harvest Rule', 'PFMC (California Current)'))) %>%
  ggplot(., aes(x = DEPL, y = FF_sq, color = REG)) +
  geom_line(lwd = 3) +
  scale_color_manual(values = mgmtPal)+
  scale_y_continuous(limits = c(0,0.6))+
  theme_sleek() +theme(axis.text = element_blank(), axis.title = element_blank())+
  facet_wrap(~REG, labeller = labeller(REG = supp.labs) ) +
  theme(legend.position = 'none',
        strip.text = element_text(face="bold", size=22)) 
ggsave(here::here("figs", "HCR_statusquo.png"), 
       width = 15, height = 6)

## plot all NPFMC ----
supp.labs <-  c('AK (using NPFMC)','BC (using NPFMC)', 'Cal Current (using NPFMC)')
names(supp.labs) <- c('NPFMC (Alaska)','BC Harvest Rule', 'PFMC (California Current)')

rbind(wc,ak,bc) %>%
  mutate(REG = factor(REG, levels = c('NPFMC (Alaska)','BC Harvest Rule', 'PFMC (California Current)'))) %>%
  ggplot(., aes(x = DEPL, y = FF_NPFMC, color = REG)) +
  geom_line(lwd = 3) +
  scale_color_manual(values = mgmtPal)+
  theme_sleek() +theme(axis.text = element_blank(), axis.title = element_blank())+
  scale_y_continuous(limits = c(0,0.6))+
  facet_wrap(~REG, labeller = labeller(REG = supp.labs) ) +
  theme(legend.position = 'none',
        strip.text = element_text(face="bold", size=22))

ggsave(here::here("figs", "HCR_NPFMC.png"), 
       width = 15, height = 6)







rbind(wc,ak,bc) %>%
  ggplot(., aes(x = DEPL, y = CATCH, color = REG)) +
  geom_line(lwd = 2) +
  scale_color_manual(values = mgmtPal)+
  # scale_y_continuous(limits = c(0,0.8))+
  # scale_x_continuous(limits = c(0,0.7))+
  # theme_sleek() +theme(axis.text = element_blank(), axis.title = element_blank())+
  ggthemes::theme_solid()+
  theme(legend.position = 'none') +
  facet_wrap(~factor(REG, levels = c('NPFMC (Alaska)','BC Harvest Rule', 'PFMC (California Current)')),
             nrow = 1)
