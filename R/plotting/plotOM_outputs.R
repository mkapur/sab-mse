# require(ggplot2)
# require(patchwork)
# library(gridExtra)
# require(dplyr)
plotOM <- function(dat, dumpfile = here('figs',paste0(Sys.Date(),"/"))){
  attach(dat)
  if(!exists(dumpfile)) dir.create(dumpfile)
  
  years <- 1960:2019
  nyear <- length(years)
  tEnd <- length(years)
  age <- 0:70 
  nage <- length(age)
  
  ## ninit ----
  png(file =paste0(dumpfile,"/",Sys.Date(),'-Ninit_ais.png'),
      width = 10, height = 8, unit = 'in', res = 420)
  Ninit_ais[,,1] %>%
    data.frame() %>%
    mutate('Age' = age) %>%
    reshape2::melt(id = c('Age')) %>%
    ggplot(., aes(x = Age, y = value, color = variable )) +
    scale_color_manual(values = subareaPal) +
    geom_line(lwd = 2) + 
    labs(x = 'Age in Initial Years',y = 'Initial Numbers', color = 'subarea') +
    ggsidekick::theme_sleek()+
    facet_wrap(~variable,scales = 'free_y')
  dev.off()
  
  
  ## N_yseason ----
  png(file =paste0(dumpfile,"/",Sys.Date(),'-N_season_iy.png'),
      width = 10, height = 8, unit = 'in', res = 420)
  par(mfrow = c(2,3))
  for(i in 1:6){
    ylt = 10*max(sum(N_yais_end[2,,i,][!is.na(N_yais_end[2,,i,])]),
                 sum(N_yais_end[10,,i,][!is.na(N_yais_end[10,,i,])]))
    
    plot(rowSums(N_yais_beg[,,i,]),
         type = 'l',
         lwd = 2, 
         col = scales::alpha(subareaPal[i],0.2),
         main = inames[i], 
         col.main = subareaPal[i], 
         ylim = c(0,ylt),
         xlim = c(0,nyear),xaxt='n',
         xlab = "Model Year", 
         ylab = 'Numbers (M+F)')
    lines(rowSums(N_yais_mid[,,i,]),
          type = 'l',
          lwd = 3,
          col = scales::alpha(subareaPal[i],0.4))
    lines(rowSums(N_yais_end[,,i,]),
          type = 'l',
          lwd = 3,
          col = scales::alpha(subareaPal[i],0.8))
    legend("topright",col = c(scales::alpha(subareaPal[i],0.2),
                              scales::alpha(subareaPal[i],0.4),
                              scales::alpha(subareaPal[i],0.8)), 
           legend = c("beg",
                      "mid (move)",
                      "end (fished)"), cex = 0.7, lty =1, lwd = 5)
    axis(1, at = seq(1,nyear,5), labels = years[seq(1,nyear,5)])
  }
  dev.off()
  ## ssb ----
  SSB_yk %>%
    data.frame() %>%
    mutate('Yr' = years[1:nrow(.)]) %>%
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
  SSB_yi <- data.frame(SSB_yi)
  names(SSB_yi) <- spmat$subarea
  SSB_ym0 <-  SSB_yi %>%
    mutate('Yr' = years[1:nrow(.)]) %>%
    reshape2::melt(id = c('Yr')) %>%
    merge(., spmat, by.x = "variable", by.y = "subarea") %>% 
    group_by(Yr, mgmt) %>%
    summarise(totSSBkg = sum(value), 
              totSSBmt = totSSBkg) %>% select(Yr, mgmt, totSSBmt)
  levels(SSB_ym0$mgmt) <- c('AK','BC','CC','WC', 'AI')
  SSB_ym0$mgmt[SSB_ym0$mgmt == 'AI'] <- 'AK'  
  SSB_ym0$mgmt[SSB_ym0$mgmt == 'CC'] <- 'WC'  
  SSB_ym <- SSB_ym0 %>% group_by(Yr, mgmt) %>%
    summarise(omSSBMT = sum(totSSBmt))
  
  merge(assSB, SSB_ym, by.x = c('Year','REG'), by.y = c('Yr','mgmt')) %>%
    select(Year, REG, CV,assSSBMT, omSSBMT) %>%
    filter(Year >1965 ) %>%
    ggplot(., aes(x = Year, y = assSSBMT, color = REG)) +
    geom_line(aes(y = omSSBMT),lwd = 1.1) +
    geom_errorbar(aes(ymin = assSSBMT-CV*assSSBMT,ymax= assSSBMT+CV*assSSBMT, color = REG)) +
    scale_color_manual(values = mgmtPal)+
    geom_point()+  ggsidekick::theme_sleek() +
    labs(x = 'Modeled Year',y = 'SSB', color = 'Mgmt Region') +
    facet_wrap(~REG,scales = 'free_y')
  ggsave(last_plot(),
         file = paste0(dumpfile,'SSB_ym-',Sys.Date(),'.png'),
  width = 10, height = 6, unit = 'in',
  dpi = 420)
  
  ## catch pred by fleet ----
  catch_yf_predt <- data.frame(catch_yf_pred)
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
    scale_x_continuous(limits = c(1960,1960+df$yRun)) +
    geom_point(pch = 1, fill = NA, col = 'black') +
    geom_errorbar(aes(ymin = lci, ymax = uci), col = 'black') +
    theme_sleek() + 
    theme(legend.position = 'none')+
    labs(y = 'catch', color = 'Fishing Fleet')+
    facet_wrap(~variable, scales = "free_y")
  ggsave(last_plot(),
         file = paste0(dumpfile,'catch_fits_TMB_',
                            'v1=',df$v1,'niter=',df$niter,'Fmax=',df$Fmax,Sys.Date(),'.png'),
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
    scale_x_continuous(limits = c(1960,1960+df$yRun)) +
    geom_point(pch = 1, fill = NA, col = 'black') +
    geom_errorbar(aes(ymin = lci, ymax = uci), col = 'black') +
    theme_sleek() + 
    theme(legend.position = 'none')+
    labs(y = 'catch', color = 'Fishing Fleet')+
    facet_wrap(~REG, scales = "free_y")
  
  ggsave(last_plot(),
         file = paste0(dumpfile,'catchm_fits_TMB_',
                            'v1=',df$v1,'niter=',df$niter,'Fmax=',df$Fmax,Sys.Date(),'.png'),
         width = 10, height = 6, unit = 'in',
         dpi = 420)
  
  ## survey preds ----
  
  survey_yf_predt <- data.frame(surv_yf_pred)
  names(survey_yf_predt) <- df$fltnames_surv
  survey_yf_predt <- survey_yf_predt %>%
    mutate(Year = years) %>%
    melt(id = 'Year') %>%
    mutate(Type = 'PRED') %>%
    mutate(REG = substr(variable,0,2)) #%>%
  filter(REG == 'BC') #%>% View()
  
  
  survey_yf_obst <- data.frame( df$surv_yf_obs)  %>%
    mutate(Year = years) %>%
    melt(id = 'Year') %>%
    ## convert CV to SD via CV = mean/sd
    mutate(Type = 'OBS', 
           lci = value - 1.96*(0.1*value),
           uci = value + 1.96*(0.1*value)) %>%
    mutate(REG = substr(variable,0,2)) %>%
    filter(value > 0) #%>% View()
  
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
         file =paste0(dumpfile,'survey_fits_selMult_',
                      'v1=',v1,'Fmax=',Fmax,Sys.Date(),'.png'),
         width = 10, height = 6, unit = 'in',
         dpi = 420)
}
