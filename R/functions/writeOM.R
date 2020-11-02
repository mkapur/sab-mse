# require(ggplot2)
# require(patchwork)
# library(gridExtra)
# require(dplyr)
writeOM <- function(dat, opt, obj, 
                    rep =NA,
                    cppname = NA,
                    runname = NA,
                    dumpfile =  here('output',paste0(Sys.Date(),"/"))){
  ## allow for extra name
  if(!is.na(runname)) dumpfile = paste0(dumpfile,runname,"/")

  if(!exists(dumpfile)) dir.create(dumpfile)
  ## save the CPP used here
  
  cppfile <- list.files(here('TMB'),
                        pattern = paste0("*",cppname,"*.cpp"),
                        full.names = TRUE)
  file.copy(cppfile,
            to = paste0(dumpfile,paste0(basename(cppfile))))
  
  ## write DF used here
  save(df, file = paste0(dumpfile,"/dfUSED.rdata"))
  
  spmat <- data.frame(subarea = c('A1',"A3","B3","B2","C2","C1"),
                      stock = c("R4","R3","R3","R2","R2","R1"),
                      mgmt = c("AI","AK", rep("BC",2), rep("CC",2)))
  inames = rev(unique(spmat$subarea))
  
  ## save stuff into dumpfile
  save(dat, file = paste0(dumpfile,"/dat.rdata"))
  save(opt, file = paste0(dumpfile,"/opt.rdata"))
  save(obj, file = paste0(dumpfile,"/obj.rdata"))
  if(!is.na(rep))   save(rep, file = paste0(dumpfile,"/rep.rdata"))
  
  years <- 1960:2019
  nyear <- length(years)
  tEnd <- length(years)
  age <- 0:70 
  nage <- length(age)
  
  ## ninit ----
  png(file = paste0(dumpfile,'Ninit_ais.png'),
      width = 10, height = 8, unit = 'in', res = 420)
  ninit0 <-  dat$Ninit_ais[,,1] %>%
    data.frame() 
  names(ninit0) <- inames
  ninit0 %>%
    mutate('Age' = age) %>%
    reshape2::melt(id = c('Age')) %>%
  ggplot(., aes(x = Age, y = value, color = variable )) +
    scale_color_manual(values = rev(subareaPal)) +
    geom_line(lwd = 2) + 
    labs(x = 'Age in Initial Years',y = 'Initial Numbers', color = 'subarea') +
    ggsidekick::theme_sleek()+
    facet_wrap(~variable,scales = 'free_y')
  dev.off()
  
  
  ## N_yseason ----
  png(file =paste0(dumpfile,"/",
                   Sys.Date(),'-N_season_iy.png'),
      width = 10, height = 8, unit = 'in', res = 420)
  par(mfrow = c(2,3))
  for(i in 1:6){
    ylt = 2*max(sum(dat$N_yais_end[2,,i,][!is.na(dat$N_yais_end[2,,i,])]),
                 sum(dat$N_yais_end[df$yRun,,i,][!is.na(dat$N_yais_end[df$yRun,,i,])]))
    
    plot(rowSums(dat$N_yais_beg[1:(df$yRun-1),,i,]),
         type = 'l',
         lwd = 2, 
         col = scales::alpha(rev(rev(subareaPal)[i]),0.2),
         main = inames[i], 
         col.main = rev(rev(subareaPal)[i]), 
         ylim = c(0,ylt),
         xlim = c(0,nyear),xaxt='n',
         xlab = "Model Year", 
         ylab = 'Numbers (M+F)')
    lines(rowSums(dat$N_yais_mid[,,i,]),
          type = 'l',
          lwd = 3,
          col = scales::alpha(rev(rev(subareaPal)[i]),0.4))
    lines(rowSums(dat$N_yais_end[,,i,]),
          type = 'l',
          lwd = 3,
          col = scales::alpha(rev(rev(subareaPal)[i]),0.8))
    legend("topright",col = c(scales::alpha(rev(rev(subareaPal)[i]),0.2),
                              scales::alpha(rev(rev(subareaPal)[i]),0.4),
                              scales::alpha(rev(rev(subareaPal)[i]),0.8)), 
           legend = c("beg",
                      "mid (move)",
                      "end (fished)"), cex = 0.7, lty =1, lwd = 5)
    axis(1, at = seq(1,nyear,5), labels = years[seq(1,nyear,5)])
  }
  dev.off()
  ## ssb ----
  dat$SSB_yk %>%
    data.frame() %>%
    mutate('Yr' = years[1:nrow(.)]) %>%
    reshape2::melt(id = c('Yr')) %>%
    ggplot(., aes(x = Yr, y = value, color = variable )) +
    scale_color_manual(values = demPal) +
    geom_line(lwd = 2) + labs(x = 'Modeled Year',y = 'SSB', color = 'stock') +
    ggsidekick::theme_sleek() + 
    theme( legend.position = c(0.8,0.8))
  ggsave(last_plot(),
         file = paste0(dumpfile,'/SSB_yk-',Sys.Date(),'.png'),
         width = 10, height = 6, unit = 'in',
         dpi = 420)
  
  
  ## ssb_ym with compare ----
  spmat <- data.frame(subarea = c('A1',"A2","B2","B1","C2","C1"),
                      stock = c("R4","R3","R3","R2","R2","R1"),
                      mgmt = c("AK","AK", rep("BC",2), rep("CC",2)))
  
  assSB <- read.csv(here('input','downloads','AssessmentDat_thru2018.csv'),fileEncoding="UTF-8-BOM") %>%
    mutate(REG = substr(Index,1,2), assSSBMT = Value) %>% filter(Type == 'SpawnBiomass')
  SSB_yi <- data.frame(dat$SSB_yi)
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
         file = paste0(dumpfile,'/SSB_ym-',Sys.Date(),'.png'),
  width = 10, height = 6, unit = 'in',
  dpi = 420)
  
  ## catch pred by fleet ----
  catch_yf_predt <- data.frame(dat$catch_yf_pred)
  names(catch_yf_predt) <- df$fltnames_fish
  catch_yf_predt <- catch_yf_predt %>%
    mutate(Year = years) %>%
    group_by(Year) %>%
    mutate("AK_FIX (aggregate)" = sum(AK_FIX_W, AK_FIX_E),
           "AK_TWL (aggregate)"= sum(AK_TWL_W, AK_TWL_E)) %>%
    select(-AK_TWL_W,-AK_TWL_E,-AK_FIX_W,-AK_FIX_E) %>%
    melt(id = 'Year') %>%
    mutate(Type = 'PRED') %>%
    mutate(REG = substr(variable,0,2)) %>% filter(value != 0)

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
    mutate(REG = substr(variable,0,2)) %>% filter(value > -1)
  
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
    facet_wrap(~variable, scales = "free_y", ncol = 2)
  ggsave(last_plot(),
         file = paste0(dumpfile,'/catch_fits_TMB_',
                            'v1=',df$v1,'niter=',df$niter,'Fmax=',df$Fmax,Sys.Date(),'.png'),
         width = 8, height = 6, unit = 'in',
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
         file = paste0(dumpfile,'/catchm_fits_TMB_',
                            'v1=',df$v1,'niter=',df$niter,'Fmax=',df$Fmax,Sys.Date(),'.png'),
         width = 10, height = 6, unit = 'in', dpi = 420)
  
  ## survey preds ----
  
  survey_yf_predt <- data.frame(dat$surv_yf_pred)
  names(survey_yf_predt) <- df$fltnames_surv
  survey_yf_predt <- survey_yf_predt %>%
    mutate(Year = years) %>%
    melt(id = 'Year') %>%
    mutate(Type = 'PRED') %>%
    mutate(REG = substr(variable,0,2)) #%>%

  
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
         file =paste0(dumpfile,'/survey_fits_selMult_',
                      'v1=',df$v1,'Fmax=',df$Fmax,Sys.Date(),'.png'),
         width = 10, height = 6, unit = 'in',
         dpi = 420)
  
  # plot FISH selex ----
  ## bring estimates out and rearrange
  inputSel <- array(exp(obj$par[grep('fsh_slx',names(obj$par))]),
                    dim = dim(df$parms$log_fsh_slx_pars),
                    dimnames = dimnames(df$parms$log_fsh_slx_pars))
  
  selP <- array(exp(opt$par[grep('fsh_slx',names(opt$par))]),
                dim = dim(df$parms$log_fsh_slx_pars),
                dimnames = dimnames(df$parms$log_fsh_slx_pars))
  dimnames(selP)[[1]] <- paste(df$fltnames_fish)
  
  fsh_sel_afs <- array(NA, dim = c(df$nage,df$nfleets_fish,2),
                       dimnames = list(c(df$age),
                                       c(paste(df$fltnames_fish)),
                                       c('Fem','Mal')))
  ## function to take the estimated parameters
  ## and info about selType, selShape
  ## and spit out vector of selx@a or selx@L
  
  for(a in 1:df$nage){
    for(s in 1:2){
      for(fish_flt in 1:df$nfleets_fish){
        fsh_sel_afs[,fish_flt,s] <- getSelec2(sex = s,
                                               selP = selP,
                                               flt_idx = fish_flt,
                                               selType = df$selType_fish[fish_flt], 
                                               selShape = df$selShape_fish[fish_flt])
          
        }
      }
    }

  
  
  png(paste0(dumpfile,'/fishery_selex.png'),
      height = 8, width = 6, unit = 'in', res = 420)
  par(mfrow = c(3,3) )
  for(flt in 1:df$nfleets_fish){
    for(s in 1:2){
      tmp <- fsh_sel_afs[,flt,s]
      if(s == 1) plot(tmp, 
                      col = sexPal[1], 
                      type = 'l', 
                      lwd = 2, 
                      xlab = ifelse(df$selType_fish[flt] == 0,
                                    'Age','Length'), 
                      ylab = 'Selectivity',
                      lty = 1,
                      ylim = c(0,1), 
                      main = df$fltnames_fish[flt], xlim = c(0,75),
                      col.main  = c(rep(mgmtPal[1],4), rep(mgmtPal[2],3),rep(mgmtPal[3],2))[flt])
      box(which = 'plot', lty = 'solid', 
          col = c(rep(mgmtPal[1],4), rep(mgmtPal[2],3),rep(mgmtPal[3],2))[flt], 
          lwd = 2)
      if(s == 2) lines(tmp, col = sexPal[2], type = 'l', lty = 2, lwd = 2)
    }
  }
  dev.off()
  
  
  # plot SURV selex ----
  ## bring estimates out and rearrange
  nsurvmod = df$nfleets_surv+(df$nfleets_acomp-5) ## bc we got selex for acomp flts too
  inputSel <- array(exp(obj$par[grep('srv_slx',names(obj$par))]),
                    dim = dim(df$parms$log_srv_slx_pars),
                    dimnames = dimnames(df$parms$log_srv_slx_pars))
  
  selP <- array(exp(opt$par[grep('srv_slx',names(opt$par))]),
                dim = dim(df$parms$log_srv_slx_pars),
                dimnames = dimnames(df$parms$log_srv_slx_pars))
  dimnames(selP)[[1]] <- dimnames(df$parms$log_srv_slx_pars)[[1]]
  
  srv_sel_afs <- array(NA, dim =  c(df$nage,nsurvmod,2),
                       dimnames = list(c(df$age),
                                                  c( dimnames(df$parms$log_srv_slx_pars)[[1]]),
                                                  c('Fem','Mal')))
  ## function to take the estimated parameters
  ## and info about selType, selShape
  ## and spit out vector of selx@a or selx@L
  
  for(a in 1:df$nage){
    for(s in 1:2){
      for(surv_flt in 1:nsurvmod){
        srv_sel_afs[,surv_flt,s] <- getSelec2(sex = s,
                                              selP = selP,
                                              flt_idx = surv_flt,
                                              selType = df$selType_surv[surv_flt], 
                                              selShape = df$selShape_surv[surv_flt])
        
      }
    }
  }
  
  
  
  png(paste0(dumpfile,'/survey_selex.png'),
      height = 8, width = 6, unit = 'in', res = 420)
  par(mfrow = c(4,2) )
  for(flt in 1:nsurvmod){
    for(s in 1:2){
      tmp <- srv_sel_afs[,flt,s]
      if(s == 1) plot(tmp, 
                      col = sexPal[1], 
                      type = 'l', 
                      lwd = 2, 
                      xlab = ifelse(df$selType_surv[flt] == 0,
                                    'Age','Length'), 
                      ylab = 'Selectivity',
                      lty = 1,
                      ylim = c(0,1), 
                      main = dimnames(df$parms$log_srv_slx_pars)[[1]][flt], xlim = c(0,75),
                      col.main  = c(rep(mgmtPal[1],4), rep(mgmtPal[2],3),rep(mgmtPal[3],2))[flt])
      box(which = 'plot', lty = 'solid', 
          col = c(rep(mgmtPal[1],4), rep(mgmtPal[2],3),rep(mgmtPal[3],2))[flt], 
          lwd = 2)
      if(s == 2) lines(tmp, col = sexPal[2], type = 'l', lty = 2, lwd = 2)
    }
  }
  dev.off()
  
  
  
} ## end func
