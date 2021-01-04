# require(ggplot2)
# require(patchwork)
# library(gridExtra)
# require(dplyr)
writeOM <- function(dat, 
                    opt,
                    obj,
                    mappy,
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
  save(df, file = paste0(dumpfile,"/df_used.rdata"))
  save(mappy, file = paste0(dumpfile,"/mappy.rdata"))
  
  sink(paste0(dumpfile,"/bounds.txt"))
  print(bounds)
  sink()  
  
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
  
  ## create and save metadata ----
  ## still in dev.
  # sink(paste0(dumpfile,"metadata.txt"))
  # cat('years \t',nyear,'\n')
  # cat('age \t',nage,'\n')
  # cat( as.data.frame(do.call(rbind, mappy)),"\n")
  # sink()
  ## ninit ----
  neqnm <- matrix(dat$NeqnR, ncol = 6, nrow = 71) %>%
    data.frame(.) %>%
    mutate(age = 0:70) %>%
    reshape2::melt(id = 'age')
  
  conservation_status <- c(
    X1 =  dimnames(df$X_ijas)[[1]][1],
    X2 = dimnames(df$X_ijas)[[1]][2],
    X3 = dimnames(df$X_ijas)[[1]][3],
    X4 =dimnames(df$X_ijas)[[1]][4],
    X5 = dimnames(df$X_ijas)[[1]][5],
    X6 = dimnames(df$X_ijas)[[1]][6]
  )
  
  ggplot(neqnm, aes(x = age, y = value, color = variable)) +
    ggsidekick::theme_sleek() +
    geom_line(lwd = 1.1) +
    scale_color_manual(values = rev(subareaPal),labels =  dimnames(df$X_ijas)[[1]]) +
    facet_wrap(~variable, scales = 'free_y', 
               labeller = labeller(variable =conservation_status)) +
    labs(x = 'Age',y = 'Numbers (M+F)', color = 'Subarea')
  
  ggsave(last_plot(),
         file =paste0(dumpfile,"/",
                      Sys.Date(),'-Neqn.png'),
         width = 8, height = 6, units = 'in', dpi = 520)
  
  ## N_yseason ----
  png(file =paste0(dumpfile,"/",
                   Sys.Date(),'-N_season_iy.png'),
      width = 12, height = 10, unit = 'in', res = 540)
  par(mfrow = c(2,3))
  for(i in 1:6){
    ylt = 20*max(sum(dat$N_yais_end[2,,i,][!is.na(dat$N_yais_end[2,,i,])]),
                sum(dat$N_yais_end[df$yRun,,i,][!is.na(dat$N_yais_end[df$yRun,,i,])]))
    ylt = 1.3*max(rowSums(dat$N_yais_beg[1:(df$yRun),,i,]))
    plot(rowSums(dat$N_yais_beg[1:(df$yRun),,i,]),
         type = 'l',
         lwd = 2, 
         col = scales::alpha(subareaPal[i],0.2),
         main = inames[i], 
         col.main =subareaPal[i], 
         ylim = c(0,ylt),
         xlim = c(0,df$yRun-1),xaxt='n',
         xlab = "Model Year", 
         ylab = 'Numbers (M+F)')
    lines(rowSums(dat$N_yais_mid[1:(df$yRun),,i,]),
          type = 'l',
          lwd = 3,
          col = scales::alpha(subareaPal[i],0.4))
    lines(rowSums(dat$N_yais_end[1:(df$yRun),,i,]),
          type = 'l',
          lwd = 3,
          col = scales::alpha(subareaPal[i],0.8))
    legend("topright",col = c(scales::alpha(subareaPal[i],0.2),
                              scales::alpha(subareaPal[i],0.4),
                              scales::alpha(subareaPal[i],0.8)), 
           legend = c("beg",
                      "mid (move + 1/2 fishing + 1/2 M)",
                      "end (1/2 fishing + 1/2 M)"), lty =1, lwd = 5)
    axis(1, at = seq(1,df$yRun,5), labels = years[seq(1,df$yRun,5)])
  }
  dev.off()
  
  ## ssb_yk ----
  dat$SSB_yk %>%
    data.frame() %>%
    mutate('Yr' = years[1:nrow(.)]) %>%
    reshape2::melt(id = c('Yr')) %>%
    ggplot(., aes(x = Yr, y = value, color = variable )) +
    scale_color_manual(values = demPal) +
    geom_line(lwd = 2) + labs(x = 'Modeled Year',y = 'SSB', color = 'stock') +
    scale_x_continuous(limits = c(1960,1959+df$yRun)) +
    ggsidekick::theme_sleek() + 
    theme( legend.position = c(0.8,0.8))
  
  ggsave(last_plot(),
         file = paste0(dumpfile,"/", Sys.Date(),'-SSB_yk.png'),
         width = 10, height = 6, unit = 'in',
         dpi = 420)
  
  ## ssb_yk gganimate ----
  # library(gganimate)
  # library(gifski)
  # myPlot <-dat$SSB_yk %>%
  #   data.frame() %>%
  #   mutate('Yr' = years[1:nrow(.)]) %>%
  #   reshape2::melt(id = c('Yr')) %>%
  #   ggplot(., aes(x = Yr, y = value, color = variable )) +
  #   scale_color_manual(values = demPal) +
  #   geom_line(lwd = 2) + 
  #   labs(x = 'Modeled Year',y = 'SSB', color = 'stock') +
  #   ggsidekick::theme_sleek() + 
  #   theme( legend.position = c(0.8,0.8))+ 
  #   transition_reveal(Yr)
  # 
  # animate(myPlot, duration = 10, fps = 20, width = 6,
  #         height = 6, unit = 'in', res = 420, renderer = gifski_renderer())
  # anim_save(paste0(dumpfile,'/SSB_yk-animate-',Sys.Date(),'.png'))
  
  # anim_save(paste0(dumpfile,"SSB_yk_animate.gif",animation = last_animation())) 
  ## ssb_ym with compare ----
  spmat <- data.frame(subarea = c('A1',"A3","B3","B2","C2","C1"),
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
    mutate('Last Assessment' = assSSBMT, 'Operating Model' =omSSBMT) %>%
    select(Year, REG, CV, 'Last Assessment', 'Operating Model') %>%
    # filter(Year >1965 ) %>%
    melt(id = c('Year','REG', 'CV')) %>%
    ggplot(., aes(x = Year, y = value, color = REG)) +
    # ggplot(., aes(x = Year, y = assSSBMT, color = REG)) +
    geom_line(aes(y = value),lwd = 1.1) +
    # geom_errorbar(aes(ymin = assSSBMT-CV*assSSBMT,ymax= assSSBMT+CV*assSSBMT, color = REG,width=0)) +
    geom_errorbar(aes(ymin = value-CV*value,ymax= value+CV*value, color = REG,width=0)) +
    scale_color_manual(values = mgmtPal)+
    geom_point()+ 
    ggsidekick::theme_sleek() +
    scale_x_continuous(limits = c(1960,1959+df$yRun)) +
    labs(x = 'Modeled Year',y = 'SSB (units vary)', color = 'Mgmt Region') +
    facet_wrap(~REG+variable,scales = 'free_y', ncol = 2)
    # facet_wrap(~REG,scales = 'free_y')
  ggsave(last_plot(),
         file = paste0(dumpfile,"/", Sys.Date(),'-SSB_ym.png'),
         width = 10, height = 6, unit = 'in',
         dpi = 420)
  
  
  ## plot SRR----
  dat$R_yk %>% 
    data.frame() %>%
    mutate('Yr' = 1:nrow(.))  %>%
    reshape2::melt(.,id = c('Yr')) %>%
    mutate(RYK = value) %>%
    select(-value) %>%
    bind_cols(.,
              dat$SSB_yk %>% 
                data.frame() %>%
                mutate('Yr' = 1:nrow(.))  %>%
                reshape2::melt(.,id = c('Yr')) %>%
                mutate(SSByk = value) %>%
                select(-value,-variable)) %>%
    ggplot(., aes(x = SSByk, y = RYK, color = variable, group = Yr...1)) +
    scale_color_manual(values = demPal) +
    geom_point() +
    labs(x = 'SSB',y = 'Recruits #', color = 'stock') +
    # scale_y_continuous(limits = c(0,25000)) +
    # scale_x_continuous(limits = c(0,400000)) +
    ggsidekick::theme_sleek() +
    facet_wrap(~ variable, scales = 'free')
  
  ## catch pred by fleet ----
  catch_yf_pred_totalt <- data.frame(dat$catch_yf_pred_total)
  names(catch_yf_pred_totalt) <- df$fltnames_fish
  catch_yf_pred_totalt <- catch_yf_pred_totalt %>%
    mutate(Year = years) %>%
    group_by(Year) %>%
    melt(id = 'Year') %>%
    mutate(Type = 'PRED') %>%
    mutate(REG = substr(variable,0,2)) %>% 
    filter(value != 0  | variable == 'AK_FIX')
  
  catch_yf_obst <- df$catch_yf_obs %>%  data.frame() %>%select(-Year) 
  
  catch_yf_obst <- catch_yf_obst %>%
    mutate(Year = years) %>%
    group_by(Year) %>%
    # mutate("AK_FIX (aggregate)" = sum(AK_FIX_W, AK_FIX_E),
    #        "AK_TWL (aggregate)"= sum(AK_TWL_W, AK_TWL_E)) %>%
    # select(-AK_TWL_W,-AK_TWL_E,-AK_FIX_W,-AK_FIX_E) %>%
    melt(id = 'Year') %>%
    ## convert CV to SD via CV = mean/sd
    mutate(Type = 'OBS', 
           lci = value - 1.96*(0.1*value),
           uci = value + 1.96*(0.1*value)) %>%
    mutate(REG = substr(variable,0,2)) %>% 
    filter(value > -1)
  
  ggplot(data = catch_yf_obst, 
         aes(x = Year, y = value, color = variable)) +
    geom_line(data = catch_yf_pred_totalt, lwd = 0.75) +
    scale_color_manual(values = fishfltPal) +
    scale_x_continuous(limits = c(1960,1959+df$yRun)) +
    geom_point(pch = 1, fill = NA, col = 'black') +
    geom_errorbar(aes(ymin = lci, ymax = uci), col = 'black',width=0) +
    theme_sleek() + 
    theme(legend.position = 'none')+
    labs(y = 'catch', color = 'Fishing Fleet')+
    facet_wrap(~variable, scales = "free_y", ncol = 2)
  ggsave(last_plot(),
         file = paste0(dumpfile,"/", Sys.Date(),'-catch_fits.png'),
         width = 8, height = 6, unit = 'in',
         dpi = 420)
  
  
  ## catch pred by m ----
  catch_yf_pred_totalm <- catch_yf_pred_totalt %>% 
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
    geom_line(data = catch_yf_pred_totalm, lwd = 1.1) +
    scale_color_manual(values = mgmtPal) +
    scale_x_continuous(limits = c(1960,1959+df$yRun)) +
    geom_point(pch = 1, fill = NA, col = 'black') +
    geom_errorbar(aes(ymin = lci, ymax = uci), col = 'black',width=0) +
    theme_sleek() + 
    theme(legend.position = 'none')+
    labs(y = 'catch', color = 'Fishing Fleet')+
    facet_wrap(~REG, scales = "free_y")
  
  ggsave(last_plot(),
         file = paste0(dumpfile,"/", Sys.Date(),'-catchM_fits.png'),
         width = 10, height = 6, unit = 'in', dpi = 420)
  
  ## survey preds ----
  
  survey_yf_predt <- data.frame(dat$surv_yf_pred)
  names(survey_yf_predt) <- df$fltnames_surv
  survey_yf_predt <- survey_yf_predt %>%
    mutate(Year = years) %>%
    melt(id = 'Year') %>%
    mutate(Type = 'PRED') %>%
    mutate(REG = substr(variable,0,2)) %>% 
    filter(value > 0)
  
  survey_yf_errt <- data.frame(df$surv_yf_err) %>%    
    mutate(Year = years) %>%
    melt(id = 'Year') %>%
    filter(!is.na(value)) 
  
  
  survey_yf_obst <- data.frame( df$surv_yf_obs)  %>%
    mutate(Year = years) %>%
    melt(id = 'Year') %>%
    filter(value > 0) %>%
    merge(survey_yf_errt, by = c('Year','variable')) %>%
    mutate(value = value.x) %>%
    ## convert CV to SD via CV = mean/sd
    ## we think these are sds in mt (input to model as log)
    mutate(Type = 'OBS',
           lci = round(value - value*(value.y)),
           uci = round(value + (value*value.y))) %>%
    mutate(REG = substr(variable,0,2)) %>%
      select(Year, variable, value, lci, uci)
  
  ggplot(data = survey_yf_obst, 
         aes(x = Year, y = value, color = variable)) +
    geom_line(data = survey_yf_predt, lwd = 0.75) +
    scale_color_manual(values = survfltPal) +
    geom_point(pch = 1, fill = NA, col = 'black') +
    geom_errorbar(aes(ymin = lci, ymax = uci), col = 'black',width=0) +
    scale_x_continuous(limits = c(1980,ifelse(df$yRun == 59,2021,1959+df$yRun)), 
                       breaks = seq(1980,1960+df$yRun,10),
                       labels = seq(1980,1960+df$yRun,10)) +
    theme_sleek() + 
    theme(legend.position = 'none')+
    labs(y = 'survey', color = 'Fishing Fleet')+
    facet_wrap(~variable, scales = "free_y")
  
  ggsave(last_plot(),
         file =paste0(dumpfile,"/", Sys.Date(),'-survey_fits.png'),
         width = 10, height = 6, unit = 'in',
         dpi = 420)
  
  # plot lengths at age ----
  
  
  # plot FISH selex ----
  ## bring estimates out and rearrange
  
  ## if everythign fixed use
  if(length(mappy$log_fsh_slx_pars) == length(df$parms$log_fsh_slx_pars) ){
    inputSel <- array(exp(df$parms$log_fsh_slx_pars),
                      dim = dim(df$parms$log_fsh_slx_pars),
                      dimnames = dimnames(df$parms$log_fsh_slx_pars))
    
    selP <- array(exp(df$parms$log_fsh_slx_pars),
                  dim = c(7,2,1,2))
  } else{
    inputSel <- array(exp(obj$par[grep('fsh_slx',names(obj$par))]),
                      dim = dim(df$parms$log_fsh_slx_pars),
                      dimnames = dimnames(df$parms$log_fsh_slx_pars))
    
    selP <- array(exp(opt$par[grep('fsh_slx',names(opt$par))]),
                  dim = c(7,2,1,2))
  }
  
  mapped_fsh_selnames <- c('AK_FIX','AK_TWL',paste(df$fltnames_fish[3:df$nfleets_fish]))

  dimnames(selP)[[1]] <- mapped_fsh_selnames
  
  fsh_sel_afs <- array(NA, dim = c(df$nage,
                                   length(mapped_fsh_selnames),#df$nfleets_fish
                                   2),
                       dimnames = list(c(df$age),
                                       mapped_fsh_selnames,# c(paste(df$fltnames_fish)),
                                       c('Fem','Mal')))
  ## function to take the estimated parameters
  ## and info about selType, selShape
  ## and spit out vector of selx@a or selx@L
  
  for(a in 1:df$nage){
    for(s in 1:2){
      for(fish_flt in 1:length(mapped_fsh_selnames)){
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
  par(mfrow = c(4,2) )
  for(flt in 1:length(mapped_fsh_selnames)){
    for(s in 1:2){
      tmp <- fsh_sel_afs[,flt,s]
      
      ## if fixed overwrite colors
      if(length(mappy$log_fsh_slx_pars) == length(df$parms$log_fsh_slx_pars) ){
        sexPal_temp = c('grey22','grey66')
      } else{
        sexPal_temp = sexPal
      }
      
      if(s == 1) plot(tmp, 
                      col = sexPal_temp[1], 
                      type = 'l', 
                      yaxt = 'n',
                      lwd = 2, 
                      xlab = ifelse(df$selType_fish[flt] == 0,
                                    'Age','Length'), 
                      ylab = 'Selectivity',
                      lty = 1,
                      ylim = c(0,1), 
                      main = mapped_fsh_selnames[flt], xlim = c(0,75),
                      col.main  = c(rep(mgmtPal[1],2), 
                                    rep(mgmtPal[2],3),
                                    rep(mgmtPal[3],2))[flt])
      
      box(which = 'plot', lty = 'solid', 
          col = c(rep(mgmtPal[1],2), rep(mgmtPal[2],3),rep(mgmtPal[3],2))[flt], 
          lwd = 2)
      axis(2, at = seq(0,1,0.2), labels = seq(0,1,0.2))
      if(s == 2) lines(tmp, col = sexPal_temp[2], type = 'l', lty = 2, lwd = 2)
    }
  }
  dev.off()
  
  
  # plot SURV selex ----
  ## bring estimates out and rearrange
  nsurvmod = df$nfleets_surv+(df$nfleets_acomp-4) ## bc we got selex for acomp flts too
  
  
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
                                              selShape = df$selShape_surv[surv_flt],
                                              fltType = 'surv')
        
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
