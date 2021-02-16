# require(ggplot2)
# require(patchwork)
# library(gridExtra)
# require(dplyr)
writeOM <- function(justPlots = FALSE,
                    dat, 
                    opt,
                    obj,
                    mappy,
                    rep =NA,
                    cppname = NA,
                    runname = NA,
                    dumpfile =  here('output',paste0(gsub(":","_",Sys.time()),"/"))){
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
  
  save(bounds, file = paste0(dumpfile,"/bounds.rdata"))
  
  # sink(paste0(dumpfile,"/bounds.txt"))
  # print(bounds)
  # sink()  
  
  sink(paste0(dumpfile,"/opt_report.txt"))
  print(opt)
  sink()
  
  spmat <- data.frame(subarea = c('A4',"A3","B3","B2","C2","C1"),
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
  tEnd <- df$tEnd
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
  
  neqm_plot <- ggplot(neqnm, aes(x = age, y = value, color = variable)) +
    ggsidekick::theme_sleek() +
    geom_line(lwd = 1.1) +
    scale_color_manual(values = rev(subareaPal),labels =  dimnames(df$X_ijas)[[1]]) +
    facet_wrap(~variable, scales = 'free_y', 
               labeller = labeller(variable =conservation_status)) +
    labs(x = 'Age',y = 'Numbers (M+F)', color = 'Subarea')
  
  if(justPlots) neqm_plot
  if(!justPlots) ggsave(neqm_plot,
                        file =paste0(dumpfile,"/",
                                     Sys.Date(),'-Neqn.png'),
                        width = 8, height = 6, units = 'in', dpi = 520)
  
  ## N_yseason ----
  if(!justPlots)  png(file =paste0(dumpfile,"/",
                                   Sys.Date(),'-N_season_iy.png'),
                      width = 12, height = 10, unit = 'in', res = 540)
  par(mfrow = c(2,3))
  for(i in 1:6){
    ylt = 20*max(sum(dat$N_yais_end[2,,i,][!is.na(dat$N_yais_end[2,,i,])]),
                 sum(dat$N_yais_end[nyear,,i,][!is.na(dat$N_yais_end[nyear,,i,])]))
    ylt = 1.3*max(rowSums(dat$N_yais_beg[1:(nyear),,i,]))
    plot(rowSums(dat$N_yais_beg[1:(nyear),,i,]),
         type = 'l',
         lwd = 2, 
         col = scales::alpha(subareaPal[i],0.2),
         main = inames[i], 
         col.main =subareaPal[i], 
         ylim = c(0,ylt),
         xlim = c(0,nyear-1),xaxt='n',
         xlab = "Model Year", 
         ylab = 'Numbers (M+F)')
    lines(rowSums(dat$N_yais_mid[1:(nyear),,i,]),
          type = 'l',
          lwd = 3,
          col = scales::alpha(subareaPal[i],0.4))
    lines(rowSums(dat$N_yais_end[1:(nyear),,i,]),
          type = 'l',
          lwd = 3,
          col = scales::alpha(subareaPal[i],0.8))
    legend("topright",col = c(scales::alpha(subareaPal[i],0.2),
                              scales::alpha(subareaPal[i],0.4),
                              scales::alpha(subareaPal[i],0.8)), 
           legend = c("beg",
                      "mid (move + 1/2 fishing + 1/2 M)",
                      "end (1/2 fishing + 1/2 M)"), lty =1, lwd = 5)
    axis(1, at = seq(1,nyear,5), labels = years[seq(1,nyear,5)])
  }
  if(!justPlots) dev.off()
  
  ## ssb_yk ----
  SSB_yk <- data.frame(dat$SSB_yk) ## 
  names(SSB_yk) <- c(paste('Stock',1:4))
  ssbyk_plot <- SSB_yk %>%
    mutate('Yr' = years[1:nrow(.)]) %>%
    reshape2::melt(id = c('Yr')) %>%
    ggplot(., aes(x = Yr, y = value, color = variable )) +
    scale_color_manual(values = rev(demPal)) +
    geom_line(lwd = 2) + 
    labs(x = 'Modeled Year',y = 'SSB', color = '') +
    scale_x_continuous(limits = c(1960,1959+nyear)) +
    ggsidekick::theme_sleek(base_size = 18) +
    theme(legend.position = c(0.8,0.8))
    # theme(legend.position = c(0.8,0.8),
    #       legend.text  = element_text(color = "white"),
    #       axis.title = element_text(color = 'white'),
    #       axis.text = element_text(color = "white"),
    #       panel.background = element_rect(fill="#01514E"),
    #       plot.background = element_rect(fill = "#01514E")) 
    
  
  if(justPlots) ssbyk_plot
  if(!justPlots) ggsave(ssbyk_plot,
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
  ## SSB_ym with compare ----
  spmat <- data.frame(subarea = c('A4',"A3","B3","B2","C2","C1"),
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
    ggplot(., aes(x = Year, y = log(value), color = REG)) +
    # ggplot(., aes(x = Year, y = assSSBMT, color = REG)) +
    geom_line(aes(y = log(value)),lwd = 1.1) +
    # geom_errorbar(aes(ymin = assSSBMT-CV*assSSBMT,ymax= assSSBMT+CV*assSSBMT, color = REG,width=0)) +
    geom_errorbar(aes(ymin = log(value-CV*value),
                      ymax = log(value+CV*value), 
                      color = REG,width=0)) +
    scale_color_manual(values = mgmtPal)+
    geom_point()+ 
    ggsidekick::theme_sleek() +
    scale_x_continuous(limits = c(1960,1959+nyear)) +
    labs(x = 'Modeled Year',y = 'log SSB (units vary)', color = 'Mgmt Region') +
    facet_wrap(~REG+variable,scales = 'free_y', ncol = 2)
  # facet_wrap(~REG,scales = 'free_y')
  ggsave(last_plot(),
         file = paste0(dumpfile,"/", Sys.Date(),'-LOG_SSB_ym.png'),
         width = 10, height = 6, unit = 'in',
         dpi = 420)
  
  merge(assSB, SSB_ym, by.x = c('Year','REG'), by.y = c('Yr','mgmt')) %>%
    mutate('Last Assessment' = assSSBMT, 'Operating Model' =omSSBMT) %>%
    select(Year, REG, CV, 'Last Assessment', 'Operating Model') %>%
    # filter(Year >1965 ) %>%
    melt(id = c('Year','REG', 'CV')) %>%
    ggplot(., aes(x = Year, y = value, color = REG)) +
    # ggplot(., aes(x = Year, y = assSSBMT, color = REG)) +
    geom_line(aes(y = value),lwd = 1.1) +
    # geom_errorbar(aes(ymin = assSSBMT-CV*assSSBMT,ymax= assSSBMT+CV*assSSBMT, color = REG,width=0)) +
    geom_errorbar(aes(ymin = (value-CV*value),
                      ymax = (value+CV*value), 
                      color = REG,width=0)) +
    scale_color_manual(values = mgmtPal) +
    geom_point()+ 
    ggsidekick::theme_sleek() +
    scale_x_continuous(limits = c(1960,1959+nyear)) +
    labs(x = 'Modeled Year',y = ' SSB (units vary)', color = 'Mgmt Region') +
    facet_wrap(~REG+variable,scales = 'free_y', ncol = 2)
  # facet_wrap(~REG,scales = 'free_y')
  ggsave(last_plot(),
         file = paste0(dumpfile,"/", Sys.Date(),'-SSB_ym.png'),
         width = 10, height = 6, unit = 'in',
         dpi = 420)
  
  ## plot bias ramp ----
  png(file =paste0(dumpfile,"/",
                       Sys.Date(),'-bias_ramp.png'),
      height = 4, width = 6, unit = 'in', res = 420)
  plot(df$parms$b_y[1:59],  type = 'l', lwd = 2, 
       xlab = 'Year',
       xaxt = 'n',
       ylab = 'B_y',
       ylim = c(0,1.2*max(df$parms$b_y,round(best[names(best) == 'b_y'],3) )))
  axis(side = 1, at = seq(0,60,5), labels= seq(1960,2020,5))

  ## blue for estimated points
  points(best[names(best) == 'b_y'][which(!is.na(mappy$b_y))], col = 'blue', pch = 19) 
  
  save(best, file = paste0(dumpfile,"/best.Rdata"))
  
  ## grey for fixed points
  points(df$parms$b_y[which(is.na(mappy$b_y))], col = 'grey44', pch = 19) 
  legend('topright',
         col = c('black','blue','grey44'),
         pch = c(NA,19,19),
         lwd = 2,
         lty = c(1,NA,NA),
         legend = c('input (start)',
                    'estimated',
                    'fixed'))
  dev.off()
  ## plot tilde_ry ----
  load( here('input','wcsummaryoutput.Rdata'))
  
  parameters <- wc$parameters
  recdev <- parameters[substring(parameters[["Label"]], 1, 12) %in% c("Main_RecrDev"), ]

  wcrd <- recdev %>% select(Value,Label) %>% mutate(Year = substr(Label,14,17)) %>%
    select(Value,Year)
  
  for(r in 1:2){
    if(r == 1){
      png(file =paste0(dumpfile,"/",
                       Sys.Date(),'-recdevs.png'),
          height = 4, width = 6, unit = 'in', res = 420)
    } else{
      png(file =paste0(dumpfile,"/",
                       Sys.Date(),'-recdevs_WC.png'),
          height = 4, width = 6, unit = 'in', res = 420)
    }
    plot(best[names(best) == 'tildeR_y'] [1:59], 
         type = 'b', 
         pch = 19,
         xlab = 'Year',
         xaxt = 'n',
         col = 'grey44',
         ylab = 'Log Recruitment Deviation',
         ylim = c(4*min(best[names(best) == 'tildeR_y'] ),1.2*max(best[names(best) == 'tildeR_y'] )))
    axis(side = 1, at = seq(0,60,5), labels= seq(1960,2020,5))
    abline(h=0,col = 'blue')
    if(r == 2){
      lines(wcrd$Value[wcrd$Year >1959], col = mgmtPal[3], pch = 19, lwd = 2)
      legend('bottomleft',
             cex = 0.8,
             legend = c('OM','WC-2019','BC-placeholder','AK-placeholder'),
             pch = c(19,NA,NA,NA),
             lty = c(NA,1,1,1),
             lwd = 2,
             col = c('grey44',rev(mgmtPal))
      )
    }
    
    dev.off()
  }
  
  # R_ym versus past assessments (?) ----
  yvec = read.csv(here('input','downloads','AssessmentDat_thru2018.csv'),fileEncoding="UTF-8-BOM") %>%
    mutate(REG = substr(Index,1,2)) %>% 
    filter(Type == 'Recruitment') %>%
    select(Year, Value, REG) %>% 
    melt(id = c('Year', 'REG')) %>% select(Year)
  assRec <- read.csv(here('input','downloads','AssessmentDat_thru2018.csv'),fileEncoding="UTF-8-BOM") %>%
    mutate(REG = substr(Index,1,2)) %>% 
    filter(Type == 'Recruitment') %>%
    select(Year, Value, REG) %>% 
    melt(id = c('Year', 'REG')) %>% 

    select(-variable) %>%
    mutate(variable = REG) %>%

    # group_by(REG) %>%
    # summarise(stdRec = value/mean(value)) %>%
    # ungroup() %>%
    
    # mutate(stdRec = value/mean(value)) %>%
    # mutate(Year = yvec$Year) %>%
    mutate(variable = REG) %>%
    select(-REG) %>%
    
    mutate(SRC = 'Assessment') %>%
    
    filter(Year > 1959) %>%
    select(Year, variable, value, SRC)
  
  R_ym <- data.frame(dat$R_ym)[1:nyear,]
  names(R_ym) <- c('WC','BC','AK')
  yvec <-   R_ym %>% data.frame() %>%
    mutate(Year = 1960:2019) %>%
    melt(id = 'Year') %>% select(Year)
  R_ym %>% data.frame() %>%
    mutate(Year = 1960:2019) %>%
    melt(id = 'Year') %>%
    # group_by(variable) %>%
    # mutate(stdRec = value/mean(value)) %>%
    # summarise(stdRec = value/mean(value)) %>%
    # select(-value) %>%
    # mutate(Year = yvec$Year) %>%
    mutate(SRC = 'Operating Model') %>%
    # mutate(variable = paste("OM_",variable)) %>%
  rbind(., assRec) %>%
    ggplot(., aes(x = Year, y = value, col = variable, linetype = SRC)) +
    geom_line(lwd = 1) +
    theme_sleek() +
    scale_linetype_manual(values = c('dotted','solid'))+
    scale_color_manual(values = rev(mgmtPal)) +
    labs(x = 'Year', y = 'Raw Recruitment',
         subtitle = 'Assessment values taken from Fenske et. al. 2018', 
         color = '',linetype = '') +
    facet_wrap(~variable+SRC, scale = 'free_y', ncol = 2)
    
  ggsave(last_plot(),
         file = paste0(dumpfile,"/", Sys.Date(),'-R_ym.png'),
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
  
  ## plot acomps by fleet----
  ## recall that the om inputs are set up in same order as fleetnames_acomp
  acomp_yafs_obs <- df$acomp_yafs_obs
  dimnames(acomp_yafs_obs) <- list(c(1960:(df$nyear+1959)),
                              c(0:70),
                              as.character(paste(unlist(fltnames_acomp))),
                              c('Fem','Mal'))
  plotacomp <- list()
  for(flt in 1:nfleets_acomp){
    obsacomp <- predacomp <- NULL
    for(s in 1:2){
      ## get position in predictions
      predpos <- ifelse(phi_ff_acomp[flt, "commacomp_pos"] != -1,
                       phi_ff_acomp[flt, "commacomp_pos"] + 1,
                       phi_ff_acomp[flt, "survacomp_pos"] + 1)
      obsacomp <- rbind(obsacomp,melt(acomp_yafs_obs[,,flt,s]) %>% 
        select("year" = Var1, "age" = Var2, "freq" = value) %>%
        mutate(fleet = fltnames_acomp[flt], 
               sex = c('Fem','Mal')[s],
               src = 'OBS'))
      if(df$acomp_flt_type[flt] == 0){
        predacomp <- rbind(predacomp,melt(comm_acomp_yafs_pred[,,predpos,s]) %>% 
          select("year" = Var1, "age" = Var2, "freq" = value) %>%
          mutate(fleet = fltnames_acomp[flt], 
                 sex = c('Fem','Mal')[s],
                 src = 'PRED'))
      } else{
        predacomp <- rbind(predacomp,melt(surv_acomp_yafs_pred[,,predpos,s]) %>% 
          select("year" = Var1, "age" = Var2, "freq" = value) %>%
          mutate(fleet = fltnames_acomp[flt], 
                 sex = c('Fem','Mal')[s],
                 src = 'PRED'))
      } ## end else
    } ## end sexes
      plotacomp[[flt]] <- rbind(obsacomp,predacomp) %>% filter(year > 2015 & freq >0)
  } ## end fleets
  
  
  ggplot(plotacomp[[1]], aes(x = age, y = freq, fill = src, color = sex)) +
    theme_sleek() +
    geom_ribbon(data = subset(plotacomp[[1]], sex == 'Fem'),
                lwd = 1.05,   aes(x=age,ymax=freq),ymin=0,alpha=0.3) +
    geom_ribbon(data = subset(plotacomp[[1]], sex == 'Mal'),  lwd = 1.05,
                aes(x=age,ymax=-freq),ymin=0,alpha=0.3) +
    scale_fill_manual(values = 'grey22') +
    scale_color_manual(values = c('red','blue')) +
    facet_wrap(~year)



  
  
  ## catch pred by fleet ----
  catch_yf_pred_totalt <- data.frame(dat$catch_yf_pred_total)[1:nyear,]
  names(catch_yf_pred_totalt) <- unlist(df$fltnames_fish)
  catch_yf_pred_totalt <- catch_yf_pred_totalt %>%
    mutate(Year = years) %>%
    group_by(Year) %>%
    melt(id = 'Year') %>%
    mutate(Type = 'PRED') %>%
    mutate(REG = substr(variable,0,2)) %>% 
    filter(value != 0  | variable == 'AK_FIX')
  
  catch_yf_obst <- df$catch_yf_obs[1:nyear,] %>%  data.frame() %>%select(-Year) 
  
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
    scale_x_continuous(limits = c(1960,1959+nyear)) +
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
    scale_x_continuous(limits = c(1960,1959+nyear)) +
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
  
  survey_yf_predt <- data.frame(dat$surv_yf_pred)[1:nyear,]
  names(survey_yf_predt) <- unlist(df$fltnames_surv)
  # names(survey_yf_predt)[3:4] <- c('BC_OFFStd','BC_StRs')
  survey_yf_predt <- survey_yf_predt %>%
    mutate(Year = years) %>%
    melt(id = 'Year') %>%
    mutate(Type = 'PRED') %>%
    mutate(REG = substr(variable,0,2)) %>% 
    filter(value > 0)
  
  survey_yf_errt <- data.frame(df$surv_yf_err)[1:nyear,] %>%    
    mutate(Year = years) %>%
    melt(id = 'Year') %>%
    filter(!is.na(value)) 
  
  
  survey_yf_obst <- data.frame( df$surv_yf_obs)[1:nyear,]   %>%
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
  
  melt(df$srv_blks) %>% 
    filter(value != 59) %>% 
    mutate(variable = Var2, blk = value+1960) %>%
    select(variable, blk) %>% 
    merge(.,
          survey_yf_obst,
          by = 'variable', all.y = TRUE) %>%
    ggplot(data = ., 
           aes(x = Year, y = value, color = variable)) +
    geom_vline(aes( xintercept = blk), col = 'grey75', 
               lwd = 1, linetype = 'dashed') +
    
    geom_line(data = survey_yf_predt, lwd = 0.75) +
    scale_color_manual(values = survfltPal) +
    geom_point(pch = 1, fill = NA, col = 'black') +
    geom_errorbar(aes(ymin = lci, ymax = uci), col = 'black',width=0) +
    scale_x_continuous(limits = c(1980,1960+nyear), 
                       breaks = seq(1980,1960+nyear,10),
                       labels = seq(1980,1960+nyear,10)) +
    theme_sleek() + 
    theme(legend.position = 'none')+
    labs(y = 'survey', color = 'Fishing Fleet')+
    facet_wrap(~variable, scales = "free_y")
  
  ggsave(last_plot(),
         file =paste0(dumpfile,"/", Sys.Date(),'-survey_fits.png'),
         width = 10, height = 6, unit = 'in',
         dpi = 420)
  
  ## plot length_yais thru time ----
  
  templyais <-  templyais_beg <-templyais_mid <-templyais_end <- NULL
  for(i in 1:6){
    for(y in 1:nyear){
      templyais_beg <- bind_rows(templyais_beg,
                                 dat$Length_yais_beg[y,,i,] %>%
                                   melt() %>% 
                                   plyr::rename(c('Var1' = 'age','Var2' = 'sex', 'value' = 'length')) %>% 
                                   mutate(year =  y, 
                                          subarea = paste0("Subarea ",i),
                                          sex = ifelse(sex == 1,'Fem','Mal')))%>%
        mutate(PD = 'BEG')
      templyais_mid <- bind_rows(templyais_mid,
                                 dat$Length_yais_mid[y,,i,] %>%
                                   melt() %>% 
                                   plyr::rename(c('Var1' = 'age','Var2' = 'sex', 'value' = 'length')) %>% 
                                   mutate(year =  y, 
                                          subarea = paste0("Subarea ",i),
                                          sex = ifelse(sex == 1,'Fem','Mal')))%>%
        mutate(PD = 'MID')
      templyais_end <- bind_rows(templyais_end,
                                 dat$Length_yais_end[y,,i,] %>%
                                   melt() %>% 
                                   plyr::rename(c('Var1' = 'age','Var2' = 'sex', 'value' = 'length')) %>% 
                                   mutate(year =  y, 
                                          subarea = paste0("Subarea ",i),
                                          sex = ifelse(sex == 1,'Fem','Mal'))) %>%
        mutate(PD = 'END')
    }
  }
  templyais <- rbind(templyais_beg,templyais_mid,templyais_end)
  templyais$PD <- factor(templyais$PD, 
                         levels = c('BEG','MID','END'))
  ggplot(filter(templyais, sex == 'Fem') , 
         aes(x = age, y = length, group = year)) +
    geom_line(aes(color = year), alpha = 0.2, lwd = 0.5) +
    scale_color_viridis_c() +
    ggsidekick::theme_sleek() +
    labs(x = 'Age (years)', y = 'Length (cm)', 
         title = 'Mean LAA in Subarea, Females') +
    facet_grid(subarea ~ PD )
  
  ggsave(last_plot(),
         file =paste0(dumpfile,"/",
                      Sys.Date(),'-meanLAA_fem.png'),
         width = 8, height = 6, units = 'in', dpi = 520)
  
  ggplot(filter(templyais, sex == 'Mal') , 
         aes(x = age, y = length, group = year)) +
    geom_line(aes(color = year), alpha = 0.2, lwd = 0.5) +
    scale_color_viridis_c() +
    ggsidekick::theme_sleek() +
    labs(x = 'Age (years)', y = 'Length (cm)', 
         title = 'Mean LAA in Subarea, Males') +
    facet_grid(subarea ~ PD )
  
  ggsave(last_plot(),
         file =paste0(dumpfile,"/",
                      Sys.Date(),'-meanLAA_mal.png'),
         width = 8, height = 6, units = 'in', dpi = 520)
  
  ## plot dist of LAA using length_alyis ----
  # pA1 <- list() ## for area 1
  # for(y in 1:5){ #45:dim(LengthAge_alyi_beg)[3]){ ## loop years
  #   a1 <- 
  #     LengthAge_alyis_end[,,y,1,1] %>%
  #     melt() %>%
  #     group_by(Var2) %>% 
  #     mutate(sumP = sum(value), ## total prob within A-L bin
  #            pbin = value/sumP) ## relative prob within A-L bin
  #   
  #   ggplot(a1,aes( x = Var2, y = pbin, color = factor(Var1))) +
  #     # geom_histogram(stat ='identity', position = 'stack') +
  #     geom_line() + 
  #     geom_area(aes(fill=factor(Var1))) +
  #     # geom_density( )+
  #     labs(x = 'len', y = 'prob(Bin)', fill = 'age', 
  #          title = paste("year ",y), 
  #          subtitle= 'subarea 1') +
  #     theme_sleek()
  #   
  #   ggplot(a1,aes( x = Var1, y = pbin, color = factor(Var2))) +
  #     # geom_histogram(stat ='identity', position = 'stack') +
  #     geom_line() + 
  #     # geom_area(aes(fill=factor(Var2))) +
  #     # geom_density( )+
  #     labs(x = 'age', y = 'prob(Bin)', fill = 'len', 
  #          title = paste("year ",y), 
  #          subtitle= 'subarea 1') +
  #     theme_sleek()
  #   
  #   pA1[[y]] <- ggplot(a1,aes(x = Var1, y = Var2, fill = pbin)) +
  #     geom_tile() +
  #     labs(x = 'age', y = 'len', 
  #          title = paste("year ",y), 
  #          subtitle= 'subarea 1') +
  #     theme_sleek()
  #   rm(a1)
  #   # p[[i]] <- qplot(1:10,10:1,main=i)
  # }
  # # png(here("figs","LAA_Dist_A1.png"), 
  # #     height = 10, width = 10, unit = 'in', res = 420)
  # do.call(gridExtra::grid.arrange,pA1)
  # dev.off()
  
  # plot FISH selex ----
  ## bring estimates out and rearrange
  nfishmod = df$nfleets_fish ## bc we got selex for acomp flts too
  Nas <- which(is.na(mappy[[grep("log_fsh_slx_pars", names(mappy))]]))
  nfixedfleets <- length(Nas)/4
  nfishmod <- nfishmod - nfixedfleets
  
  ## if at least some were est:
  if( length(obj$par[grep('log_fsh_slx_pars',names(obj$par))]) != 0 ){
    
    ## use map to match input pars which were actually used
    inputSel <- array(exp(obj$par[grep('log_fsh_slx_pars',names(obj$par))]),
                      dim = dim(df$parms$log_fsh_slx_pars),
                      dimnames = dimnames(df$parms$log_fsh_slx_pars))

    map_fshslx <- array(as.numeric(mappy$log_fsh_slx_pars), 
                        dim = c(df$nfleets_fish,2,max(df$fsh_blks_size),2),
                        dimnames = dimnames(df$parms$log_fsh_slx_pars))
    ## replace non-fixed values with starting pars
    map_fshslx[!is.na(map_fshslx)] <-  array(as.numeric(exp(best[grep('log_fsh_slx_pars',names(best))])))
    selP <- map_fshslx
    fsh_sel_afsb <- array(NA, dim =  c(df$nage,
                                       df$nfleets_fish,
                                       2,
                                       max(df$fsh_blks_size)),
                          dimnames = list(c(df$age),
                                          c(dimnames(df$parms$log_fsh_slx_pars)[[1]]),
                                          c('Fem','Mal'),
                                          c(dimnames(df$parms$log_fsh_slx_pars)[[3]])))
    
    for(fish_flt in 1:df$nfleets_fish){
      for(blk in 1:df$fsh_blks_size[fish_flt]){
        for(s in 1:2){
          if(!is.na(map_fshslx[fish_flt,1,blk,s])){
            fsh_sel_afsb[,fish_flt,s,blk] <- getSelec2(sex = s,
                                                       selP = selP,
                                                       flt_idx = fish_flt,
                                                       selType = df$selType_fish[fish_flt], 
                                                       selShape = df$selShape_fish[fish_flt],
                                                       fltType = 'fish')
          } else if(is.na(map_fshslx[fish_flt,1,blk,s])){
            fsh_sel_afsb[,fish_flt,s,blk] <- getSelec2(sex = s,
                                                       selP = inputSel,
                                                       flt_idx = fish_flt,
                                                       selType = df$selType_fish[fish_flt], 
                                                       selShape = df$selShape_fish[fish_flt],
                                                       fltType = 'fish')
          }
        }  ## end sex
      } ## end blk
    }## end fish fleet
    
    for(blk in 1:max(df$fsh_blks_size)){
      png(paste0(dumpfile,'/fishey_selex_blk',blk,".png"),
      height = 8, width = 6, unit = 'in', res = 420)
      par(mfrow = c(4,2) )
      for(flt in 1:df$nfleets_fish){
        ## if not in this block, skip
        if(is.na(fsh_sel_afsb[1,flt,1,blk])) next()
        ## if fixed overwrite colors
        if(is.na(map_fshslx[flt,1,blk,1])){
          sexPal_temp = c('grey22','grey66')
        } else{
          sexPal_temp = sexPal
        }
        for(s in 1:2){
          tmp <- fsh_sel_afsb[,flt,s,blk]
          if(s == 1) plot(tmp, 
                          col = sexPal_temp[1], 
                          type = 'l', 
                          lwd = 2, 
                          xlab = ifelse(df$selType_fish[flt] == 0,
                                        'Age','Length'), 
                          ylab = 'Selectivity',
                          lty = 1,
                          ylim = c(0,1), 
                          main = paste0(dimnames(df$parms$log_fsh_slx_pars)[[1]][flt],
                                        " ", df$fsh_blks[blk,flt]+1960),
                          xlim = c(0,75),
                          col.main  = c(rep(mgmtPal[1],4), rep(mgmtPal[2],3),rep(mgmtPal[3],2))[flt])
          text(x=60, y = 0.2, cex = 1.1,label = 
                 paste0(ifelse(df$selType_fish[flt] == 0,
                               'A','L'),"50=", 
                        round(selP[flt,1,blk,s],1),"\n",
                        paste0(ifelse(df$selType_fish[flt] == 0,
                                      'A','L'),"95=", 
                               round(selP[flt,2,blk,s],2))))
              box(which = 'plot', lty = 'solid', 
              col = c(rep(mgmtPal[1],4), rep(mgmtPal[2],3),rep(mgmtPal[3],2))[flt], 
              lwd = 2)
          if(s == 2) lines(tmp, col = sexPal_temp[2], type = 'l', lty = 2, lwd = 2)
        }
      }
      dev.off()
    } ## end blks
    
  } else{
    inputSel <- selP <- exp(df$parms$log_fsh_slx_pars)
    fsh_sel_afs <- array(NA, dim =  c(df$nage,df$nfleets_fish,2),
                         dimnames = list(c(df$age),
                                         c(dimnames(df$parms$log_fsh_slx_pars)[[1]]),
                                         c('Fem','Mal')))
    for(s in 1:2){
      for(fish_flt in 1:df$nfleets_fish){
        fsh_sel_afs[,fish_flt,s] <- getSelec2(sex = s,
                                              selP = selP,
                                              flt_idx = fish_flt,
                                              selType = df$selType_fish[fish_flt], 
                                              selShape = df$selShape_fish[fish_flt],
                                              fltType = 'fish')
        
      } ## end fish fleet
    } ## end sex
    
    png(paste0(dumpfile,'/fishery_selex.png'),
        height = 8, width = 6, unit = 'in', res = 420)
    par(mfrow = c(4,2) )
    sexPal_temp = c('grey22','grey66')
    for(flt in 1:df$nfleets_fish){
      for(s in 1:2){
        tmp <- fsh_sel_afs[,flt,s]
        if(s == 1) plot(tmp, 
                        col = sexPal_temp[1], 
                        type = 'l', 
                        lwd = 2, 
                        xlab = ifelse(df$selType_fish[flt] == 0,
                                      'Age','Length'), 
                        ylab = 'Selectivity',
                        lty = 1,
                        ylim = c(0,1), 
                        main = dimnames(df$parms$log_fsh_slx_pars)[[1]][flt], xlim = c(0,75),
                        col.main  = c(rep(mgmtPal[1],4), rep(mgmtPal[2],3),rep(mgmtPal[3],2))[flt])
        box(which = 'plot', lty = 'solid', 
            col = c(rep(mgmtPal[1],4), rep(mgmtPal[2],3),rep(mgmtPal[3],2))[flt], 
            lwd = 2)
        if(s == 2) lines(tmp, col = sexPal_temp[2], type = 'l', lty = 2, lwd = 2)
      }
    }
    dev.off()
    
  } ## end if all fixed
  
  # plot SURV selex ----
  ## bring estimates out and rearrange
  nsurvmod = length(df$selType_surv) ## bc we got selex for acomp flts too
  Nas <- which(is.na(mappy[[grep("log_srv_slx_pars", names(mappy))]]))
  nfixedfleets <- length(Nas)/4
  nsurvmod <- nsurvmod - nfixedfleets
  
  ## if at least some were est:
  if( length(obj$par[grep('log_srv_slx_pars',names(obj$par))]) != 0 ){
    
    ## use map to match input pars which were actually used
    inputSel <- array(exp(df$parms$log_srv_slx_pars),
                      dim = dim(df$parms$log_srv_slx_pars),
                      dimnames = dimnames(df$parms$log_srv_slx_pars))

    map_srvslx <- array(as.numeric(mappy$log_srv_slx_pars), 
                        dim = c(length(df$selType_surv),2,max(df$srv_blks_size),2),
                        dimnames = dimnames(df$parms$log_srv_slx_pars))
    ## replace non-fixed values with starting pars
    map_srvslx[!is.na(map_srvslx)] <-  array(as.numeric(exp(best[grep('log_srv_slx_pars',names(best))])))
    selP <- map_srvslx
    
    srv_sel_afsb <- array(NA, dim =  c(df$nage,
                                   length(df$selType_surv),
                                      2,
                                      max(df$srv_blks_size)),
                         dimnames = list(c(df$age),
                                         c(dimnames(df$parms$log_srv_slx_pars)[[1]]),
                                         c('Fem','Mal'),
                                         c(dimnames(df$parms$log_srv_slx_pars)[[3]])))

    for(surv_flt in 1:length(df$selType_surv)){
      for(blk in 1:df$srv_blks_size[surv_flt]){
        for(s in 1:2){
          if(!is.na(map_srvslx[surv_flt,1,blk,s])){
            srv_sel_afsb[,surv_flt,s,blk] <- getSelec2(sex = s,
                                                      selP = selP,
                                                      blk = blk,
                                                      flt_idx = surv_flt,
                                                      selType = df$selType_surv[surv_flt], 
                                                      selShape = df$selShape_surv[surv_flt],
                                                      fltType = 'surv')
          } else if(is.na(map_srvslx[surv_flt,1,blk,s])){
            srv_sel_afsb[,surv_flt,s,blk] <- getSelec2(sex = s,
                                                  selP = inputSel,
                                                  flt_idx = surv_flt,
                                                  blk= blk,
                                                  selType = df$selType_surv[surv_flt], 
                                                  selShape = df$selShape_surv[surv_flt],
                                                  fltType = 'surv')
          }
        }  ## end sex
      } ## end blk
    }## end surv fleet
    
    for(blk in 1:max(df$srv_blks_size)){
      png(paste0(dumpfile,'/survey_selex_blk',blk,".png"),
          height = 8, width = 6, unit = 'in', res = 420)
      par(mfrow = c(4,2) )
      for(flt in 1:length(df$selType_surv)){
        ## if not in this block, skip
        if(is.na(srv_sel_afsb[1,flt,1,blk])) next()
        ## if fixed overwrite colors
        if(is.na(map_srvslx[flt,1,blk,1])){
          sexPal_temp = c('grey22','grey66')
        } else{
          sexPal_temp = sexPal
        }
        for(s in 1:2){
          tmp <- srv_sel_afsb[,flt,s,blk]
          if(s == 1) plot(tmp, 
                          col = sexPal_temp[1], 
                          type = 'l', 
                          lwd = 2, 
                          xlab = ifelse(df$selType_surv[flt] == 0,
                                        'Age','Length'), 
                          ylab = 'Selectivity',
                          lty = 1,
                          ylim = c(0,1), 
                          main = paste0(dimnames(df$parms$log_srv_slx_pars)[[1]][flt],
                                       " ", df$srv_blks[blk,flt]+1960),
                          xlim = c(0,75),
                          col.main  = c(survfltPal,rep('black',1))[flt])
                            # c(rep(mgmtPal[1],4), rep(mgmtPal[2],3),rep(mgmtPal[3],2))[flt])
          text(x=60, y = 0.2, cex = 1.1,label = 
                 paste0(ifelse(df$selType_surv[flt] == 0,
                               'A','L'),"50=", 
                        round(selP[flt,1,blk,s],1),"\n",
                        paste0(ifelse(df$selType_surv[flt] == 0,
                                      'A','L'),"95=", 
                               round(selP[flt,2,blk,s],2))))
          box(which = 'plot', lty = 'solid', 
              col = survfltPal[flt],
              # col = c(rep(mgmtPal[1],4), rep(mgmtPal[2],3),rep(mgmtPal[3],2))[flt], 
              lwd = 2)
          if(s == 2) lines(tmp, col = sexPal_temp[2],
                           type = 'l', lty = 2, lwd = 2)
          # seeddim = nestflts
          # lines(bounds$lower type = 'v', col = sexPal[2])
        }
      }
      dev.off()
    } ## end blks

  } else{ ## all were fixed
    inputSel <- selP <- array(exp(df$parms$log_srv_slx_pars), dim = c(df$nfleets_surv+df$nfleets_acomp-(df$nsurvflts_acomp+df$nfishflts_acomp),2,
                                                                      max(df$srv_blks_size),2),
                              dimnames = dimnames(df$parms$log_srv_slx_pars))
    # array(exp(best[grep('log_srv_slx_pars',names(best))]),
    #                         dim = c(df$nfleets_surv+df$nfleets_acomp-4,2,max(df$srv_blks_size),2),
    #                         dimnames = dimnames(df$parms$log_srv_slx_pars))
    
    srv_sel_afs <- array(NA, dim =  c(df$nage,length(df$selType_surv),2),
                         dimnames = list(c(df$age),
                                         c(dimnames(df$parms$log_srv_slx_pars)[[1]]),
                                         c('Fem','Mal')))
    for(s in 1:2){
      for(surv_flt in 1:length(df$selType_surv)){
        srv_sel_afs[,surv_flt,s] <- getSelec2(sex = s,
                                              selP = selP,
                                              flt_idx = surv_flt,
                                              selType = df$selType_surv[surv_flt], 
                                              selShape = df$selShape_surv[surv_flt],
                                              fltType = 'surv')
        
      } ## end surv fleet
    } ## end sex
    
    png(paste0(dumpfile,'/survey_selex.png'),
        height = 8, width = 6, unit = 'in', res = 420)
    par(mfrow = c(4,2) )
    sexPal_temp = c('grey22','grey66')
    for(flt in 1:length(df$selType_surv)){
      for(s in 1:2){
        tmp <- srv_sel_afs[,flt,s]
        if(s == 1) plot(tmp, 
                        col = sexPal_temp[1], 
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
        if(s == 2) lines(tmp, col = sexPal_temp[2], type = 'l', lty = 2, lwd = 2)
      }
    }
    dev.off()
    
  } ## end if all fixed
  

  
  
  
  
  
  
} ## end func
