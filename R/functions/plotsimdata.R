

## this can take either a list of replicates or single value
plotSimData<- function(simdat, new_strata =6){
  savedir <- here('output','om_simulations',paste0(Sys.time(),"_n=", ifelse(is.null(dim(simdat)[2]),1,dim(simdat)[2])))
  dir.create(savedir)
  save(simdat, file = paste0(savedir, 'sim.rdata'))
  if(class(simdat) == 'list'){
    ## not replicates, simple plot
    ctch <- simdat$catch_yf_obs
    dimnames(ctch) <- dimnames(df$catch_yf_obs)
    ggplot(melt(data.frame(ctch), id = c('Year')), aes(x = Year,y = value, color = variable)) +
      geom_rect(data = NULL, aes(xmin =  df$years[df$yRun], xmax = Inf,
                                 ymin = -0.5, ymax = Inf), fill = 'grey 88', color = NA, 
                alpha = 0.5) +
      geom_line(lwd = 1) +
      theme_sleek() +
      scale_color_manual(values = fishfltPal)
    ggsave(last_plot(), file = paste0(savedir, 'catch.png'), height = 4, width = 6, dpi = 520 )
    
    surv <- simdat$surv_yf_obs
    dimnames(surv) <- dimnames(df$surv_yf_obs)
    ggplot(melt(data.frame(surv, 'Year' =df$catch_yf_obs[,1]), 
                id = c('Year')), aes(x = Year,y = value, color = variable)) +
      geom_rect(data = NULL, aes(xmin =  df$years[df$yRun], xmax =  Inf,
                                 ymin = -0.5, ymax = Inf), fill = 'grey 88', color = NA, 
                alpha = 0.5) +
      geom_line(lwd = 1) +
      scale_x_continuous(limits = c(1980, max(df$years)))+
      theme_sleek() +
      scale_color_manual(values = survfltPal)
    
    
    
  } else{
    nsim = dim(simdat)[2]
    lapply(sim,mutate(src = nsim))
    
    ## only show squid for proj years
    
  }
  
  
}
