## a function to generate simulated data (ideally from a conditioned OM, but pars are inputs)
## and create plots & return datasets usable in the format aWCepted by shire.
## likely will need to do this for phi objects as well.

## passing it a seed will return a unique dataset; can embed in loop if you want many replicates
## at this stratification
## options for new_strata are 6 (return at OM dims), 4 (stock dims), 3 (mgmt dims) or 1 (panmictic)
## this only changes the colors for plotting, not the actual input data, which is fleet-specific anyway
format_phi <- function(df = NULL, 
                       no_strata = 3,
                       survey, catch,acomp_flt_type){
nstocks = 4
    ## get idx of all dfs with phi
    spmat6 <- data.frame(subarea = c('A4',"A3","B3","B2","C2","C1"),
                        stock = c("R4","R3","R3","R2","R2","R1"),
                        mgmt = c("AK","AK", rep("BC",2), rep("WC",2)))
    
    spmat4 <- data.frame(subarea = c('S4',"S3","S3","S2","S2","S1"),
                         stock = c("R4","R3","R3","R2","R2","R1"),
                         mgmt = c("AK","AK", rep("BC",2), rep("WC",2)))
    
    spmat3 <- data.frame(subarea = c('M3',"M3","M2","M2","M1","M1"),
                         stock = c("R4","R3","R3","R2","R2","R1"),
                         mgmt = c("AK","AK", rep("BC",2), rep("WC",2)))
    
    if(no_strata == 6){
      spmat <- spmat6
    } else  if(no_strata == 4){
      spmat <- spmat4
    }else{
      spmat <- spmat3
    }

    ## phi_survey
    phi_if_surv <- matrix(0, nrow = nfleets_surv, ncol = no_strata)
    rownames(phi_if_surv) <- names(survey)
    colnames(phi_if_surv) <- unique(rev(spmat$subarea))
    for(i in 1:nrow(phi_if_surv)){
      reg = substr(rownames(phi_if_surv)[i],1,2)
      SA <- unique(spmat[which(spmat$mgmt == reg),"subarea"])
      phi_if_surv[i,which(colnames(phi_if_surv)  %in% SA)] <- 1
    }
    
    phi_if_acomp <- matrix(0, nrow = nfleets_acomp, ncol = no_strata)
    rownames(phi_if_acomp) <- fltnames_acomp
    colnames(phi_if_acomp) <-unique(rev(spmat$subarea))
    for(i in 1:nrow(phi_if_acomp)){
      reg = substr(rownames(phi_if_acomp)[i],1,2)
      SA <- unique(spmat[which(spmat$mgmt == reg),"subarea"])
      phi_if_acomp[i,which(colnames(phi_if_acomp)  %in% SA)] <- 1
    }
    

    
    ## indicates the position of acomp fleet (shouldn't change with subareas)
    phi_ff_acomp <- matrix(-1, nrow = nfleets_acomp, ncol = 5) 
    rownames(phi_ff_acomp) <- fltnames_acomp
    colnames(phi_ff_acomp) <- c('fsh_slx_pos','srv_slx_pos',"nsamp_pos","commacomp_pos","survacomp_pos")
    
    phi_ff_acomp[which(rownames(phi_ff_acomp) %in% paste(fltnames_fish)),1] <- 
      which(grepl(paste(rownames(phi_ff_acomp), collapse = "|"), paste(fltnames_fish)))-1
    
    # if(tolower(x) == 'n'){
    #   phi_ff_acomp[Wnsurvflts_acomp,'srv_slx_pos'] <- 
    #     which(fltnames_surv %in% fltnames_acomp )-1## Pos in survey
    #   
    #   phi_ff_acomp[acomp_flt_type == 1,"survacomp_pos"] <- 
    #     which(fltnames_surv %in% fltnames_acomp[acomp_flt_type == 1])-1
    # } else{
      ## if a standalone survey, fill with last idx
      phi_ff_acomp[acomp_flt_type == 1,"srv_slx_pos"] <-  
        (nfleets_surv +1):(nfleets_surv + sum(acomp_flt_type==1))-1
      
      phi_ff_acomp[acomp_flt_type == 1,"survacomp_pos"] <-  
        (nfleets_surv +1):(nfleets_surv + sum(acomp_flt_type==1))-1
    # }
    
    phi_ff_acomp[acomp_flt_type == 0,"commacomp_pos"] <- 
      which(fltnames_fish %in% fltnames_acomp[acomp_flt_type == 0])-1
    
    phi_ff_acomp[,"nsamp_pos"] <- 0:(nrow(phi_ff_acomp)-1)
    
    
    ## phi_fish
    phi_if_fish <- matrix(0, nrow = nfleets_fish, ncol = no_strata) ## placeholder for fishing fleets
    rownames(phi_if_fish) <- names(catch)[2:ncol(catch)]
    colnames(phi_if_fish) <-  unique(rev(spmat$subarea))
    for(i in 1:nrow(phi_if_fish)){
      reg = substr(rownames(phi_if_fish)[i],1,2)
      SA <- unique(spmat[which(spmat$mgmt == reg),"subarea"])
      phi_if_fish[i,which(colnames(phi_if_fish)  %in% SA)] <- 1
    }

    ## tau_ki
    tau_ki <-  matrix(0, ncol = no_strata, nrow = nstocks) ## nesting of subareas within stocks, for recruitment purposes
    rownames(tau_ki) <- rev(unique(spmat$stock))
    colnames(tau_ki) <- unique(rev(spmat$subarea))
    if(no_strata == 4){
      diag(tau_ki) <- 1 ## all self-seeding
    }
    if(no_strata == 3){
      tau_ki['R1','M1'] <- 1## all s 36 in a1
      tau_ki['R2','M2'] <-  0.75; tau_ki['R3','M2'] <- 0.25 ## most of R2 from CC
      tau_ki['R2','M2'] <-  0.75; tau_ki['R3','M2'] <- 0.25 ## most of R2 from CC
      tau_ki['R3','M3'] <- 0.25 ; tau_ki['R4','M3'] <- 0.75 ## most of R4 from AK
    }
    if(no_strata == 6){
      tau_ki[1,1] <-   tau_ki[4,6]  <- 1 ## 100% of recruitment in stock
      tau_ki[2,2:3] <-  tau_ki[3,5:4] <-  c(0.75,0.25) ## A2 and C1 are larger
    }
    
    
    ## phi_im
    phi_im <- matrix(0, ncol = 3, nrow = no_strata)
    colnames(phi_im) <- rev(unique(spmat$mgmt))
    rownames(phi_im) <- unique(rev(spmat$subarea))
    if (no_strata == 6) {
      phi_im[1:2, 1] <- phi_im[3:4, 2] <- phi_im[5:6, 3] <- 1
    } else if (no_strata == 3) {
      diag(phi_im) <- 1
    } else{
      phi_im['S1','WC'] <- 1; phi_im['S2','WC'] <- 0.25## all s 36 in a1
      phi_im['S2','BC'] <-  0.75; phi_im['S3','BC'] <- 0.25 ## most of S2 fSom CC
      phi_im['S3','BC'] <- 0.25 ## most of S2 fSom CC
      phi_im['S3','AK'] <- 0.75 ; phi_im['S4','AK'] <- 1 ## most of S4 fSom AK## use the same upscaling as SecSuitment to allocate eveSything else ???
      ## might break because cpp expecting IMATRIX...
    }
    ## phi_ki
    phi_ki <-  matrix(0, ncol = no_strata, nrow = nstocks) ## nesting of subareas within stocks, for recruitment purposes
    rownames(phi_ki) <- rev(unique(spmat$stock))
    colnames(phi_ki) <- unique(rev(spmat$subarea))
    
    for(i in 1:nrow(phi_ki)){
      st <- rownames(phi_ki)[i]
      SA <- unique(spmat[which(spmat$stock == st),"subarea"])
      phi_ki[i,which(colnames(phi_ki)  %in% SA)] <- 1
    }
   
    phi_ik2 <- matrix(apply(phi_ki,2, function(x)which(x == 1))-1) ## a vector for par subsetting, the columns are subareas
    
    ## phi_ij [eq 6]
    phi_ij <-  matrix(1, ncol = no_strata, nrow = no_strata) ## 0 indicates subareas comprise THE SAME stock
    rownames(phi_ij) = colnames(phi_ij) = unique(rev(spmat$subarea))
    if(no_strata == 3) phi_ij[upper.tri(phi_ij)] <- phi_ij[lower.tri(phi_ij)] <- 0  ## not convinced what to do for n = 3
    if(no_strata == 4) diag(phi_ij) <- 0 
    if(no_strata == 6) diag(phi_ij) <- phi_ij[2:3,2:3] <- phi_ij[4:5,4:5] <- 0

    ## phi_fm
    phi_fm <- matrix(0, nrow = nfleets_fish, ncol = 3)
    rownames(phi_fm) = names(catch)[2:ncol(catch)]
    colnames(phi_fm) = rev(unique(spmat$mgmt))
    for(i in 1:nrow(phi_fm)){
      reg = substr(rownames(phi_fm)[i],1,2)
      if(reg == 'AK') {
        phi_fm[i,3] <- 1
      } else  if(reg == 'BC'){
        phi_fm[i,2] <- 1
      }else{
        phi_fm[i,1] <- 1
      }
    }
    ## same as above but for comps (mix of fisheries & surveys)
    phi_fm_acomp <- matrix(0, nrow = nfleets_acomp, ncol = 3)
    rownames(phi_fm_acomp) = fltnames_acomp
    colnames(phi_fm_acomp) = rev(unique(spmat$mgmt))
    for(i in 1:nrow(phi_fm_acomp)){
      reg = substr(rownames(phi_fm_acomp)[i],1,2)
      if(reg == 'AK') {
        phi_fm_acomp[i,3] <- 1
      } else  if(reg == 'BC'){
        phi_fm_acomp[i,2] <- 1
      }else{
        phi_fm_acomp[i,1] <- 1
      }
    }
    # phi_fm_acomp[1:3,3] <- phi_fm_acomp[4:6,2]  <- phi_fm_acomp[7:8,1]  <- 1
    phi_fm_acomp2 <- matrix(apply(phi_fm_acomp,1, function(x)which(x == 1))-1) ## a vector for par subsetting, the columns are survey fleets

    phi_lcomp_fm <- matrix(0, nrow = nfleets_lcomp, ncol = 3)
    rownames(phi_lcomp_fm) = fltnames_lcomp
    colnames(phi_lcomp_fm) = rev(unique(spmat$mgmt))
    # phi_lcomp_fm[1:6,3] <- phi_lcomp_fm[7:9,2]  <- phi_lcomp_fm[10,1]  <- 1
    for(i in 1:nrow(phi_lcomp_fm)){
      reg = substr(rownames(phi_lcomp_fm)[i],1,2)
      if(reg == 'AK') {
        phi_lcomp_fm[i,3] <- 1
      } else  if(reg == 'BC'){
        phi_lcomp_fm[i,2] <- 1
      }else{
        phi_lcomp_fm[i,1] <- 1
      }
    }
   
    return(
      list(
        "phi_if_surv" = phi_if_surv,
        "phi_if_fish"=phi_if_fish,
        'phi_if_acomp'=phi_if_acomp,
        'phi_ff_acomp'=phi_ff_acomp,
        'phi_ki'=phi_ki,
        'phi_im'=phi_im,
        'phi_ik2'=phi_ik2,
        'phi_ij'=phi_ij,
        'phi_fm'=phi_fm,
        'tau_ki'=tau_ki,
        'phi_fm_acomp'=phi_fm_acomp,
        'phi_fm_acomp2'=phi_fm_acomp2,
        'phi_lcomp_fm'=phi_lcomp_fm
      )
    )
} ## end func


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
