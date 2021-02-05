getSelec2 <- function(sex, flt_idx, selP, selShape, selType,blk, fltType = 'fish'){

  age = df$age

  if(selType ==0){
    if(selShape == -1){
      selec = 1
    } else if(selShape==0){
      selec = 1/ ( 1 + exp(-log(19) * (age -  selP[flt_idx,1,blk,sex]) / ( selP[flt_idx,2,blk,sex] -
                                                                          selP[flt_idx,1,blk,sex])))
    } else if(selShape==1){
      cat('not ready for age sel shape == 1')
      # // Logistic with a50 and slope, where  selP(flt_idx,0,0,s) = a50 and  selP(flt_idx,1,0,s) = slope.
      # //  *This is the preferred logistic parameterization b/c it reduces parameter correlation*
      #   for (int a= 0; a < nage; a++){
      #     fsh_slx_yafs(i,a,flt_idx,s)  = Type(1.0) / ( Type(1.0) + exp( Type(-1.0) *
      #                                                                      selP(flt_idx,1,0,s) * (a -  selP(flt_idx,0,0,s)) ) );
      #   } // end ages
      # break;
    } else if(selShape==2){
      # // Dome Normal with alpha (mean) and beta (sd)
      selec  = exp(-(0.5 * (age -    selP[flt_idx,1,blk,sex])/  selP[flt_idx,2,blk,sex])^2);
    } else if(selShape==3){
      selec0 =    age ^( selP[flt_idx,1,blk,sex] - 1) * exp(-age/selP[flt_idx,2,blk,sex])
      selec <- selec0/ max(selec0)
    } ## end selshape for age sel
  } else if(selType ==1){
    ## get mla_yais for fleet x area
   if(fltType == 'fish') i = which(df$phi_if_fish[flt_idx,]==1)[1]
   
   if(fltType == 'surv') {
     temp_phi_if_surv <- rbind(df$phi_if_surv,df$phi_if_acomp)
     rown <- which(rownames(temp_phi_if_surv) ==  dimnames(df$parms$log_srv_slx_pars)[[1]][flt_idx])
     i = which(temp_phi_if_surv[rown,]==1)[1]
   }
    len = df$mla_yais[df$yRun,,i,sex]
    if(selShape == -1){
      selec = 1
    } else if(selShape==0){
      selec = 1/ ( 1 + exp(-log(19) * (len -  selP[flt_idx,1,blk,sex]) / ( selP[flt_idx,2,blk,sex] -
                                                                           selP[flt_idx,1,blk,sex])))
    } else if(selShape==1){
      cat('not ready for len sel shape == 1')
      # // Logistic with a50 and slope, where  selP(flt_idx,0,0,s) = a50 and  selP(flt_idx,1,0,s) = slope.
      # //  *This is the preferred logistic parameterization b/c it reduces parameter correlation*
      #   for (int a= 0; a < nlen; a++){
      #     fsh_slx_yafs(i,a,flt_idx,s)  = Type(1.0) / ( Type(1.0) + exp( Type(-1.0) *
      #                                                                      selP(flt_idx,1,0,s) * (a -  selP(flt_idx,0,0,s)) ) );
      #   } // end lens
      # break;
    } else if(selShape==2){
      # // Dome Normal with alpha (mean) and beta (sd)
      selec  = exp(-(0.5 * (len -    selP[flt_idx,1,blk,sex])/  selP[flt_idx,2,blk,sex])^2);
    } else if(selShape==3){
      ## k is first parm, theta is second (shape, scale)
      ## the mean is given by k*theta
      selec0 =    len^(selP[flt_idx,1,blk,sex] - 1) * exp(-len/selP[flt_idx,2,blk,sex])
      selec <- selec0/max(selec0)
    } ## end selshape for len sel
  } ## end leng or age sel
  # if(max(selec) > 1.2) cat('selec high')
  
  return(selec)
}
