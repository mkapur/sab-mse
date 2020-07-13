catch_afk_TEMP <- array(0, dim = c(nage, nfleets_fish, niter+1)) ## storage for intermediate guesses

## change this to be FI-indexed, remove summations over nspace

Adj <- NULL; Z_a_TEMP <- Z_a_TEMP2 <- NULL

for(fish_flt in 1:nfleets_fish){
  catch_yaf_pred[y,,fish_flt] <- catch_yf_pred[y,fish_flt] <-  0
  ## make an initial guess for Ff using obs catch - need to update selex whihc is 1.0 now
  denom = 0
  for(i in 1:nspace){ 
    denom <- denom +   (phi_if_fish[fish_flt, i] * 
                          sum(N_yai_beg[y,,i])*
                          sum(wage_catch[,y]) * 1.0 + catch_yf_obs[y, fish_flt+1])
  }
  # F1_yf[y, fish_flt, 1] <-     F1_yf[y, fish_flt, 1] + catch_yf_obs[y, fish_flt+1]/denom
  F1_yf[y, fish_flt, 1] <-    catch_yf_obs[y, fish_flt+1]/denom
  
  latest_guess <-    F1_yf[y, fish_flt, 1]
  ## iterative tuning
  for(k in 2:(niter+1)){
    ## modify the guess Eq 20
    term0 = 1/(1+exp(v2*( latest_guess - v1)));
    term1 = latest_guess*term0;
    term2 = v1*(1-term0);
    F1_yf[y,fish_flt,k] = -log(1-(term1+term2))
    
    # Predicted catches @ F Eq 21; need to add SELEX
    for(i in 1:nspace){ 
      for(a in 1:nage){
        Z_a_TEMP[a] <- F1_yf[y,fish_flt,k] + M[a]
        
        catch_afk_TEMP[a,fish_flt,k] <-     catch_afk_TEMP[a,fish_flt,k] +
          (F1_yf[y,fish_flt,k]/(Z_a_TEMP[a]))*(1-exp(-Z_a_TEMP[a]))*
          phi_if_fish[fish_flt, i]* 
          N_yai_beg[y,a,i]*1.0*
          wage_catch[a,i]
        
        
        
      } ## end ages
    } ## end space
    
    
    
    ## Calc Adj Eq 22
    Adj[k] <- catch_yf_obs[y, fish_flt+1]/sum(catch_afk_TEMP[,fish_flt,k])
    
    ## Get new Z given ADJ - need to add discard and selex here
    for(a in 1:nage) Z_a_TEMP2[a] <-  Adj[k]*F1_yf[y,fish_flt,k] +  M[a]
    
    ## Second Guess for F (EQ 24)
    denom = 0
    for(i in 1:nspace){ 
      for(a in 1:nage){
        denom <- denom + phi_if_fish[fish_flt, i] * 
          N_yai_beg[y,a,i]*
          wage_catch[a,y] * 1.0 *(1-exp(-Z_a_TEMP2[a])) * (F1_yf[y,fish_flt,k]/(Z_a_TEMP2[a]))
      }
    }
    # F2_yf[y, fish_flt, k] <- F2_yf[y, fish_flt, k-1] + catch_yf_obs[y, fish_flt+1]/denom
    F2_yf[y, fish_flt, k] <- catch_yf_obs[y, fish_flt+1]/denom
    
    ## Modify the guess again Eq 25
    term0 = 1/(1+exp(v2*( F2_yf[y,fish_flt,k] - v1)));
    term1 = F2_yf[y,fish_flt,k]*term0;
    term2 = v1*(1-term0);
    F2_yf[y,fish_flt,k] = -log(1-(term1+term2))
    # cat(F2_yf[y,fish_flt,k],"\n")
    latest_guess <-     F2_yf[y,fish_flt,k]
    
    
  } ## end hybrid F iterations
  
  ## Define F, Z and predicted catches ----
  Freal_yf[y, fish_flt] <- latest_guess # F2_yf[y,fish_flt,niter] ## final as Freal_yf
  for(i in 1:nspace){ 
    for(a in 1:nage){
      Zreal_ya[y,a] <-   Freal_yf[y, fish_flt] + M[a]
      catch_yaf_pred[y,a,fish_flt] <- catch_yaf_pred[y,a,fish_flt] +
        (Freal_yf[y, fish_flt]/(Zreal_ya[y,a]))*(1-exp(-Zreal_ya[y,a]))*
        phi_if_fish[fish_flt, i]* 
        N_yai_beg[y,a,i]*1.0*
        wage_catch[a,i]
    } ## end ages for predicted catch
  } ## end nspace for predicted catch
  # if(catch_yaf_pred[y,a,fish_flt] == 0) stop(y," ",fish_flt)
  catch_yf_pred[y,fish_flt] <- sum(catch_yaf_pred[y,,fish_flt])
} ## end fishery fleets
# cat( Freal_yf[y, fish_flt],fish_flt,"\n")
