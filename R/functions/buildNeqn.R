buildNeqn <- function(df){
  
  nage = df$nage ## #just trying ages 5:7 for now
  nspace =  df$nspace
  Ndim <- nage*nspace
  S <- matrix(0, nrow = Ndim, ncol = Ndim)
  for(i in 1:nspace){
    srt = ifelse(i == 1,1, 1+71*(i-1))
    #1+71*(i-1)
    stp = 71*(i)
    # cat(srt," ",stp,"\n")
    diag(S[srt:stp,srt:stp]) <- exp(-df$parms$mort_k[df$phi_ik2[i]+1]/2)
    #   diag(S[,srt:stp]) <- exp(-df$parms$mort_k[df$phi_ik2[i]+1]/2)
   
  }
  # diag(S) <- exp(-0.2/2)
  # cat( diag(S),"\n")
  # S[1:25,1:25]
  ## Grwoth aka age transition
  X <- matrix(0, nrow = Ndim, ncol = Ndim) 
  ## fill lower diag with 1s
  for(i in 1:Ndim){
    if(i %% nage != 0) {X[i+1,i] <- 1}#;  cat(i,"\t",i+1,"\n")}
  }
  for(i in seq(nage, Ndim,nage) ) X[i,i] <- 1 ## accumulate at plus group age 
  # colSums(X)
  
  I <- matrix(0,ncol=Ndim,nrow=Ndim); diag(I) = 1
  
  A_fem <- A_mal <- matrix(0, nrow = Ndim, ncol = Ndim)
  ## easier to make simple space x age matrix for each, then stich
  for(a in seq_along(1:nage)){ ## loop ages 5:nage
    for(i in 1:nspace){
      for(j in 1:nspace){
        ## working with one source area at a time, make a vertical vector for the age at hand
        ## the position (row) for the infill is a recurisve function of the age and block position
        ## tis is the same position to be used to index columns of final product
        colpos = a+(i-1)*nage ## doens't change with J; position in vector
        rowpos = a+(j-1)*nage
        # cat(a,"\t",i,"\t",j,"\t rowpos = ",rowpos,"\t colpos = ",colpos,'\n')
        A_fem[rowpos,colpos] <- round(df$X_ijas[i,j,a,1],2)
        A_mal[rowpos,colpos] <-  round(df$X_ijas[i,j,a,2],2)
      }
    }
  }
  A <- A_fem

  H <- matrix(0,ncol=Ndim,nrow=Ndim); diag(H) = 1
  # Specify the H matrix, which is literally the proportion remaining post-fishing
  # H[nage,nage] <- 1.0 - 0.2/0.8*FF ## 1- selex * f slide 7
  H[nage,nage] <- 1.0 - 0 ## 1- selex * f slide 7
  H[nage+nage,nage+nage] <- 1.0 - 0 ## next stage bin is just 1-FF
  # Multiply the matrices
  Mat1 <- S ## natural survival e^-M
  Mat2 <- (X %*% (A %*% (S %*% (H %*% S)))) 
  Neqn <- solve(I-Mat2) 
  return(Neqn)
  
  ## some testing stuff
  # R_init_temp <-  df$parms$logR_0k*df$tau_ki ## fixed R0 in each earea
  # r0vect <- matrix(0, nrow = 1, ncol = Ndim)
  # r0vect[seq(1,Ndim,nage)] <- exp(R_init_temp[R_init_temp !=0])
  # r0vect = c(r0vect)
  # neqnm <- matrix(Neqn%*% r0vect, ncol = 6, nrow = length(0:70)) %>%
  #   data.frame(.)
  # neqnm %>%  mutate(age = 0:70) %>%
  #   reshape2::melt(id = 'age') %>%
  #   ggplot(., aes(x = age, y = value, color = variable)) +
  #   ggsidekick::theme_sleek() +
  #   geom_line(lwd = 1.1) +
  #   scale_color_manual(values = rev(subareaPal),labels =  dimnames(df$X_ijas)[[1]]) +
  #   facet_wrap(~variable, scales = 'free_y' )
}
