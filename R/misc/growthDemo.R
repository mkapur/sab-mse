## example growth module

## incremental growth
vonB <- function(Lnow, K, Linf){
  lthen = Lnow + (Linf - Lnow)*(1-exp(-K))
  return(lthen)
}
vonBmid <- function(Lnow, K, Linf){
  lthen = Lnow + (Linf - Lnow)*(1-exp(-0.5*K))
  return(lthen)
}
Plus_LI <- function(Natage, Lmid, Lnow){
  lthen0 <- Natage[1]*Lmid[ages]+Natage[2]*vonB(Lnow, K,Linf)
  lthen1 <- sum(Natage)
  lthen = lthen0/lthen1
  return(lthen)
  
}
Linit <- function()
## setup tester params
ages = 95
Lstart = 4.5 ## length at age 1
K = 0.3
Linf = 150

dat <- matrix(NA, ncol = ages, nrow = 100)
row.names(dat) = seq(0.5,50,0.5)
colnames(dat) <- 1:ages
dat[1,1] <- Lstart

for(i in 1:(nrow(dat)-2)){ ## loop years
  for(j in 2:ncol(dat)){ ## loop ages
    dat[1,j] <- Linf + (dat[1,j-1]- Linf)*exp(-K)
  }
  for(j in 1:ncol(dat)){
    dat[i+1,j] <- vonBmid(Lnow = dat[i,j],K, Linf)
    dat[i+2,j] <- vonB(Lnow = dat[i,j], K, Linf)
  }
  }


dat %>%
  reshape2::melt() %>%
ggplot(., aes(x = Var2, y = value, color = Var1, group = Var1)) +
  geom_line()


dat <- data.frame("age" = 1:ages, 
                  
                  "len_y1" = c(Lstart, rep(NA, ages-1)),
                  "len_y0.5" = c(Lstart, rep(NA, ages-1)))


for(i in 2:(ages)){
  ## do them all
  dat$len_y1[i] <-  vonB(Lnow = dat[i-1,"len_y1"],K, Linf) ## start of next year
  dat$len_y0.5[i] <-  vonBmid(Lnow = dat[i-1,"len_y1"],K, Linf) ## midpoint of this year
  if(i != ages){  ## plus group
    dat$len_y1[i] <- Plus_LI(Natage = c(100,96), Lmid =   dat$len_y0.5[i],
                             Lnow = dat[i-1,"len"])
  }
  
  
}



## get plus group growht