


template <class Type>
Type ssbi(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, int i){
  int stateDimN=logN.dim[0];
  Type ssb=0;
  for(int j=0; j<stateDimN; ++j){
    if(conf.keyLogFsta(0,j)>(-1)){
      ssb += exp(logN(j,i))*exp(-exp(logF(conf.keyLogFsta(0,j),i))*dat.propF(i,j)-dat.natMor(i,j)*dat.propM(i,j))*dat.propMat(i,j)*dat.stockMeanWeight(i,j);
    }else{
      ssb += exp(logN(j,i))*exp(-dat.natMor(i,j)*dat.propM(i,j))*dat.propMat(i,j)*dat.stockMeanWeight(i,j);
    }
  }
  return ssb;
}

template <class Type>
vector<Type> ssbFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF){
  int timeSteps=logF.dim[1];
  vector<Type> ssb(timeSteps);
  ssb.setZero();
  for(int i=0;i<timeSteps;i++){
    ssb(i)=ssbi(dat,conf,logN,logF,i);
  }
  return ssb;
}

template <class Type>
vector<Type> catchFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF){
  int len=dat.catchMeanWeight.dim(0);
  vector<Type> cat(len);
  cat.setZero();
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      Type z=dat.natMor(y,a-conf.minAge);
      if(conf.keyLogFsta(0,a-conf.minAge)>(-1)){
        z+=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y));
        cat(y)+=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y))/z*exp(logN(a-conf.minAge,y))*(Type(1.0)-exp(-z))*dat.catchMeanWeight(y,a-conf.minAge);
      }
    }
  }
  return cat;
}

template <class Type>
vector<Type> varLogCatchFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, paraSet<Type> par){
  int len=dat.catchMeanWeight.dim(0);
  vector<Type> cat(len);
  cat.setZero();
  vector<Type> varLogCat(len);
  varLogCat.setZero();
  Type CW=0;
  Type Ca=0;
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      Type z=dat.natMor(y,a-conf.minAge);
      if(conf.keyLogFsta(0,a-conf.minAge)>(-1)){
        z+=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y));
        CW=dat.catchMeanWeight(y,a-conf.minAge);
        Ca=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y))/z*exp(logN(a-conf.minAge,y))*(Type(1.0)-exp(-z));
        cat(y)+=Ca*CW;
        varLogCat(y)+=exp(2.0*par.logSdLogObs(conf.keyVarObs(0,a-conf.minAge)))*CW*CW*Ca*Ca; 
      }
    }
    varLogCat(y)/=cat(y)*cat(y);
  }
  return varLogCat;
}

template <class Type>
vector<Type> landFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF){
  int len=dat.landMeanWeight.dim(0);
  vector<Type> land(len);
  land.setZero();
  Type LW=0;
  Type LF=0;
  Type LWLF=0;
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){
      Type z=dat.natMor(y,a-conf.minAge);
      if(conf.keyLogFsta(0,a-conf.minAge)>(-1)){
        z+=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y));
        LW=dat.landMeanWeight(y,a-conf.minAge);
        if(LW<0){LW=0;}
	LF=dat.landFrac(y,a-conf.minAge);
        if(LF<0){LF=0;}
        LWLF=LW*LF;
        land(y)+=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y))/z*exp(logN(a-conf.minAge,y))*(Type(1.0)-exp(-z))*LWLF;
      }
    }
  }
  return land;
}

template <class Type>
vector<Type> varLogLandFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF, paraSet<Type> par){
  int len=dat.landMeanWeight.dim(0);
  vector<Type> land(len);
  land.setZero();
  vector<Type> varLogLand(len);
  varLogLand.setZero();
  Type LW=0;
  Type LF=0;
  Type LWLF=0;
  Type Ca=0;
  for(int y=0;y<len;y++){
    for(int a=conf.minAge;a<=conf.maxAge;a++){
      Type z=dat.natMor(y,a-conf.minAge);
      if(conf.keyLogFsta(0,a-conf.minAge)>(-1)){
        z+=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y));
        LW=dat.landMeanWeight(y,a-conf.minAge);
        if(LW<0){LW=0;}
        LF=dat.landFrac(y,a-conf.minAge);
        if(LF<0){LF=0;}
        LWLF=LW*LF;
        Ca=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y))/z*exp(logN(a-conf.minAge,y))*(Type(1.0)-exp(-z));
        land(y)+=Ca*LWLF;
        varLogLand(y)+=exp(2.0*par.logSdLogObs(conf.keyVarObs(0,a-conf.minAge)))*LWLF*LWLF*Ca*Ca;
      }
    }
    varLogLand(y)/=land(y)*land(y);
  }
  return varLogLand;
}

template <class Type>
vector<Type> fsbFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN, array<Type> &logF){
  int len=dat.catchMeanWeight.dim(0);
  vector<Type> fsb(len);
  fsb.setZero();
  Type sumF;
  for(int y=0;y<len;y++){  // calc logfsb
    sumF=Type(0);
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      if(conf.keyLogFsta(0,a-conf.minAge)>(-1)){
        sumF+=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y));
      }
    }
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      Type z=dat.natMor(y,a-conf.minAge);
      if(conf.keyLogFsta(0,a-conf.minAge)>(-1)){
        z+=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y));
        fsb(y)+=(exp(logF(conf.keyLogFsta(0,a-conf.minAge),y))/sumF)*exp(logN(a-conf.minAge,y))*exp(-Type(0.5)*z)*dat.catchMeanWeight(y,a-conf.minAge);
      }
    }
  }
  return fsb;
}

template <class Type>
vector<Type> tsbFun(dataSet<Type> &dat, confSet &conf, array<Type> &logN){
  int timeSteps=logN.dim[1];
  vector<Type> tsb(timeSteps);
  tsb.setZero();  
  for(int y=0;y<timeSteps;y++){  
    for(int a=conf.minAge;a<=conf.maxAge;a++){  
      tsb(y)+=exp(logN(a-conf.minAge,y))*dat.stockMeanWeight(y,a-conf.minAge);
    }
  }
  return tsb;
}

template <class Type>
vector<Type> rFun(array<Type> &logN){
  int timeSteps=logN.dim[1];
  vector<Type> R(timeSteps);
  R.setZero();
  for(int y=0;y<timeSteps;y++){  
    R(y)=exp(logN(0,y));
  }
  return R;
}

template <class Type>
vector<Type> fbarFun(confSet &conf, array<Type> &logF){
  int timeSteps=logF.dim[1];
  vector<Type> fbar(timeSteps);
  fbar.setZero();
  for(int y=0;y<timeSteps;y++){  
    for(int a=conf.fbarRange(0);a<=conf.fbarRange(1);a++){  
      fbar(y)+=exp(logF(conf.keyLogFsta(0,a-conf.minAge),y));
    }
    fbar(y)/=Type(conf.fbarRange(1)-conf.fbarRange(0)+1);
  }
  return fbar;
}
