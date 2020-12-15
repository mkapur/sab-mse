// Hybrid F tuning inputs & temp storage
Type v2 = 30;
for(int fish_flt =0;fish_flt<(nfleets_fish);fish_flt++){
  if(catch_yf_obs(y,fish_flt+1) != Type(-1.0)){
    // std::cout << fish_flt << " F TUNING" << "\n";
    // std::cout << y << 	" " << fish_flt << 	" THIS # IS NOT -1 " << catch_yf_obs(y,fish_flt+1) << std::endl;
    Type denom = 0;
    for(int s=0;s<nsex;s++){
      for(int i=0;i<(nspace);i++){
        for(int a=0;a<(nage);a++){
          switch(selType_fish(fish_flt)){
            case 0: // age sel
            denom += phi_if_fish(fish_flt,i)*
              fsh_slx_yafs(y,a,fish_flt,s)*
              N_yais_mid(y,a,i,s)*
              wtatlen_kab(phi_ik2(i),0)*
              pow( mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1));//+
              // catch_yf_obs(y,fish_flt+1);
            break;
            case 1: // length sel
            denom += phi_if_fish(fish_flt,i)*
              fsh_slx_yafs(y, mla_yais(y,a,i,s),fish_flt,s)*
              N_yais_mid(y,a,i,s)*
              mla_yais(y,a,i,s)*
              wtatlen_kab(phi_ik2(i),0)*
              pow( mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1));//+
              // std::cout << y  <<"\t area"<< i <<"\t flt"<<  fish_flt << "\t age" << a <<
              //   "\t aL^b \t" <<  wtatlen_kab(phi_ik2(i),0)*pow( mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1))  << "\n";
            // std::cout << y  <<"\t area"<< i <<"\t flt"<<  fish_flt << "\t age" << a <<
              //   "\t N_yais_mid(y,a,i,s) \t" <<  N_yais_mid(y,a,i,s)  << "\n";
            // std::cout << y  <<"\t area"<< i <<"\t flt"<<  fish_flt << "\t age" << a <<
              //   "\t  phi*fsh_slx_yafs(y,a,i,s) \t" <<   
              //     phi_if_fish(fish_flt,i)*fsh_slx_yafs(y, mla_yais(y,a,i,s),fish_flt,s)  << "\n";
            
            // catch_yf_obs(y,fish_flt+1);
            break;
          } // end selType_fish
        } // end age
        // std::cout << y  <<"\t area"<< i <<"\t flt"<<  fish_flt << "\t first loop denom \t" << denom  << "\n";
      } // end space
    } // end sex
    F1_yf(y,fish_flt,0) = catch_yf_obs(y, fish_flt+1)/(denom+ catch_yf_obs(y,fish_flt+1));
    Type latest_guess = F1_yf(y,fish_flt,0);
    // std::cout << y  <<"\t  catch_yf_obs = " <<  catch_yf_obs(y, fish_flt+1)  << "\n";
    // std::cout << y  <<"\t" <<fish_flt<<"\t  pre-iter denom  = " <<  denom  << "\n";
    // std::cout << y  <<"\t  pre-iter F1_yf(y,fish_flt,k) aka latest guess  = " <<  latest_guess  << "\n";
    // k iterations
    for(int k=1;k<(niter+1);k++){
      // std::cout << y << "\t" << k <<  "\t" << fish_flt <<  "\t" << "doing k iters" << "\n";
      // modify the guess Eq 20
      Type term0 = 1/(1+exp(v2*( latest_guess - v1)));
      Type term1 = latest_guess*term0;
      Type term2 = v1*(1-term0);
      F1_yf(y,fish_flt,k) = -log(1-(term1+term2));
      // std::cout << y << "\t" << k << "\t within iter F1_yf(y,fish_flt,k)  = " <<   F1_yf(y,fish_flt,k)  << "\n";
      Z_a_TEMP.setZero();
      for(int i=0;i<(nspace);i++){
        for(int a=0;a<(nage);a++){
          for(int s=0;s<nsex;s++){
            switch(selType_fish(fish_flt)){
              case 0: // age sel
              for(int s=0;s<nsex;s++){
                Z_a_TEMP[a] += fsh_slx_yafs(y, a, fish_flt, s)*F1_yf(y,fish_flt,k) + mort_k(phi_ik2(i));
              } // end sex for z a temp
              catch_afk_TEMP(a,fish_flt,k) +=
                F1_yf(y,fish_flt,k)/Z_a_TEMP[a]*
                (1-exp(-Z_a_TEMP[a]))*
                phi_if_fish(fish_flt,i)*
                fsh_slx_yafs(y,a,fish_flt,s)*
                N_yais_mid(y,a,i,s) *
                wtatlen_kab(phi_ik2(i),0) *
                pow(mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1));
              break;
              case 1: // length sel
              for(int s=0;s<nsex;s++){
                Z_a_TEMP[a] += fsh_slx_yafs(y, mla_yais(y,a,i,s), fish_flt, s)*F1_yf(y,fish_flt,k) + mort_k(phi_ik2(i));
              } // end sex for z a temp
              catch_afk_TEMP(a,fish_flt,k) +=
                F1_yf(y,fish_flt,k)/Z_a_TEMP[a]*
                (1-exp(-Z_a_TEMP[a]))*
                phi_if_fish(fish_flt,i)*
                fsh_slx_yafs(y,mla_yais(y,a,i,s),fish_flt,s)*
                N_yais_mid(y,a,i,s)*
                mla_yais(y,a,i,s)*
                wtatlen_kab(phi_ik2(i),0)*
                pow(mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1));
              break;
            } // end selType_fish
          } // end sex
          // std::cout << y << "\t" << k << "\t" << "\t  catch_afk_TEMP[a]  = " << catch_afk_TEMP(a,fish_flt,k) << "\n";
        } // end age
      } // end space
      // std::cout << y << "\t" << k << "\t" << "\t  Z_a_TEMP[0]  = " << Z_a_TEMP(0) << "\n";
      // std::cout << y << "\t" << k << "\t" << "\t  Z_a_TEMP[a+]  = " << Z_a_TEMP(nage-1) << "\n";
      // std::cout << y << "\t" << k << "\t" << "\t  catch_afk_TEMP[a+]  = " << catch_afk_TEMP(nage-1,fish_flt,k) << "\n";
      // std::cout << y << "\t" << k << "\t" << "\t  catch_afk_TEMP[a0]  = " << catch_afk_TEMP(0,fish_flt,k) << "\n";
      // std::cout << y << "\t" << k << "\t" << "\t  catch_yf_obs(y,fish_flt+1) = " << catch_yf_obs(y,fish_flt+1) << "\n";
      
      vector<Type>Adj(niter+1);
      Adj.setZero();
      Type denom2 = 0.0;
      for(int a=0;a<(nage);a++){
        // buffer in case catch is straight up zero for any age, leading to NaN
        denom2 += catch_afk_TEMP(a,fish_flt,k);
      }
      Adj(k) += catch_yf_obs(y,fish_flt+1)/(denom2+1e-9);
      // std::cout << y << "\t" << k << "\t" << fish_flt << "\t  Adjk  = " <<   Adj(k)  << "\n";
      
      // Get new Z given ADJ - need to add discard here and switch selex
      Z_a_TEMP2.setZero();
      for(int a=0;a<(nage);a++){
        for(int s=0;s<nsex;s++){
          Z_a_TEMP2(a) += Adj(k)  *
            fsh_slx_yafs(y, a, fish_flt, s) * F1_yf(y, fish_flt, k) +
            0.125; // standin mean M
        } // end sex
      } // end age
      // std::cout << y << "\t" << k << "\t" << "\t  Z_a_TEMP2[0]  = " <<   Z_a_TEMP2(0)  << "\n";
      
      // Second Guess for F (EQ 24)
      Type denom = 0;
      for(int s=0;s<nsex;s++){
        for(int i=0;i<(nspace);i++){
          for(int a=0;a<(nage);a++){
            switch(selType_fish(fish_flt)){
              case 0: // age sel
              denom += phi_if_fish(fish_flt,i)*
                fsh_slx_yafs(y,a,fish_flt,s)*
                N_yais_mid(y,a,i,s)*
                wtatlen_kab(phi_ik2(i),0)*
                pow(mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1))*
                (1-exp(-Z_a_TEMP2(a))) * (F1_yf(y,fish_flt,k)/(Z_a_TEMP2(a)));
              break;
              case 1: // length sel
              denom += phi_if_fish(fish_flt,i)*
                fsh_slx_yafs(y,mla_yais(y,a,i,s),fish_flt,s)*
                N_yais_mid(y,a,i,s)*
                mla_yais(y,a,i,s)*
                wtatlen_kab(phi_ik2(i),0)*
                pow(mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1))*
                (1-exp(-Z_a_TEMP2(a))) * (F1_yf(y,fish_flt,k)/(Z_a_TEMP2(a)));
              break;
            } // end selType_fish
          } // end age
          // std::cout << y << "\t" << k << "\t" << i << "\t  space accum denom for F2  = " <<   denom  << "\n";
        } // end space
      } // end sex
      F2_yf(y, fish_flt, k) = catch_yf_obs(y, fish_flt+1)/denom;
      // std::cout << y << "\t" << k << "\t" << "\t  denom for F2  = " <<   denom  << "\n";
      
      // std::cout << y << "\t" << k << "\t    F2_yf(y, fish_flt, k) = " <<    F2_yf(y, fish_flt, k)   << "\n";
      // Modify the guess again Eq 25
      term0 = 1/(1+exp(v2*( F2_yf(y,fish_flt,k )- v1*Fmax)));
      term1 = F2_yf(y,fish_flt,k)*term0;
      term2 = v1*(1-term0);
      F2_yf(y, fish_flt, k) = -log(1-(term1+term2));
      latest_guess =    F2_yf(y, fish_flt, k);
      // std::cout << y << "\t" << k << "\t latest guess (mod   F2_yf(y, fish_flt, k) )= " << latest_guess << "\n";
    } // end k hybrid F iterations
    
    // std::cout << y << "\t" << "end k iters" << "\n";
    // Define F, Z and predicted catches
    Freal_yf(y, fish_flt) = F2_yf(y, fish_flt, niter); //final as Freal_yf
    // annoying multi-loops for F in area  get total N exploitable by this fleet
    for(int a=0;a<(nage);a++){
      for(int i=0;i<(nspace);i++){
        for(int s=0;s<nsex;s++){
          N_avail_yf(y,fish_flt) += phi_if_fish(fish_flt, i)*N_yais_mid(y,a,i,s);
          // std::cout << y << "\t" << fish_flt << "\t N_avail_yf" <<  N_avail_yf(y,fish_flt) << "\n";
        } // end sex
      } // end nspace
    } // end nage
    // std::cout << y << "\t" << "end navail" << "\n";
    // get ratio of N in area & reweight F; will just return Freal and 1 for single-area fisheries
    for(int i=0;i<(nspace);i++){
      for(int a=0;a<(nage);a++){
        for(int s=0;s<nsex;s++){
          N_weight_yfi(y,fish_flt, i) += (phi_if_fish(fish_flt, i)* N_yais_mid(y,a,i,s))/N_avail_yf(y,fish_flt);
        } // end sex
      } // end age
      // std::cout << y << "\t" << fish_flt << "\t" << i << "\t" <<   N_weight_yfi(y,fish_flt, i) << "\n";
      F_area_yfi(y,fish_flt,i) = Freal_yf(y, fish_flt) * N_weight_yfi(y,fish_flt, i);
    } // end space
    // std::cout << y << "\t" << "end N_weight_yfi" << "\n";
    // add together for mgmt regions
    for(int m=1;m<(nmgmt_reg);m++){
      F_ym(y,m) += phi_fm(fish_flt,m)*Freal_yf(y, fish_flt);
    } // end mgmt regions
    // generate predicted catches
    // std::cout << y << "\t"<<fish_flt << "\t Freal + half mort (Zreal) \t" << Freal_yf(y, fish_flt)+0.2/2 << "\n";
    for(int a=0;a<(nage);a++){
      for(int i=0;i<(nspace);i++){
        for(int s=0;s<nsex;s++){
          switch(selType_fish(fish_flt)){
            case 0: // age sel
            
            // Zreal_ya(y,a) += Freal_yf(y, fish_flt) + mort_k(phi_ik2(i))/2;
            // Zreal_yai(y,a,i) += F_area_yfi(y, fish_flt,i) + mort_k(phi_ik2(i))/2;
            // catch_yaf_pred(y,a,fish_flt) +=
              //   Freal_yf(y, fish_flt)/ Zreal_ya(y,a) *
              //   (1-exp(- Zreal_ya(y,a) ))*
              //   phi_if_fish(fish_flt,i)*
              //   fsh_slx_yafs(y,a,fish_flt,s)*
              //   N_yais_mid(y,a,i,s)*
              //   wtatlen_kab(phi_ik2(i),0)*
              //   pow(mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1));
            
            // instantaneous (midyear) version
            instF_yafs(y,a,fish_flt,s) = fsh_slx_yafs(y,a,fish_flt,s)* Freal_yf(y, fish_flt);
            catch_yaf_pred(y,a,fish_flt) +=
              phi_if_fish(fish_flt,i)*
              instF_yafs(y,a,fish_flt,s) *
              N_yais_mid(y,a,i,s)*
              wtatlen_kab(phi_ik2(i),0)*
              pow(mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1));
            
            // catch_yaif_pred(y,a,i,fish_flt) += (F_area_yfi(y,fish_flt,i)/
                                                     //   ( Zreal_yai(y,a,i)))*(1-exp(- Zreal_yai(y,a,i)  ))*
              //     phi_if_fish(fish_flt, i)*
              //     fsh_slx_yafs(y,a,fish_flt,s)*
              //     N_yais_mid(y,a,i,s)*
              //     wtatlen_kab(phi_ik2(i),0)*
              //     pow(mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1));
            
            // instantaneous (midyear) version
            catch_yaif_pred(y,a,i,fish_flt) +=
              phi_if_fish(fish_flt, i)*
              fsh_slx_yafs(y,a,fish_flt,s)*
              F_area_yfi(y,fish_flt,i)*
              N_yais_mid(y,a,i,s)*
              wtatlen_kab(phi_ik2(i),0)*
              pow(mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1));
            // std::cout << y << "\t" << fish_flt <<"\t catch_yaf_pred \t" <<    catch_yaf_pred(y,a,fish_flt)<< "\n";
            break;
            case 1: // length sel
            // Zreal_ya(y,a) += Freal_yf(y, fish_flt) + mort_k(phi_ik2(i))/2;
            // Zreal_yai(y,a,i) += F_area_yfi(y, fish_flt,i) + mort_k(phi_ik2(i))/2;
            
            // catch_yaf_pred(y,a,fish_flt) +=
              //   Freal_yf(y, fish_flt)/ Zreal_ya(y,a) *
              //   (1-exp(- Zreal_ya(y,a) ))*
              //   phi_if_fish(fish_flt,i)*
              //   fsh_slx_yafs(y,mla_yais(y,a,i,s),fish_flt,s)*
              //   N_yais_mid(y,a,i,s)*
              //   mla_yais(y,a,i,s)*
              //   wtatlen_kab(phi_ik2(i),0)*
              //   pow(mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1));
            
            // instantaneous (midyear) version
            instF_yafs(y,a,fish_flt,s) = fsh_slx_yafs(y,mla_yais(y,a,i,s),fish_flt,s)* Freal_yf(y, fish_flt);
            catch_yaf_pred(y,a,fish_flt) +=
              phi_if_fish(fish_flt,i)*
              instF_yafs(y,a,fish_flt,s)*
              N_yais_mid(y,a,i,s)*
              mla_yais(y,a,i,s)*
              wtatlen_kab(phi_ik2(i),0)*
              pow(mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1));
            
            // catch_yaif_pred(y,a,i,fish_flt) += (F_area_yfi(y,fish_flt,i)/
                                                     //   ( Zreal_yai(y,a,i)))*(1-exp(- Zreal_yai(y,a,i)  ))*
              //     phi_if_fish(fish_flt,i)*
              //     fsh_slx_yafs(y, mla_yais(y,a,i,s),fish_flt,s)*
              //     N_yais_mid(y,a,i,s)*
              //     mla_yais(y,a,i,s)*
              //     wtatlen_kab(phi_ik2(i),0)*
              //     pow(mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1));
            
            // instantaneous (midyear) version
            catch_yaif_pred(y,a,i,fish_flt) +=
              phi_if_fish(fish_flt,i)*
              F_area_yfi(y,fish_flt,i)*
              fsh_slx_yafs(y, mla_yais(y,a,i,s),fish_flt,s)*
              N_yais_mid(y,a,i,s)*
              mla_yais(y,a,i,s)*
              wtatlen_kab(phi_ik2(i),0)*
              pow(mla_yais(y,a,i,s),wtatlen_kab(phi_ik2(i),1));
            
            
            break;
          } // end selType_fish
        } // end sex
        catch_yfi_pred(y,fish_flt,i) += catch_yaif_pred(y,a,i,fish_flt);
      } // end space
      catch_yf_pred(y,fish_flt) += catch_yaf_pred(y,a,fish_flt);
    } // end age
  } // end -1 NA trap
  // std::cout << y << "\t" << fish_flt <<"\t catch_yf_pred \t" << catch_yf_pred(y,fish_flt)<< "\n";
} // end nfleets_fish