#include <TMB.hpp>

// Do not transform sigma parameters in the likelhood
// Include variable transformation with bph in the likelihood
//get rid of beta0  
//If you want the non-i model to work you have to estimate a random effect for every observation

/* Parameter transform */
  template <class Type>
  Type f(Type x, int mod){
    if(mod==1)
      return Type(1)/(Type(1) + exp(x));
    if(mod==2)
      return atan(x)/(PI/2);
    }
  //Type f(Type x){return atan(x)/(PI/2);}
  //Type f(Type x){return Type(2)/(Type(1) + exp(-Type(2) * x)) - Type(1);}
  //Type f(Type x){return (exp(x)-exp(-x))/(exp(x)+exp(-x));}
  
  template <class Type>
  Type xt(Type par_xinf, Type par_x0, Type par_a, Type par_k)
  {
    Type tmpX0 = ((par_xinf-par_x0)/par_xinf);
    return par_xinf*(1-tmpX0*exp(-par_k*(par_a)));
  }

  template <class Type>
  Type gomp_xt(Type par_xinf, Type par_x0, Type par_a, Type par_k)
  {
    Type tmpX0 = ((log(par_xinf) - log(par_x0))/log(par_xinf));
    return exp(log(par_xinf)*(1-tmpX0*exp(-par_k*(par_a))));
  }

  template <class Type>
  Type logistic_xt(Type par_xinf, Type par_x0, Type par_a, Type par_k)
  {
    //See Schoener and Schoener 1978
    Type tmpb = 1 - par_x0/par_xinf;
    return ((par_xinf*(1-tmpb)*exp(par_k*par_a))/(tmpb+(1-tmpb)*exp(par_k*par_a)));
  }
  
  template <class Type>
  Type dxdt(Type x, Type par_xinf, Type delT, Type par_k)
  {return x + (par_xinf - x)*(1-exp(-par_k*delT));}

  template <class Type>
  Type gomp_dxdt(Type x, Type par_xinf, Type delT, Type par_k)
  {return exp(log(x) + (log(par_xinf) - log(x))*(1-exp(-par_k*delT)));}
  
  template <class Type>
  Type logistic_dxdt(Type x, Type par_xinf, Type delT, Type par_k)
  {
    //See Schoener and Schoener 1978
    return (1/(1/x + (1/par_xinf - 1/x)*(1-exp(-1.*par_k*delT))));
  }

  template <class Type>
  Type inv_VB(Type x, Type par_xinf, Type par_x0, Type par_k)
  {
    Type tmpX0 = ((par_xinf-par_x0)/par_xinf);
    return log(-1*(x/par_xinf - 1)/tmpX0)/(-par_k);}

  template <class Type>
  Type inv_gomp(Type x, Type par_xinf, Type par_x0, Type par_k)
  {
    Type tmpX0 = ((log(par_xinf) - log(par_x0))/log(par_xinf));
    return log(-1*((log(x)/log(par_xinf) - 1)/tmpX0))/(-par_k);}
  
  template <class Type>
  Type inv_logistic(Type x, Type par_xinf, Type par_x0, Type par_k)
  {
    Type tmpb = 1 - par_x0/par_xinf;
    return log((par_xinf*tmpb - x*(1-tmpb))/(x*tmpb))/par_k;}
  
  template <class Type>
  Type g(Type x,Type par_b1)
  {
    return par_b1 * x;}
  
//Integrated model of growth for species with hard parts and mark-recapture
template<class Type>
Type objective_function<Type>::operator() ()
{
  DATA_INTEGER(ni); //Number of stranded turtles
  DATA_INTEGER(nia);//Number of stranded turtles without an annulus 
  DATA_IVECTOR(a0); //Does the turtle possess an annulus
  DATA_INTEGER(growthModel); //the type of growth model used
  DATA_INTEGER(rho_mod); //the type of growth model used
  DATA_VECTOR(nLAG); //Number of LAG for each turtle
  DATA_VECTOR(Xobs); //humerus observations
  DATA_VECTOR(Xmax); //maximum humerus diameter of stranded turtle
  DATA_VECTOR(SCLmax);
  DATA_VECTOR(del_xT); //Time at liberty
  DATA_IVECTOR(x_iID); //unique identifier for stranded turtles
  DATA_IVECTOR(x_ijID); //identifiers for the LAG of in each humerus
  DATA_IVECTOR(x_yrID); //year id for the different LAGs
  DATA_IVECTOR(x_stateID);
  DATA_IVECTOR(hij);
  DATA_INTEGER(hum_obs_model);
  DATA_INTEGER(mr_obs_model);
  DATA_IVECTOR(sex_ijID);

    
  DATA_VECTOR(Lobs);
  DATA_IVECTOR(tij);
  DATA_VECTOR(n_mr);
  DATA_VECTOR(del_mrT);
  DATA_IVECTOR(Lmr_iID); //unique identifier for stranded turtles
  DATA_IVECTOR(Lmr_ijID); //identifiers for the LAG of in each humerus
  DATA_VECTOR(Lvec);
  DATA_INTEGER(ni_mr);
  

  //humerus time-series simulation
  PARAMETER(lnLinf);   //Average max humerus
  PARAMETER(lnk);    //Average Brody growth coefficient
  PARAMETER(lnL0);
  PARAMETER(lnmu_a);//average age at stranding
  PARAMETER(lnmu_a_mr);//average age at stranding
  PARAMETER(lnsig_x);//observation error of humerus predictions
  PARAMETER(lnsig_Linfi);//standard deviation of Dinf random effects
  PARAMETER(lnsig_kij);//standard deviation of transient k-effects
  PARAMETER(lnsig_ki);//strandard deviation of individual k-effects
  PARAMETER(lnsig_A); //stanard deviation of age random effect
  PARAMETER(lnsig_A_mr); //stanard deviation of age random effect
  
  PARAMETER_VECTOR(re_Linfi);
  PARAMETER_VECTOR(re_ki);
  PARAMETER_VECTOR(re_kij);
  PARAMETER_VECTOR(re_ai);
  PARAMETER(lnbeta1);
  PARAMETER(lnsig_bph);
  
  //mark_recapture time-series simulation
  PARAMETER(lnsig_L);
  PARAMETER_VECTOR(re_Linfi_mr);
  PARAMETER_VECTOR(re_ki_mr);
  PARAMETER_VECTOR(re_kij_mr);
  PARAMETER_VECTOR(re_ai_mr);
  PARAMETER_VECTOR(sex_ki);
  PARAMETER_VECTOR(sex_Linfi);
  
  Type res = 0.;
  
  Type Linf = exp(lnLinf);   //Average max humerus
  Type L0 = exp(lnL0);
  Type k = exp(lnk);
  Type mu_a = exp(lnmu_a);
  Type mu_a_mr = exp(lnmu_a_mr);
  Type beta1 = exp(lnbeta1);

  Type sig_x = exp(lnsig_x);
  Type sig_Linfi = exp(lnsig_Linfi);
  Type sig_ki = exp(lnsig_ki);
  Type sig_kij = exp(lnsig_kij);
  Type sig_A = exp(lnsig_A);
  Type sig_A_mr = exp(lnsig_A_mr);
  
  Type sig_L = exp(lnsig_L);
  Type sig_bph = exp(lnsig_bph);
  
  int nij = Xobs.size();
  vector<Type> tmp_kij(nij-ni);
  vector<Type> k_ij_noYr(nij);
  vector<Type> k_ij(nij);
  vector<Type> Linf_ij(nij);
  vector<Type> L_ij(nij);
  vector<Type> a_ij(nij);


  // * These are for remaining humerus observations
  //bph likelihood
  Type nll7 = 0.;
  vector<Type> predDmax = beta1 * SCLmax;
  nll7 = -1*sum(dnorm(Xmax, beta1 * SCLmax,exp(lnsig_bph),true));

  Type X0 = L0; //This is meaningless except for some of the previous plotting functions used X0 instead of L0
  
  vector<Type> a_i(ni);
  int ii = 0;
  int icnt = 0;
  int jcnt = 0;
  int jcnt_mr = 0;
  Type nll1 = 0.;
  for(int ij=0; ij<nij; ij++){ //loop over all of the observations
    Linf_ij(ij) = Linf + (re_Linfi(x_iID(ij))); //Individual asymp. size
    if(sex_ijID(ij)>0) //Males are 1 and females are 2; zero is unknown
      Linf_ij(ij) += sex_Linfi(sex_ijID(ij)-1); //additive sex effect
    
    if(x_ijID(ij)==0){ //If it's the first observation
      if(a0(icnt)==0){ //If the first observation is not an annulus
        a_ij(ij) = exp(lnmu_a + re_ai(ii) -0.5*sig_A*sig_A); //age of first observation
        a_i(icnt) = exp(lnmu_a + re_ai(ii) -0.5*sig_A*sig_A);//This is just a variable to store the age of the first observation
        ii++;
      }
      if(a0(icnt)==1){ //Turtle with an annulus.
        a_ij(ij) = 0.75; //Avens et al. The first LAG is set to 0.75
        a_i(icnt) = 0.75; //Avens et al. The first LAG is set to 0.75
      }

      k_ij(ij) = k * exp(re_ki(x_iID(ij))); //growth coefficient for the first observation

      //sex effects for growth coefficeint
      if(sex_ijID(ij)>0)//male = 1, unknown = 2
        k_ij(ij) *= exp(sex_ki(sex_ijID(ij)-1));

      //Size at age of the first observations of the humerus
      if(growthModel==1)
        L_ij(ij) = xt(Linf_ij(ij), L0, a_ij(ij), k_ij(ij));
      if(growthModel==2)
        L_ij(ij) = gomp_xt(Linf_ij(ij), L0, a_ij(ij), k_ij(ij));
      if(growthModel==3)
        L_ij(ij) = logistic_xt(Linf_ij(ij), L0, a_ij(ij), k_ij(ij));
      //increment counter
      icnt++;
    }//end if
    
    //subsequent observations of the same indiviudal
    if(x_ijID(ij)>0){ //If the increment ID is greater than 0
      a_ij(ij) = a_ij(ij-1) + del_xT(ij); //at the number of years between LAG increment         
      
      if(hij(ij)==0){
        k_ij(ij) = k * exp(re_ki(x_iID(ij)));
        if(sex_ijID(ij)>0) //males and unknown effects of sex
          k_ij(ij) *= exp(sex_ki(sex_ijID(ij)-1));
      }
      
      if(hij(ij)==1){ //This was an index that was added to aid in convergence, 
        k_ij(ij) = k * exp(re_ki(x_iID(ij)) + re_kij(jcnt)); 
        if(sex_ijID(ij)>0) //males and unknown effects of sex
          k_ij(ij) *= exp(sex_ki(sex_ijID(ij)-1)); //It has to be sex_ijID - 1, because males are 1 and unknowns are 2.
        jcnt++;
      }
      
      //size at age, after growth increment del_xT
      if(growthModel==1)
        L_ij(ij) = dxdt(L_ij(ij-1), Linf_ij(ij), del_xT(ij), k_ij(ij));
      if(growthModel==2)
        L_ij(ij) = gomp_dxdt(L_ij(ij-1), Linf_ij(ij), del_xT(ij), k_ij(ij));
      if(growthModel==3)
        L_ij(ij) = logistic_dxdt(L_ij(ij-1), Linf_ij(ij), del_xT(ij), k_ij(ij));
    }//end if
    
    //Likelihood
    if(hum_obs_model==1)
      nll1 -= dnorm(log(Xobs(ij)),log(beta1*L_ij(ij))-0.5*sig_x*sig_x, sig_x, true);
   }// end i
  
  int nij_mr = Lobs.size();
  vector<Type> Lpred(nij_mr);
  vector<Type> k_ij_mr(nij_mr);
  vector<Type> Linfi(ni_mr);
  vector<Type> a_ij_mr(nij_mr);
  
  icnt = 0;
  jcnt_mr = 0;
  Type nll8 = 0.;
  vector<Type> a_i_mr(ni_mr);
  //mark -recapture data
  for(int ij=0; ij<nij_mr; ij++)
  {
    Linfi(Lmr_iID(ij)) = Linf + (re_Linfi_mr(Lmr_iID(ij))); //Individual Linf
    Linfi(Lmr_iID(ij)) += sex_Linfi(1); //sex effect. All mark-recapture turtles are unknown sex. ID = 1

    if(Lmr_ijID(ij)==0){
      //growth coefficent of the first observation
      k_ij_mr(ij) = k*exp(re_ki_mr(Lmr_iID(ij))+sex_ki(1));
      //age of the first observation
      a_ij_mr(ij) = exp(lnmu_a_mr + re_ai_mr(Lmr_iID(ij)) -0.5*sig_A_mr*sig_A_mr);
      a_i_mr(icnt) = exp(lnmu_a_mr + re_ai_mr(Lmr_iID(ij)) -0.5*sig_A_mr*sig_A_mr);

      icnt++;
      //Size of the initial observation
      if(growthModel==1)
        Lpred(ij) = xt(Linfi(Lmr_iID(ij)), L0, a_ij_mr(ij), k_ij_mr(ij));
      if(growthModel==2)
        Lpred(ij) = gomp_xt(Linfi(Lmr_iID(ij)), L0, a_ij_mr(ij), k_ij_mr(ij));
      if(growthModel==3)
        Lpred(ij) = logistic_xt(Linfi(Lmr_iID(ij)), L0, a_ij_mr(ij), k_ij_mr(ij));
    }

    //Size of subsequent observations
    if(Lmr_ijID(ij)>0){
      if(tij(ij)!=1){//This index is created by the data file
        //growth coefficient with persistent effects
        k_ij_mr(ij) = k*exp(re_ki_mr(Lmr_iID(ij))+sex_ki(1));
      }
      if(tij(ij)==1){//This index is created by the data file
        //growth coefficent with transient effect
        k_ij_mr(ij) = k*exp(re_ki_mr(Lmr_iID(ij))+re_kij_mr(jcnt_mr)+sex_ki(1));
        jcnt_mr++;
      }
      //Age of subsequent observations
      a_ij_mr(ij) = a_ij_mr(ij-1) + del_mrT(ij);
      //Size at age at a funciton of del_mrT
      if(growthModel==1)
        Lpred(ij) = dxdt(Lpred(ij-1), Linfi(Lmr_iID(ij)), del_mrT(ij), k_ij_mr(ij));
      if(growthModel==2)
        Lpred(ij) = gomp_dxdt(Lpred(ij-1), Linfi(Lmr_iID(ij)), del_mrT(ij), k_ij_mr(ij));
      if(growthModel==3)
        Lpred(ij) = logistic_dxdt(Lpred(ij-1), Linfi(Lmr_iID(ij)), del_mrT(ij), k_ij_mr(ij));
    }
    
    if(mr_obs_model==1)
     nll8 -= dnorm(log(Lobs(ij)),log(Lpred(ij))-0.5*sig_L*sig_L, sig_L, true);
    //nll8 += pow(log(Lobs(ij))-log(Lpred(ij)),2);//
  }
  
  // nll1 -= sum(dnorm(log(Xobs),log(beta1*L_ij)-0.5*sig_x*sig_x, sig_x, true));
  // nll8 -= sum(dnorm(log(Lobs),log(Lpred)-0.5*sig_L*sig_L, sig_L, true));
  
  //Humerus likelihoods
  Type nll2 = 0.;
  Type nll3 = 0.;
  Type nll4 = 0.;
  Type nll5 = 0.;
  if(hum_obs_model==1){//deprecated boolean flag
    nll2 += -1.*sum(dnorm(re_Linfi,Type(0.),(sig_Linfi),true));
    nll3 += -1.*sum(dnorm(re_ki,Type(-0.5*sig_ki*sig_ki),sig_ki,true));
    nll4 += -1.*sum(dnorm(re_kij,Type(-0.5*sig_kij*sig_kij),sig_kij,true));
    nll5 += -1.*sum(dnorm(re_ai,0,exp(lnsig_A),true));
  }
  
  //Mark recapture data streams and random effects
  Type nll9 = 0.;
  Type nll10 = 0.;
  Type nll11 = 0.;
  Type nll12 = 0.;
  if(mr_obs_model==1){//deprecated boolean flag
    nll9 += -1.*sum(dnorm(re_ki_mr,Type(-0.5*sig_ki*sig_ki),sig_ki,true));
    nll10 += -1.*sum(dnorm(re_kij_mr,Type(-0.5*sig_kij*sig_kij),sig_kij,true));
    nll11 += -1.*sum(dnorm(re_Linfi_mr,Type(0.),(sig_Linfi),true));
    nll12 += -1.*sum(dnorm(re_ai_mr,0,(sig_A_mr),true));
  }
  
  res = nll1 + 
    nll2 + 
    nll3 + 
    nll4 + 
    nll5 + 
    nll7 +
    nll8 +
    nll9 +
    nll10 +
    nll11 +
    nll12;

  Type resFixed = nll1 + nll7 + nll8;
  
  //Estimated total age  
  vector<Type> estA = nLAG + a_i;

  vector<Type> Linf_i = Linf  + (re_Linfi);
  vector<Type> Linf_i_mr = Linf  + (re_Linfi_mr);
  
  vector<Type> k_i = k * exp(re_ki);
  vector<Type> k_i_mr = k * exp(re_ki_mr);

  REPORT(re_ki_mr);  
  REPORT(Linf);
  REPORT(Linf_i_mr);
  REPORT(k);
  REPORT(k_i);
  REPORT(k_i_mr);
  REPORT(Linf_i);
  REPORT(Linf_ij);
  REPORT(re_Linfi);
  REPORT(re_ki);
  REPORT(re_kij);
  REPORT(tmp_kij);
  REPORT(X0);
  REPORT(L0);
  REPORT(re_kij);
  REPORT(re_ki);
  REPORT(re_ai);
  REPORT(a_i);
  REPORT(estA);
  REPORT(mu_a);
  REPORT(sig_x);
  REPORT(a_ij_mr);
  REPORT(sig_Linfi);
  REPORT(mu_a_mr);
  REPORT(re_ai_mr);
  REPORT(a_i_mr);
  REPORT(sig_ki);
  REPORT(sig_kij);
  REPORT(jcnt);
  REPORT(L_ij);
  REPORT(k_ij);
  REPORT(sig_A);
  REPORT(sig_A_mr);
  REPORT(re_ai);
  REPORT(nll1);
  REPORT(nll2);
  REPORT(nll3);
  REPORT(nll4);
  REPORT(nll5);
  REPORT(nll7);
  REPORT(nll8);
  REPORT(nll9);
  REPORT(nll10);
  REPORT(nll11);
  REPORT(res);
  REPORT(beta1);
  REPORT(lnsig_bph);
  REPORT(predDmax); 
  REPORT(Lobs);
  REPORT(Lpred);
  REPORT(n_mr);
  REPORT(sig_L);
  REPORT(k_ij_mr);
  REPORT(re_Linfi_mr);
  REPORT(re_kij_mr);
  //REPORT(sig_kij_mr);
  REPORT(Linfi);
  REPORT(a_ij);
  REPORT(hum_obs_model);
  REPORT(mr_obs_model);
  REPORT(resFixed);
  REPORT(lnL0);
  REPORT(jcnt);
  REPORT(jcnt_mr);
  REPORT(sex_ki);
  REPORT(sex_Linfi);
  
  ADREPORT(Linf);
  ADREPORT(k);
  ADREPORT(Linf_i);
  ADREPORT(k_i);
  ADREPORT(Linf_i_mr);
  ADREPORT(k_i_mr);
  ADREPORT(L0);
  ADREPORT(beta1)
  ADREPORT(mu_a)
  ADREPORT(mu_a_mr)  
  ADREPORT(sig_x)
  ADREPORT(sig_A)
  ADREPORT(sig_A_mr)
  ADREPORT(sig_kij)
  //ADREPORT(sig_kij_mr)
  ADREPORT(sig_bph)
  ADREPORT(sig_ki)
  ADREPORT(sig_Linfi)  
  ADREPORT(sig_L)
  ADREPORT(a_i)
  ADREPORT(a_i_mr)
  return res;
    
}


