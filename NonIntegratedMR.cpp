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
  Type g(Type x, Type par_b0, Type par_b1)
  {
    return par_b0 + par_b1 * x;}
  
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
  PARAMETER(lnmu_a_mr);//average age at stranding
  PARAMETER(lnsig_Linfi);//standard deviation of Dinf random effects
  PARAMETER(lnsig_kij);//standard deviation of transient k-effects
  PARAMETER(lnsig_ki);//strandard deviation of individual k-effects
  PARAMETER(lnsig_A_mr); //stanard deviation of age random effect

  //mark_recapture time-series simulation
  PARAMETER(lnsig_L);
  PARAMETER_VECTOR(re_Linfi_mr);
  PARAMETER_VECTOR(re_ki_mr);
  PARAMETER_VECTOR(re_kij_mr);
  PARAMETER_VECTOR(re_ai_mr);

  Type res = 0.;
  
  Type Linf = exp(lnLinf);   //Average max humerus
  Type L0 = exp(lnL0);
  Type k = exp(lnk);
  Type mu_a_mr = exp(lnmu_a_mr);
  Type sig_A_mr = exp(lnsig_A_mr);
  
  Type sig_L = exp(lnsig_L);
  Type sig_kij = exp(lnsig_kij);
  Type sig_Linfi = exp(lnsig_Linfi);
  Type sig_ki = exp(lnsig_ki);
  
  Type X0 = L0;
  int nij_mr = Lobs.size();
  vector<Type> Lpred(nij_mr);
  vector<Type> k_ij_mr(nij_mr);
  vector<Type> Linfi(ni_mr);
  vector<Type> a_ij_mr(nij_mr);
  
  int icnt = 0;
  int jcnt = 0;
  int jcnt_mr = 0;
  Type nll8 = 0.;
  vector<Type> a_i_mr(ni_mr);
  for(int ij=0; ij<nij_mr; ij++)
  {
    Linfi(Lmr_iID(ij)) = Linf + (re_Linfi_mr(Lmr_iID(ij)));      
    if(Lmr_ijID(ij)==0){
      k_ij_mr(ij) = k*exp(re_ki_mr(Lmr_iID(ij)));
      a_ij_mr(ij) = exp(lnmu_a_mr + re_ai_mr(Lmr_iID(ij)) -0.5*sig_A_mr*sig_A_mr);
      a_i_mr(icnt) = exp(lnmu_a_mr + re_ai_mr(Lmr_iID(ij)) -0.5*sig_A_mr*sig_A_mr);
      icnt++;
      if(growthModel==1)
        Lpred(ij) = xt(Linfi(Lmr_iID(ij)), L0, a_ij_mr(ij), k_ij_mr(ij));
      if(growthModel==2)
        Lpred(ij) = gomp_xt(Linfi(Lmr_iID(ij)), L0, a_ij_mr(ij), k_ij_mr(ij));
      if(growthModel==3)
        Lpred(ij) = logistic_xt(Linfi(Lmr_iID(ij)), L0, a_ij_mr(ij), k_ij_mr(ij));
    }
      
    if(Lmr_ijID(ij)>0){
      if(tij(ij)!=1){
        k_ij_mr(ij) = k*exp(re_ki_mr(Lmr_iID(ij)));
      }
      if(tij(ij)==1){
        k_ij_mr(ij) = k*exp(re_ki_mr(Lmr_iID(ij))+re_kij_mr(jcnt_mr));
        jcnt_mr++;
      }
      a_ij_mr(ij) = a_ij_mr(ij-1) + del_mrT(ij);
      if(growthModel==1)
        Lpred(ij) = dxdt(Lpred(ij-1), Linfi(Lmr_iID(ij)), del_mrT(ij), k_ij_mr(ij));
      if(growthModel==2)
        Lpred(ij) = gomp_dxdt(Lpred(ij-1), Linfi(Lmr_iID(ij)), del_mrT(ij), k_ij_mr(ij));
      if(growthModel==3)
        Lpred(ij) = logistic_dxdt(Lpred(ij-1), Linfi(Lmr_iID(ij)), del_mrT(ij), k_ij_mr(ij));
    }
    
    if(mr_obs_model==1)
     nll8 -= dnorm(log(Lobs(ij)),log(Lpred(ij))-0.5*sig_L*sig_L, sig_L, true);
  }
  


  //Mark recapture data streams and random effects
  Type nll9 = 0.;
  Type nll10 = 0.;
  Type nll11 = 0.;
  Type nll12 = 0.;
  if(mr_obs_model==1){
    nll9 += -1.*sum(dnorm(re_ki_mr,Type(-0.5*sig_ki*sig_ki),sig_ki,true));
    nll10 += -1.*sum(dnorm(re_kij_mr,Type(-0.5*sig_kij*sig_kij),sig_kij,true));
    nll11 += -1.*sum(dnorm(re_Linfi_mr,Type(0.),(sig_Linfi),true));
    nll12 += -1.*sum(dnorm(re_ai_mr,0,(sig_A_mr),true));
  }
  
  res = 
    nll8 +
    nll9 +
    nll10 +
    nll11 +
    nll12;
  vector<Type> Linf_i_mr = Linf  + (re_Linfi_mr);
  vector<Type> k_i_mr = k * exp(re_ki_mr);

  REPORT(Linf);
  REPORT(re_ki_mr);  
  REPORT(Linf_i_mr);
  REPORT(k);
  REPORT(k_i_mr);
  REPORT(X0);
  REPORT(L0);
  REPORT(a_ij_mr);
  REPORT(lnmu_a_mr);
  REPORT(mu_a_mr);
  REPORT(re_ai_mr);
  REPORT(a_i_mr);
  REPORT(sig_Linfi);
  REPORT(sig_ki);
  REPORT(sig_kij);
  REPORT(sig_A_mr);
  REPORT(nll8);
  REPORT(nll9);
  REPORT(nll10);
  REPORT(nll11);
  REPORT(nll12);
  REPORT(res);
  REPORT(Lpred);
  REPORT(n_mr);
  REPORT(sig_L);
  REPORT(k_ij_mr);
  REPORT(re_Linfi_mr);
  REPORT(re_kij_mr);
  REPORT(jcnt_mr);

  ADREPORT(Linf);
  ADREPORT(k);
  ADREPORT(Linf_i_mr);
  ADREPORT(k_i_mr);
  ADREPORT(L0);
    ADREPORT(sig_A_mr)
  return res;
    
}


