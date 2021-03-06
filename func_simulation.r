
x_func = function(par_k, par_x0, par_Linf, par_a, arg_model){
  if(arg_model==1)
    x = par_Linf *(1-((par_Linf-par_x0)/par_Linf)*exp(-par_k*par_a))
  if(arg_model==2)
    x = exp(log(par_Linf) *(1-((log(par_Linf)-log(par_x0))/log(par_Linf))*exp(-par_k*par_a)))
  if(arg_model==3)
  {
    tmpb = 1 - par_x0/par_Linf
    x = par_Linf*(1-tmpb)*exp(par_k*par_a)/((tmpb+(1-tmpb)*exp(par_k*par_a)))
  }
  return(x)
}
dx_func = function(par_k, par_Linf, arg_x, arg_delT, arg_model){
  if(arg_model==1)
    x = arg_x + (par_Linf - arg_x)*(1-exp(-par_k*arg_delT))
  if(arg_model==2)
    x = exp(log(arg_x) + (log(par_Linf) - log(arg_x))*(1-exp(-par_k*arg_delT)))
  if(arg_model==3)
    x = (1/(1/arg_x + (1/par_Linf - 1/arg_x)*(1-exp(-1.*par_k*arg_delT))))
  return(x)
}
link_func = function(par_b1, x){
  y = par_b1 * x
  return(y)
}

func_simulation = function(
  sim_k = true_k,
  sim_Linf = true_Linf,
  sim_x0 = true_x0,
  sim_rho = true_rho,
  sim_mu = true_mu, 
  sim_B1 = true_B1,
  sim_sig = true_sig,
  sim_siglink = true_siglink,  
  sim_phi_ki = true_phi_ki,
  sim_phi_kij = true_phi_kij,   
  sim_phi_Linf = true_phi_Linf,  
  sim_phi_a = true_phi_a,
  sim_shape = true_shape,
  sim_scale = true_scale,
  sim_ni_per_year = true_ni_per_year,
  sim_nyears = true_nyears,
  sim_nij = true_nij,
  sim_fixA = fixA,
  arg_model = trueModel
)
{
  ni = c(sim_ni_per_year[1]*sim_nyears[1],sim_ni_per_year[2]*sim_nyears[2]) #number of strandings
  #Create array for holding information
  x = list(h = NA, mr = NA)
  obsx = list(h = NA, mr = NA)
  a0 = NA
  xmax = NA
  obsxmax = NA
  SCL = NA
  Linfi = list(h = NA, mr = NA)
  ki = list(h = NA, mr = NA)
  kij = list(h = NA, mr = NA)
  ai = list(h = NA, mr = NA)
  aij = list(h = NA, mr = NA)
  re_ai = list(h = NA, mr = NA)
  dT = list(h = NA, mr = NA)
  xi = list(h = NA, mr = NA)
  xij = list(h = NA, mr = NA)
  yr_ef = NA
  xi_yr = rep(1:sim_nyears[1],each=sim_ni_per_year[1]) 
  xij_yr = NA

  #loop over the two data types
  for(d in 1:2){
    jcnt = 1
    ii = 1
    for(i in 1:ni[d]){
      #True age of first LAG
      if(d==1)
      {
        a0[i] = 0
        if(runif(1)<0.05){
          ai[[d]][i] = 0.75
          a0[i] = 1
        }
        if(a0[i]==0){
          re_ai[[d]][ii] = rnorm(1,0,sim_phi_a[d])
          ai[[d]][i] = sim_mu[d] * exp(re_ai[[d]][ii]-0.5*(sim_phi_a[d])^2)
          ii = ii + 1
        }
      }
      
      #Mark recapture age at tagging
      if(d==2)
      {
        re_ai[[d]][i] = rnorm(1,0,sim_phi_a[d])
        ai[[d]][i] = sim_mu[d] * exp(re_ai[[d]][i]-0.5*(sim_phi_a[d])^2)
      }
      
      aij[[d]][jcnt] = ai[[d]][i]
      
      #Growth parameters
      Linfi[[d]][i] = sim_Linf + rnorm(1,0,sim_phi_Linf)
      ki[[d]][i] = sim_k[trueModel] * exp(rnorm(1,-0.5*(sim_phi_ki)^2,sim_phi_ki))
      #First kij is just the ki
      kij[[d]][jcnt] = ki[[d]][i]
      
      #Index IDs
      xi[[d]][jcnt] = i
      xij[[d]][jcnt] = 1
      
      #Time-at-liberty
      dT[[d]][jcnt] = 0
      
      #Initial size-at-age
      x[[d]][jcnt] = x_func(par_k = ki[[d]][i],
                            par_x0 = sim_x0,
                            par_Linf = Linfi[[d]][i],
                            par_a = ai[[d]][i],
                            arg_model = arg_model)

      jcnt = jcnt + 1
      #"recaptures"
      for(j in 1:sim_nij[[d]][i]){
        #Index IDs
        xi[[d]][jcnt] = i
        xij[[d]][jcnt] = j+1
        
        #time at liberty
        if(d==1)
          dT[[d]][jcnt] = 1
        if(d==2)
          dT[[d]][jcnt] = rgamma(1,shape=sim_shape, scale=sim_scale)
        
        #age at nex time step
        aij[[d]][jcnt] = aij[[d]][jcnt-1] + dT[[d]][jcnt]
        

        #combination of persistent and transient k
        if(d==1)
          kij[[d]][jcnt] = ki[[d]][i] * exp(rnorm(1,-0.5*(sim_phi_kij)^2,sim_phi_kij))
        if(d==2)
          kij[[d]][jcnt] = ki[[d]][i] * exp(rnorm(1,-0.5*(sim_phi_kij)^2,sim_phi_kij))
        
        #subsequent sizes at age  
        x[[d]][jcnt] = dx_func(par_k = kij[[d]][jcnt],
                               par_Linf = Linfi[[d]][i],
                               arg_x = x[[d]][jcnt-1],
                               arg_delT = dT[[d]][jcnt],
                               arg_model = arg_model)
        jcnt = jcnt + 1
      }
      
      #Create link dataset
      if(d == 1){
        SCL[i] = dx_func(par_k = sim_k[trueModel],
                         par_Linf = Linfi[[d]][i],
                         arg_x = x[[d]][jcnt-1],
                         arg_delT = runif(1), #Add some random amount of additional growth to the last LAG
                         arg_model = arg_model)
        xmax[i] = link_func(sim_B1, SCL[i])
        obsxmax[i] = xmax[i] + rnorm(1,0,sim_siglink)
      }
    }
    #Body proportional hypothesis with observation error
     if(d==1){
        obsx[[d]] = link_func(sim_B1,x[[d]]) * exp(rnorm(length(x[[d]]),-0.5*sim_sig[d]*sim_sig[d],sim_sig[d]))
     }
     #Mark-recapture with observation error
     if(d==2){
        obsx[[d]] = x[[d]] * exp(rnorm(length(x[[d]]),-0.5*sim_sig[d]*sim_sig[d],sim_sig[d]))
     }
  }
  
  hij = rep(1,length(obsx[[1]]))
  hij[(xij[[1]] - 1)==0] = 0

  tij = rep(1,length(obsx[[2]]))
  tij[(xij[[2]] - 1)==0] = 0
  
  dataList = list(
    ni = ni[1],
    nia = sum(a0),
    a0 = a0,
    growthModel = trueModel, #1 is VB, #2 is gompertz
    rho_mod = 1, #Deprecated
    nLAG = sim_nij[[1]],
    Xobs = obsx[[1]],
    Xmax = obsxmax,
    SCLmax = SCL,
    del_xT = dT[[1]],
    x_iID = xi[[1]] - 1,
    x_ijID = xij[[1]] - 1,
    x_yrID = xij_yr - 1,
    x_stateID = 1, #Deprecated
    hij = hij,
    hum_obs_model = 1, #1 = log-normal, 2 = normal with cv * mean, 3 = normal with sigma
    mr_obs_model = 1, #1 = log-normal, 2 = normal with cv * mean, 3 = normal with sigma
    sex_ijID = rep(2,sum(sim_nij[[1]]+1)),
    Lobs = obsx[[2]],
    tij = tij,
    n_mr = sim_nij[[2]]+1,
    del_mrT = dT[[2]],
    Lmr_iID = xi[[2]]-1,
    Lmr_ijID = xij[[2]]-1,
    Lvec = seq(75,95,1),
    ni_mr = ni[2]
  )

  parList = list(lnLinf = log(95),
                 lnk = log(0.1),
                 lnL0 = log(5),
                 lnmu_a = log(sim_mu[1]),
                 lnmu_a_mr = log(sim_mu[2]),
                 lnsig_x = log(0.1),
                 lnsig_Linfi = log(2),
                 lnsig_kij = log(0.3),
                 lnsig_ki = log(0.2),
                 lnsig_A = log(0.2),
                 lnsig_A_mr = log(0.2),
                 re_Linfi = rep(0,ni[1]),
                 re_ki = rep(0,ni[1]),
                 re_kij = rep(0,sum(hij)),
                 re_ai = rep(0,dataList$nia),
                 lnbeta1 = log(sim_B1),
                 lnsig_bph = log(sim_siglink),
                 lnsig_L = log(0.5),
                 re_Linfi_mr = rep(0,dataList$ni_mr),
                 re_ki_mr = rep(0,dataList$ni_mr),
                 re_kij_mr = rep(0,length(tij)),
                 re_ai_mr = rep(0,dataList$ni_mr)
  )
  simList = list(aij= aij, 
                 x = x)
  return(output=list(simList=simList,dataList=dataList, parList=parList))
}
