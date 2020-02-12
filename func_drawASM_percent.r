drawASM_percent = function(mu,cov,offset,nn,mySample = 100, mu_rSSM = 88, sig_rSSM=2.5){

  bmean= 0.95
  
  bvar = (0.95*0.02)^2
  alpha = ((1-bmean)/bvar - 1/bmean)*bmean^2
  beta = alpha*(1/bmean -1)
  
  scnt = 1
  draw = matrix(0,nn,7)
  ASM = NA
  ASMmu = NA
  for(s in 1:nn){
    i1 = 2 + offset + s
    i2 = 2 + offset + nn + s
    i3 = which(names(mu)=="L0")
    indi = c(i1,i2,i3)
    indmu = c(1,2,i3) #deterministic individual
    mui = mu[indi]
    Vari = cov[indi,indi]
    mu_ = mu[indmu]
    Var_ = cov[indmu,indmu]
    smu = rmvnorm(mySample,mu_,Var_)
    si = rmvnorm(mySample,mui,Vari)
    
    Linf = si[,1]
    k = si[,2]
    X0 = si[,3]
    Linfmu = smu[,1]
    kmu = smu[,2]
    
    rSSM = rnorm(mySample,mu_rSSM,sig_rSSM)
    
    if(model==1){
      x0tmp = (Linf - X0)/Linf
      x0tmpmu = (Linfmu - X0)/Linfmu
    }
    if(model==2){
      x0tmp = (log(Linf) - log(X0))/log(Linf)
      x0tmpmu = (log(Linfmu) - log(X0))/log(Linfmu)
    }
    # if(model==3)
    #   x0tmp = tmpb = 1 - X0/Linf
    
    if(model==1){
      ASM = c(ASM,-1.*log((1-((Linf*rbeta(mySample,alpha,beta))/Linf)) / x0tmp)/k)
      ASMmu = c(ASMmu,-1.*log((1-((Linfmu*rbeta(mySample,alpha,beta))/Linfmu)) / x0tmpmu)/kmu)
    }
    if(model==2){
      ASM = c(ASM,-1.*log((1-(log(Linf*rbeta(mySample,alpha,beta))/log(Linf))) / x0tmp)/k)
      ASMmu = c(ASMmu,-1.*log((1-(log(Linfmu*rbeta(mySample,alpha,beta))/log(Linfmu))) / x0tmpmu)/kmu)
    }
  }
  
  return(list(ASMi=ASM,ASMmu=ASMmu))
}
