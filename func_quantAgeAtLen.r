quantAgeAtLen = function(mu,cov,offset,nn){
  
  mu = sd[[map]][[model]]$value
  cov = sd[[map]][[model]]$cov
  offset=0
  nn = dataList$ni
  scnt = 1
  draw = matrix(0,nn,7)
  AatLsim = array(0,c(nn*100,length(1:100)))
  for(s in 1:nn){
    i1 = 2 + offset + s
    i2 = 2 + offset + nn + s
    i3 = which(names(mu)=="L0")
    indi = c(i1,i2,i3)
    mui = mu[indi]
    Vari = cov[indi,indi]
    for(ss in 1:100){ #sample each individual
      si = rmvnorm(1,mui,Vari)
      sLinfi = si[1]
      ski = si[2]
      sX0 = si[3]
      if(model==1){
        AatLsim[scnt,] = -log((1 - (1:100)/sLinfi)/((sLinfi-sX0)/sLinfi))/ski
      }
      if(model==2){
        tmpX0i = ((log(sLinfi) - log(sX0))/log(sLinfi))
      }
      if(model==3){
        sbi = 1-sX0/sLinfi
      }
      scnt = scnt + 1
    }
  }
  qAaL = apply(AatLsim,2,function(x){quantile(na.omit(x[x>0]),probs=c(0.025,0.5,0.975))})
  return(qAaL)
}
