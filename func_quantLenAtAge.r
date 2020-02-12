quantLenAtAge = function(mu,cov,offset,nn){
  scnt = 1
  draw = matrix(0,nn,7)
  LaAsim = array(0,c(nn,length(0:80)))
  for(s in 1:nn){
    i1 = 2 + offset + s
    i2 = 2 + offset + nn + s
    i3 = which(names(mu)=="L0")
    indi = c(i1,i2,i3)
    mui = mu[indi]
    Vari = cov[indi,indi]
    si = rmvnorm(1,mui,Vari)

    draw[scnt,] = c(s,si,mui)
    
    sLinfi = si[1]
    ski = si[2]
    sX0 = si[3]
    
    if(model==1){
      LaAsim[scnt,] = sLinfi*(1-((sLinfi-sX0)/sLinfi)*exp(-ski*0:80))
    }
    if(model==2){
      tmpX0i = ((log(sLinfi) - log(sX0))/log(sLinfi))
      LaAsim[scnt,] = exp(log(sLinfi)*(1-tmpX0i*exp(-ski*(0:80))))
    }
    if(model==3){
      sbi = 1-sX0/sLinfi
      LaAsim[scnt,] = (sLinfi*(1-sbi)*exp(ski*0:80))/(sbi+(1-sbi)*exp(ski*0:80))
    }
    scnt = scnt + 1
  }
  qLaA = apply(LaAsim,2,function(x){quantile(x,probs=c(0.05,0.5,0.95))})
  return(qLaA)
}
