#Things that affect performance
# set sig_ki = 0.1 at a minimum
# set X0 to 0.75 for the von Bertalanffy model.
rm(list=ls())

update.packages("Matrix")
library(TMB)
library(TMBhelper)
try(dyn.unload("NonIntegratedMR"))
compile("NonIntegratedMR.cpp")
dyn.load("NonIntegratedMR")

rep = rep(list(rep(list(NA),3)),4)
sd = rep(list(rep(list(NA),3)),4)
res = matrix(0,4,3)
fres = matrix(0,4,3)
maxgr = matrix(0,4,3)
conv = matrix(0,4,3)
mapList = list()

runModel = TRUE

for(models in 1:2){
  for(maps in 1:4){
    if(maps==1) #full model -  - set this to two if you want 100% convergence
      minObs = 2 #This refers to whether or not k_ij is estimated for observations > minObs
    if(maps==2) #no ki - set this to two if you want 100% convergence and estimates of age
      minObs = 2
    if(maps==3) #no linfi -  - set this to two if you want 100% convergence and estimates of age
      minObs = 2
    if(maps==4) #process only, no hierarchy -  - you can set this to 0 and still get 100% convergence
      minObs = 2
    source("readInRealData.r")
    
    dataList$growthModel = models
    
    mapList[[1]] = list()
    
    mapList[[2]] = list(re_ki_mr = as.factor(rep(NA,dataList$ni_mr))
                        ,lnsig_ki = as.factor(NA)
    )
    
    mapList[[3]] = list(re_Linfi_mr = as.factor(rep(NA,dataList$ni_mr))
                        ,lnsig_Linfi = as.factor(NA)
    )
    
    mapList[[4]] = list(re_ki_mr = as.factor(rep(NA,dataList$ni_mr))
                        ,re_Linfi_mr = as.factor(rep(NA,dataList$ni_mr))
                        ,lnsig_Linfi = as.factor(NA)
                        ,lnsig_ki = as.factor(NA)
    )
    
    parList = list(lnLinf = log(95),
                   lnk = log(0.1),
                   lnL0 = log(0.5),
                   lnmu_a_mr = log(12),
                   lnsig_Linfi = log(2),
                   lnsig_kij = log(0.3),
                   lnsig_ki = log(0.2),
                   lnsig_A_mr = log(0.5),
                   lnsig_L = log(0.5),
                   re_Linfi_mr = rep(0,dataList$ni_mr),
                   re_ki_mr = rep(0,dataList$ni_mr),
                   re_kij_mr = rep(0,sum(dataList$tij)),
                   re_ai_mr = rep(0,dataList$ni_mr)
    )
    
    if(runModel==TRUE){
      obj = MakeADFun(data=dataList, 
                      parameters=parList, 
                      map = mapList[[maps]],
                      DLL="NonIntegratedMR",
                      silent=FALSE,
                      #
                      random=c("re_Linfi_mr", "re_ki_mr", "re_ai_mr","re_kij_mr"))
      
#      Opt = try(nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr))
      Opt = tryCatch(TMBhelper::Optimize(obj, 
                                         start=obj$par, 
                                         objective=obj$fn, 
                                         gradient=obj$gr,
                                         # loopnum = 3,
                                         newtonsteps = 8,
                                         getsd=FALSE), error=function(e) NULL)
      
      if(is.null(Opt)==FALSE){
        rep[[maps]][[models]] = obj$report()
        tmpsd = sdreport(obj)
        sd[[maps]][[models]] =  tmpsd #list(value=tmpsd$value,sd=tmpsd$sd,cov.fixed=tmpsd$cov.fixed)  
        fres[maps,models] = TMBAIC(Opt)
        conv[maps,models] = sd[[maps]][[models]]$pdHess
        maxgr[maps,models] = max(sd[[maps]][[models]]$gradient.fixed)
      }
    }
  }
}

#Opt = TMBhelper::Optimize(obj, newtonsteps = 3)

bestModel = apply(res,1,function(x){which(x==min(x))})

#save(rep, sd, dataList, fres,conv,maxgr,bestModel, file="Rep_NonIntegratedMR.Rdata")
