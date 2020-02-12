#Things that affect performance
# set sig_ki = 0.1 at a minimum
# set X0 to 0.75 for the von Bertalanffy model.
rm(list=ls())

update.packages("Matrix")
library(TMB)
try(dyn.unload("NonIntegratedHumerus"))
compile("NonIntegratedHumerus.cpp")
dyn.load("NonIntegratedHumerus")

rep = rep(list(rep(list(NA),3)),4)
sd = rep(list(rep(list(NA),3)),4)
res = matrix(0,4,3)
fres = matrix(0,4,3)
maxgr = matrix(0,4,3)
conv = matrix(0,4,3)
mapList = list()

runModel = TRUE
for(sex in 0:0){ #0 = No sex, 1 = Sex.
  for(models in 1:3){
    for(maps in 1:4){
      if(maps==1) #full model -  - set this to two if you want 100% convergence
        minObs = 2
      if(maps==2) #no ki - set this to two if you want 100% convergence and estimates of age
        minObs = 2
      if(maps==3) #no linfi -  - set this to two if you want 100% convergence and estimates of age
        minObs = 2
      if(maps==4) #process only, no hierarchy -  - you can set this to 0 and still get 100% convergence
        minObs = 2
      source("readInRealData.r")
      
      dataList$growthModel = models
      parList = list(lnLinf = log(95),
                     lnk = log(0.1),
                     lnL0 = log(0.5),
                     lnmu_a = log(5),
                     lnsig_x = log(0.1),
                     lnsig_Linfi = log(2),
                     lnsig_kij = log(0.3),
                     lnsig_ki = log(0.2),
                     lnsig_A = log(0.5),
                     re_Linfi = rep(0,dataList$ni),
                     re_ki = rep(0,dataList$ni),
                     re_kij = rep(0,sum(dataList$hij)),
                     re_ai = rep(0,dataList$ni - dataList$nia),
                     lnbeta1 = log(0.33),
                     lnsig_bph = log(2),
                     sex_ki = rep(0,2),
                     sex_Linfi = rep(0,2)
      )
      
      mapList[[1]] = list(
      )
      
      mapList[[2]] = list(re_ki = as.factor(rep(NA,dataList$ni))
                          ,lnsig_ki = as.factor(NA)
      )
      
      mapList[[3]] = list(re_Linfi = as.factor(rep(NA,dataList$ni))
                          ,lnsig_Linfi = as.factor(NA)
      )
      
      mapList[[4]] = list(re_ki = as.factor(rep(NA,dataList$ni))
                          ,re_Linfi = as.factor(rep(NA,dataList$ni))
                          ,lnsig_Linfi = as.factor(NA)
                          ,lnsig_ki = as.factor(NA)
      )
      
      if(sex == 0){
        mapList[[maps]]$sex_ki = as.factor(rep(NA,2))
        mapList[[maps]]$sex_Linfi = as.factor(rep(NA,2))
      }
      
      if(runModel==TRUE){
        obj = MakeADFun(data=dataList, 
                        parameters=parList, 
                        map = mapList[[maps]],
                        DLL="NonIntegratedHumerus",
                        silent=FALSE,
                        #
                        random=c("re_Linfi","re_ki", "re_ai","re_kij"))
        
        Opt = tryCatch(TMBhelper::Optimize(obj,
                                           start=obj$par,
                                           objective=obj$fn,
                                           gradient=obj$gr,
                                           loopnum = 3,
                                           newtonsteps = 8,
                                           getsd=FALSE),
                       error=function(e) NULL)
        
        if(is.null(Opt)==FALSE){
          rep[[maps]][[models]] = obj$report()
          tmpsd = sdreport(obj)
          sd[[maps]][[models]] =  tmpsd 
          fres[maps,models] = Opt$AIC[1]
          conv[maps,models] = sd[[maps]][[models]]$pdHess
          maxgr[maps,models] = max(sd[[maps]][[models]]$gradient.fixed)
        }
      }
    }
  }
}#sex

bestModel = apply(res,1,function(x){which(x==min(x))})
save(rep, sd, dataList, fres,conv,maxgr,bestModel, file="Rep_NonIntegratedHumerus.Rdata")
