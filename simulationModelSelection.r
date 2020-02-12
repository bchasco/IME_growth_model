library(viridis)
library(mvtnorm)
library(TMB)
library(TMBhelper)
rm(list=ls())

try(dyn.unload("IntegratedHumerusMR"))
compile("IntegratedHumerusMR.cpp")
dyn.load("IntegratedHumerusMR")

try(dyn.unload("NonIntegratedMR"))
compile("NonIntegratedMR.cpp")
dyn.load("NonIntegratedMR")

try(dyn.unload("NonIntegratedHumerus"))
compile("NonIntegratedHumerus.cpp")
dyn.load("NonIntegratedHumerus")

source("func_simulation.r")

#set.seed(10)
#true pars
true_k = c(0.07,0.09,0.13)
true_Linf = 95 
true_x0 = 5
true_mu = c(8,15)  
true_B0 = 0.53  
true_B1 = 0.38  
true_sig = c(0.001,0.01)
true_siglink = 2   
true_phi_ki = 0.1  #low sigma is 0.001
true_phi_kij = 0.7  #low sigma is 0.007 
true_phi_Linf = 5  #low sigma is 0.003
true_phi_a = c(0.8,0.3)
true_shape = 0.45
true_scale = 2.2
true_rho = 0.5
true_phi_yr = 0.2
runOptim = TRUE
myPlot = FALSE
mySilent = TRUE
callSD = TRUE
fixSeed = FALSE
mySimCnt = 1
nSim = 1
#Run the simulations
ptm = proc.time()

nSim = 1
true_hum_obs_model = 1
true_mr_obs_model = 1

fixA = FALSE
true_nyears = c(1,1) #Leave set to 1 since there are no temporal effects.
tx = read.table("tagDat.dat", header=TRUE)
ttx = table(tx$McNeillID) - 1
sx = read.table("strandingData.dat", header=TRUE)
tsx = sx$nLAG - 1

myOut = paste("trueMod","sim","simMod","map","integrated","sampleSize","TMBAIC","X0","k","Linf","a","a_mr", "sig_A","sig_A_mr", "sig_ki", "sig_kij", "sig_Linfi", "sig_x", "sig_L", "sig_bph", "beta1", "gr","conv",sep="\t")
write.table(t(myOut), "simulationModelSelection.out", quote=FALSE, row.names = FALSE, col.names = FALSE)
mapList = list()
st = proc.time()

for(sampleSize in c(50,150,300)){
  true_ni_per_year = c(sampleSize,sampleSize)
  true_ni = true_ni_per_year * true_nyears
  true_nij = list(sample(tsx,true_ni[1],replace = TRUE),sample(ttx,true_ni[2],replace = TRUE)) #number of observations after the initial observation 
  for(trueModel in 1:3){
    for(sim in 1:100){
      x = func_simulation(arg_model=trueModel)
      
      for(modi in 1:3){ #model 1=VB, 2=gomp, 3=log
        for(inti in 1:1){ #1 = Hum/MR, 2 = Hum, 3 = MR: Integrated or non-integrated data sets
          for(mapi in 1:1){ #which random effects do you want
            if(inti==1){ #Fully integrated model
              
              mapList[[1]] = list()
              
              mapList[[2]] = list(re_ki_mr= as.factor(rep(NA,x$dataList$ni_mr))
                                  ,re_ki = as.factor(rep(NA,x$dataList$ni))
                                  ,lnsig_ki = as.factor(NA)
              )
              
              mapList[[3]] = list(re_Linfi_mr= as.factor(rep(NA,x$dataList$ni_mr))
                                  ,re_Linfi = as.factor(rep(NA,x$dataList$ni))
                                  ,lnsig_Linfi = as.factor(NA)
              )
              
              mapList[[4]] = list(re_ki_mr= as.factor(rep(NA,x$dataList$ni_mr))
                                  ,re_Linfi_mr= as.factor(rep(NA,x$dataList$ni_mr))
                                  ,re_ki = as.factor(rep(NA,x$dataList$ni))
                                  ,re_Linfi = as.factor(rep(NA,x$dataList$ni))
                                  ,lnsig_Linfi = as.factor(NA)
                                  ,lnsig_ki = as.factor(NA)
              )
              parList = list(lnLinf = log(95),
                             lnk = log(0.1),
                             lnL0 = log(0.5),
                             lnmu_a = log(5),
                             lnmu_a_mr = log(12),
                             lnsig_x = log(0.01),
                             lnsig_Linfi = log(2),
                             lnsig_kij = log(0.3),
                             lnsig_ki = log(0.2),
                             lnsig_A = log(0.5),
                             lnsig_A_mr = log(0.5),
                             re_Linfi = rep(0,x$dataList$ni),
                             re_ki = rep(0,x$dataList$ni),
                             re_kij = rep(0,sum(x$dataList$hij)),
                             re_ai = rep(0,x$dataList$ni - x$dataList$nia),
                             lnbeta1 = log(0.33),
                             lnsig_bph = log(2),
                             lnsig_L = log(0.5),
                             re_Linfi_mr = rep(0,x$dataList$ni_mr),
                             re_ki_mr = rep(0,x$dataList$ni_mr),
                             re_kij_mr = rep(0,sum(x$dataList$tij)),
                             re_ai_mr = rep(0,x$dataList$ni_mr)
              )
              
              x$dataList$growthModel = modi
              obj = MakeADFun(data=x$dataList,
                              parameters=x$parList,
                              map = mapList[[mapi]],
                              DLL="IntegratedHumerusMR",
                              silent=TRUE,
                              random=c("re_Linfi",
                                       "re_Linfi_mr",
                                       "re_ki",
                                       "re_kij",
                                       "re_ki_mr",
                                       "re_kij_mr",
                                       "re_ai",
                                       "re_ai_mr"))
              simOpt = tryCatch(TMBhelper::Optimize(obj, 
                                                    start=obj$par, 
                                                    objective=obj$fn, 
                                                    gradient=obj$gr,
                                                    loopnum = 3,
                                                    newtonsteps = 5,
                                                    getsd=FALSE), error=function(e) NULL)
              #print(paste(trueModel,modi,sim, simOpt$max_gradient, simOpt$max_gradient<0.0001))
              rep = obj$report()
            }#end if init ==1
            if(is.null(simOpt))
              print(paste(trueModel,modi,sim,"err"))
            if(is.null(simOpt)==FALSE){
              if(inti==1) myOut = c(trueModel,sim,modi,mapi,inti,sampleSize,TMBAIC(simOpt),rep$X0,rep$k,rep$Linf,rep$mu_a,rep$mu_a_mr,rep$sig_A,rep$sig_A_mr,rep$sig_ki, rep$sig_kij, rep$sig_Linfi, rep$sig_x, rep$sig_L, exp(rep$lnsig_bph), rep$beta1, max(obj$gr()),simOpt$max_gradient<0.0001)
              if(inti==2) myOut = c(trueModel,sim,modi,mapi,inti,sampleSize,TMBAIC(simOpt),rep$X0,rep$k,rep$Linf,rep$mu_a,NA,rep$sig_A,NA,rep$sig_ki, rep$sig_kij, rep$sig_Linfi, rep$sig_x, NA, NA, NA, max(obj$gr()),simOpt$max_gradient<0.0001)
              if(inti==3) myOut = c(trueModel,sim,modi,mapi,inti,sampleSize,TMBAIC(simOpt),rep$X0,rep$k,rep$Linf,NA,rep$mu_a_mr,NA,rep$sig_A_mr,rep$sig_ki, rep$sig_kij, rep$sig_Linfi, NA, rep$sig_L, NA, NA, max(obj$gr()),simOpt$max_gradient<0.0001)
              print(paste(trueModel,modi,sampleSize,sim,simOpt$max_gradient<0.0001))
              write.table(t(myOut), "simulationModelSelection.out", quote=FALSE, append = TRUE, row.names = FALSE, col.names = FALSE)
            }
            
          }#end mapi
        }#end inti
      }#end modi
    }#end simi
  }#end trueModel
}#end sample size

tt = proc.time() - st
print(tt)
