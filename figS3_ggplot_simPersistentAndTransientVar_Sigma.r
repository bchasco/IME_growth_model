#https://plot.ly/ggplot2/box-plots/
#http://www.cookbook-r.com/Graphs/Lines_(ggplot2)/ for the horizontal lines
require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)
require(tidyr)

tiff("ggplot_persistentAndTransientVar_Sigma.tiff", width = 620, height=920)
sim = read.table("simulationPersistentAndTransientVar_copy2.out", header=TRUE)
sim$conv = sim$gr * 0
sim$conv[sim$gr<0.0001] = 1

simSigma = sim[,c(4,6,13:20)]
simSigma = simSigma[simSigma$varSim<=3,]
simSigma = simSigma[simSigma$sig_L<=0.03,]
simSigma = simSigma[simSigma$sig_A_mr<=1,]
#simSigma = simSigma[simSigma$sig_x<=0.003,]

simSigma_long <- simSigma %>% 
  gather(key   = map,
         value = varSim)
names(simSigma_long) = c("mapi", "varSim", "var", "val")

sim$bestAIC = 0
for(svar in unique(sim$varSim))
  for(ss in unique(sim$sim)){
    tf = (sim$TMBAIC[sim$varSim==svar & sim$sim==ss] == min(sim$TMBAIC[sim$varSim==svar & sim$sim==ss]))
    sim$bestAIC[sim$varSim==svar & sim$sim==ss] = tf
  }

#
simSigma_long$trueVal = 1 #dummy place column
simSigma_long$trueVal[simSigma_long$var=="sig_A"] = 0.8
simSigma_long$trueVal[simSigma_long$var=="sig_A_mr"] = 0.3
simSigma_long$trueVal[simSigma_long$var=="sig_bph"] = 2
simSigma_long$trueVal[simSigma_long$var=="sig_ki" & simSigma_long$varSim==1] = 0.1
simSigma_long$trueVal[simSigma_long$var=="sig_ki" & simSigma_long$varSim==3] = 0.7
simSigma_long$trueVal[simSigma_long$var=="sig_ki" & simSigma_long$varSim==2] = 0.4
simSigma_long$trueVal[simSigma_long$var=="sig_kij" & simSigma_long$varSim==1] = 0.7
simSigma_long$trueVal[simSigma_long$var=="sig_kij" & simSigma_long$varSim==3] = 0.1
simSigma_long$trueVal[simSigma_long$var=="sig_kij" & simSigma_long$varSim==2] = 0.4
simSigma_long$trueVal[simSigma_long$var=="sig_Linfi"] = 5
simSigma_long$trueVal[simSigma_long$var=="sig_L"] = 0.01
simSigma_long$trueVal[simSigma_long$var=="sig_x"] = 0.001

for(v in 1:3){
  simSigma_long$val[simSigma_long$varSim==v & simSigma_long$mapi==2 & simSigma_long$var=="sig_ki"] = NA
  simSigma_long$val[simSigma_long$varSim==v & simSigma_long$mapi==3 & simSigma_long$var=="sig_Linfi"] = NA
  simSigma_long$val[simSigma_long$varSim==v & simSigma_long$mapi==4 & simSigma_long$var=="sig_ki"] = NA
  simSigma_long$val[simSigma_long$varSim==v & simSigma_long$mapi==4 & simSigma_long$var=="sig_Linfi"] = NA
}

simSigma_long$err = (simSigma_long$trueVal-simSigma_long$val)/simSigma_long$trueVal*100 

varsimSigmas = matrix(c(c(0.1,0.7,0.4),c(0.7,0.1,0.4)),3,2)

simSigma_long$varSim[simSigma_long$varSim==1] = "ki0.1_kij0.7"
simSigma_long$varSim[simSigma_long$varSim==2] = "ki0.7_kij0.1"
simSigma_long$varSim[simSigma_long$varSim==3] = "ki0.4_kij0.4"

simSigma_long$varSim = as.factor(simSigma_long$varSim)

# simSigma_long$varSim[simSigma_long$varSim==4] = "ki0.7_kij0.1_Linf1"
simSigma_long$mapi[simSigma_long$mapi==1] = "all"
simSigma_long$mapi[simSigma_long$mapi==2] = "no_ki"
simSigma_long$mapi[simSigma_long$mapi==3] = "no_Linfi"
simSigma_long$mapi[simSigma_long$mapi==4] = "no_ki_Linfi"

simSigma_long$var[simSigma_long$var=="sig_A"] = "phi[A]^hum"
simSigma_long$var[simSigma_long$var=="sig_A_mr"] = "phi[A]^mr"
simSigma_long$var[simSigma_long$var=="sig_bph"] = "sigma[bph]"
simSigma_long$var[simSigma_long$var=="sig_ki"] = "phi[ki]"
simSigma_long$var[simSigma_long$var=="sig_kij"] = "phi[kij]"
simSigma_long$var[simSigma_long$var=="sig_Linfi"] = "phi[L[infinity,i]]"
simSigma_long$var[simSigma_long$var=="sig_L"] = "sigma[L]"
simSigma_long$var[simSigma_long$var=="sig_x"] = "sigma[X]"
#

levels(simSigma_long$varSim) <- c(expression(paste(phi[k[i]]  ==  0.1, " , ", phi[k[ij]]  ==  0.7)),
                                  expression(paste(phi[k[i]]  ==  0.4, " , ", phi[k[ij]]  ==  0.4)),
                                  expression(paste(phi[k[i]]  ==  0.7, " , ", phi[k[ij]]  ==  0.1)))

convTable = aggregate(sim$conv,list(varSim=sim$varSim, mapi=sim$map), sum)
convTable$var = rep('sigma[X]',nrow(convTable))
#
convTable$varSim[convTable$varSim==1] = levels(simSigma_long$varSim)[1]
convTable$varSim[convTable$varSim==2] = levels(simSigma_long$varSim)[2]
convTable$varSim[convTable$varSim==3] = levels(simSigma_long$varSim)[3]

bestAIC = aggregate(sim$bestAIC, list(varSim=sim$varSim, mapi=sim$map), sum)
bestAIC$var = rep('sigma[X]',nrow(bestAIC))
bestAIC$varSim[bestAIC$varSim==1] = levels(simSigma_long$varSim)[1]
bestAIC$varSim[bestAIC$varSim==2] = levels(simSigma_long$varSim)[2]
bestAIC$varSim[bestAIC$varSim==3] = levels(simSigma_long$varSim)[3]

convTable$x = paste(convTable$x,"/",bestAIC$x, sep="")

myLabs = aggregate(list(val=simSigma_long$val),by=list(var=simSigma_long$var, varSim=simSigma_long$varSim, mapi=simSigma_long$mapi), max)
myLabs = myLabs[myLabs$mapi=="no_ki_Linfi",]
myLabs = myLabs[order(myLabs$var),]
myMax = aggregate(list(val=simSigma_long$val),by=list(var=simSigma_long$var), max)
myMin = aggregate(list(val=simSigma_long$val),by=list(var=simSigma_long$var), min)
for(i in unique(myLabs$var))
  myLabs$val[myLabs$var==i] = myMax$val[myMax$var==i] - 0.1*(myMax$val[myMax$var==i] - myMin$val[myMin$var==i])
myLabs$labels = LETTERS[1:nrow(myLabs)]

bp = ggplot(simSigma_long, aes(x = mapi, y = val, fill=as.factor(mapi))) +
  geom_violin(size = .75) +
  facet_grid(var ~ varSim, scales="free_y", labeller=label_parsed) +
  xlab("\n\n Random effects included") +
  ylab("Estimates based on simulated data \n") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
  scale_fill_brewer(palette="Blues") +
  theme(legend.position="none")

# # Draw with separate lines for each bar
bp = bp + geom_errorbar(aes(ymax=trueVal, ymin=trueVal), colour="#AA0000") +
  scale_x_discrete(labels = c('all' = expression(paste(epsilon[k[i]],", ",epsilon[k[ij]],", ",epsilon[L[infinity]])),
                              'no_ki'   = expression(paste(epsilon[k[ij]],", ",epsilon[L[infinity]])),
                              'no_Linfi' = expression(paste(epsilon[k[i]],", ",epsilon[k[ij]])),
                              'no_ki_Linfi' = expression(epsilon[k[ij]]))) +
  theme(axis.text.x=element_text(size=18), 
        axis.text.y=element_text(size=12),
        strip.text.x = element_text(size = 14), 
        strip.text.y = element_text(size = 16), 
        axis.title=element_text(size=14)) +
  geom_text(data=myLabs, aes(label=labels, y = val), 
            x = Inf,  hjust=2, vjust=0,
            inherit.aes = FALSE)+
  geom_text(data=convTable, aes(label=x, x = mapi, y=-0.007),
            hjust=0.5, vjust=0.,
            inherit.aes = FALSE)
print(bp)

dev.off()
