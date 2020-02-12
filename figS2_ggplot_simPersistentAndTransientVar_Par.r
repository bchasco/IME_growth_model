#https://plot.ly/ggplot2/box-plots/
#http://www.cookbook-r.com/Graphs/Lines_(ggplot2)/ for the horizontal lines
require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)
require(tidyr)

#tiff("ggplot_persistentAndTransientVar_Sigma2.tiff", width=620, height=920)
sim = read.table("simulationPersistentAndTransientVar_copy2.out", header=TRUE)
sim = sim[sim$varSim<=3,]
sim = sim[sim$k<=0.2,]
sim$conv = sim$gr * 0
sim$conv[sim$gr<0.0001] = 1
sim$map = as.factor(sim$map)
sim$varSim = as.factor(sim$varSim)

simSigma = sim[,c(4,6,8:12,21)]
simSigma = simSigma[simSigma$a_mr<=30,]

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


simSigma_long$trueVal = 1 #dummy place column
simSigma_long$trueVal[simSigma_long$var=="sig_A"] = 0.8
simSigma_long$trueVal[simSigma_long$var=="sig_A_mr"] = 0.3
simSigma_long$trueVal[simSigma_long$var=="sig_bph"] = 2
simSigma_long$trueVal[simSigma_long$var=="sig_ki" & simSigma_long$varSim==1] = 0.1
simSigma_long$trueVal[simSigma_long$var=="sig_ki" & simSigma_long$varSim==2] = 0.7
simSigma_long$trueVal[simSigma_long$var=="sig_ki" & simSigma_long$varSim==3] = 0.4
simSigma_long$trueVal[simSigma_long$var=="sig_kij" & simSigma_long$varSim==1] = 0.7
simSigma_long$trueVal[simSigma_long$var=="sig_kij" & simSigma_long$varSim==2] = 0.1
simSigma_long$trueVal[simSigma_long$var=="sig_kij" & simSigma_long$varSim==3] = 0.4
simSigma_long$trueVal[simSigma_long$var=="sig_Linfi"] = 5
simSigma_long$trueVal[simSigma_long$var=="sig_L"] = 0.01
simSigma_long$trueVal[simSigma_long$var=="sig_x"] = 0.001

simSigma_long$var[simSigma_long$var=="sig_A"] = "phi[A]^hum"
simSigma_long$var[simSigma_long$var=="sig_A_mr"] = "phi[A]^mr"
simSigma_long$var[simSigma_long$var=="sig_bph"] = "sigma[bph]"
simSigma_long$var[simSigma_long$var=="sig_ki"] = "phi[ki]"
simSigma_long$var[simSigma_long$var=="sig_kij"] = "phi[kij]"
simSigma_long$var[simSigma_long$var=="sig_Linfi"] = "phi[L[infinity,i]]"
simSigma_long$var[simSigma_long$var=="sig_L"] = "sigma[L]"
simSigma_long$var[simSigma_long$var=="sig_x"] = "sigma[X]"
#

#Convergence and best AIC text at the bottom  
convTable = aggregate(sim$conv,list(varSim=sim$varSim, mapi=sim$map), sum)
convTable$var = rep('X[0]',nrow(convTable))
bestAIC = aggregate(sim$bestAIC, list(varSim=sim$varSim, mapi=sim$map), sum)
bestAIC$var = rep('X[0]',nrow(bestAIC))
convTable$x = paste(convTable$x,"/",bestAIC$x, sep="")


#Labels for each box
myLabs = aggregate(list(val=simSigma_long$val),
                   by=list(var=simSigma_long$var, 
                           varSim=simSigma_long$varSim, 
                           mapi=simSigma_long$mapi), max)
myLabs = myLabs[myLabs$varSim==levels(simSigma_long$varSim)[3],]
myLabs = myLabs[order(myLabs$var),]
myMax = aggregate(list(val=simSigma_long$val),by=list(var=simSigma_long$var), max)
myMin = aggregate(list(val=simSigma_long$val),by=list(var=simSigma_long$var), min)
for(i in unique(myLabs$var))
  myLabs$val[myLabs$var==i] = myMax$val[myMax$var==i] - 0.1*(myMax$val[myMax$var==i] - myMin$val[myMin$var==i])
myLabs$labels = LETTERS[1:nrow(myLabs)]


#Build the dumb ass plot
bp = ggplot(simSigma_long, aes(x = varSim, y = val, fill=as.factor(varSim))) +
  geom_violin(size = .75) +
  facet_grid(var ~ mapi, scales="free_y", labeller=label_parsed) +
  #  facet_grid(var ~ mapi, scales="free_y", labeller=bquote(c("1" = "x", "2" = "x", "3" = "x", "4" = "x"))) +
  xlab("\n\n Persistent and transient variability") +
  ylab("Estimates based on simulated data \n") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
  scale_fill_brewer(palette="Blues") +
  theme(legend.position="none") +
  geom_text(data=convTable, aes(label=x, x = varSim, y=3.1),
            hjust=0.5, vjust=0,
            inherit.aes = FALSE)


# Draw with separate lines for each bar
bp = bp + 
  #TRUE parameters
  geom_errorbar(aes(ymax=trueVal, ymin=trueVal), colour="#AA0000") +
  #Labels for the x-axis
  scale_x_discrete(labels = c('1' = expression(paste(phi[k[i]]  ==  0.1, " , ", phi[k[ij]]  ==  0.7)),
                              '2' = expression(paste(phi[k[i]]  ==  0.4, " , ", phi[k[ij]]  ==  0.4)),
                              '3' = expression(paste(phi[k[i]]  ==  0.7, " , ", phi[k[ij]]  ==  0.1)))
  ) +
  theme(axis.text.x=element_text(size=18),
        axis.text.y=element_text(size=12),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 16),
        axis.title=element_text(size=14)) 
#Labels for the boxes
geom_text(data=myLabs, aes(label=labels, y = val),
          x = Inf,  hjust=2, vjust=0,
          inherit.aes = FALSE)
print(bp)

#dev.off()
