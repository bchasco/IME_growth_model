#https://plot.ly/ggplot2/box-plots/
#http://www.cookbook-r.com/Graphs/Lines_(ggplot2)/ for the horizontal lines
require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)
require(tidyr)

tiff("ggplot_persistentAndTransientVar_Sigma.tiff", width=620, height=820)
sim = read.table("simulationPersistentAndTransientVar_copy2.out", header=TRUE)
sim = sim[sim$varSim<=3,]
sim = sim[sim$k<=0.2,]
sim$conv = sim$gr * 0
sim$conv[sim$gr<0.0001] = 1

simSigma = sim[,c(4,6,13:20)]
simSigma = simSigma[simSigma$sig_L<=0.03,]
simSigma = simSigma[simSigma$sig_A_mr<=1,]

simSigma_long <- simSigma %>% 
  gather(key   = map,
         value = varSim)
names(simSigma_long) = c("mapi", "varSim", "var", "val")
#change the map labels
simSigma_long$mapi = as.factor(simSigma_long$map)
simSigma_long$varSim = as.factor(simSigma_long$varSim)

#Remove the estimates for the simulation where the parameters were not estimated because of mapping
for(v in 1:3){
  simSigma_long$val[simSigma_long$varSim==v & simSigma_long$mapi==2 & simSigma_long$var=="sig_ki"] = NA
  simSigma_long$val[simSigma_long$varSim==v & simSigma_long$mapi==3 & simSigma_long$var=="sig_Linfi"] = NA
  simSigma_long$val[simSigma_long$varSim==v & simSigma_long$mapi==4 & simSigma_long$var=="sig_ki"] = NA
  simSigma_long$val[simSigma_long$varSim==v & simSigma_long$mapi==4 & simSigma_long$var=="sig_Linfi"] = NA
}


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
simSigma_long$trueVal[simSigma_long$var=="sig_ki" & simSigma_long$varSim==3] = 0.7
simSigma_long$trueVal[simSigma_long$var=="sig_ki" & simSigma_long$varSim==2] = 0.4
simSigma_long$trueVal[simSigma_long$var=="sig_kij" & simSigma_long$varSim==1] = 0.7
simSigma_long$trueVal[simSigma_long$var=="sig_kij" & simSigma_long$varSim==3] = 0.1
simSigma_long$trueVal[simSigma_long$var=="sig_kij" & simSigma_long$varSim==2] = 0.4
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
levels(simSigma_long$mapi) = c(expression(""[1]~"both"),
                               expression(""[2]~epsilon[L[infinity,i]]),
                               expression(""[3]~epsilon[k[i]]),
                               expression(""[4]~"none"))

#Convergence and best AIC text at the bottom where the sigma[X] panel is  
convTable = aggregate(sim$conv,list(varSim=sim$varSim, mapi=sim$map), sum)
convTable$var = rep('sigma[X]',nrow(convTable))
bestAIC = aggregate(sim$bestAIC, list(varSim=sim$varSim, mapi=sim$map), sum)
bestAIC$var = rep('sigma[X]',nrow(bestAIC))
convTable$x = paste(convTable$x,"\n/",bestAIC$x, sep="")
convTable$bestAIC = bestAIC$x
convTable$mapi[convTable$mapi==1] = levels(simSigma_long$mapi)[1]
convTable$mapi[convTable$mapi==2] = levels(simSigma_long$mapi)[2]
convTable$mapi[convTable$mapi==3] = levels(simSigma_long$mapi)[3]
convTable$mapi[convTable$mapi==4] = levels(simSigma_long$mapi)[4]

#Labels for each box
myLabs = aggregate(list(val=simSigma_long$val),
                   by=list(var=simSigma_long$var,
                           varSim=simSigma_long$varSim,
                           mapi=simSigma_long$mapi), max)

myLabs = myLabs[myLabs$varSim==levels(simSigma_long$varSim)[1],]
myLabs = myLabs[order(myLabs$var),]
myMax = aggregate(list(val=(na.omit(simSigma_long$val))),by=list(var=simSigma_long$var[is.na(simSigma_long$val)==FALSE]), max)
myMin = aggregate(list(val=(na.omit(simSigma_long$val))),by=list(var=simSigma_long$var[is.na(simSigma_long$val)==FALSE]), min)
for(i in unique(myLabs$var))
   myLabs$val[myLabs$var==i] = myMax$val[myMax$var==i] - 0.1*(myMax$val[myMax$var==i] - myMin$val[myMin$var==i])
myLabs$labels = c(LETTERS,letters)[1:nrow(myLabs)]


#Build the dumb ass plot
bp = ggplot(simSigma_long, aes(x = varSim, y = val, fill=as.factor(varSim))) +
  geom_violin(size = .75) +
  facet_grid(var ~ mapi, scales="free_y", labeller=label_parsed) +
  xlab("\n\nRatio of persistent to transient variability") +
  ylab("Estimates based on simulated data \n") +
  labs(title="Random effects in estimation model") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues") +
  theme(legend.position="none") +
  geom_text(data=convTable, aes(label=bestAIC, x = varSim, y=-0.006),
            hjust=0.5, vjust=0,
            inherit.aes = FALSE)


# # Draw with separate lines for each bar
bp = bp +
  #TRUE parameters
  geom_errorbar(aes(ymax=trueVal, ymin=trueVal), colour="#AA0000") +
  #Labels for the x-axis
  # scale_x_discrete(labels = c('1' = expression(paste(phi[k[i]]  ==  0.1, " \n, ", phi[k[ij]]  ==  0.7)),
  #                             '2' = expression(paste(phi[k[i]]  ==  0.4, " \n, ", phi[k[ij]]  ==  0.4)),
  #                             '3' = expression(paste(phi[k[i]]  ==  0.7, " \n, ", phi[k[ij]]  ==  0.1)))
  # ) +
  scale_x_discrete(labels = c('1' = "lower",
                              '2' = "equal",
                              '3' = "higher"))+
  theme(axis.text.x=element_text(size=12),
        
        axis.text.y=element_text(size=12),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        axis.title=element_text(size=14))+ 
  #Labels for the boxes
  geom_text(data=myLabs, aes(label=labels, y = val),
            x = -Inf,  hjust=-2, vjust=0,
            inherit.aes = FALSE)

print(bp)

dev.off()
