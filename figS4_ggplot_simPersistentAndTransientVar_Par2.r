#https://plot.ly/ggplot2/box-plots/
#http://www.cookbook-r.com/Graphs/Lines_(ggplot2)/ for the horizontal lines
require(foreign)
require(ggplot2)
require(MASS)
require(Hmisc)
require(reshape2)
require(tidyr)

tiff("ggplot_persistentAndTransientVar_Par.tiff", width=620, height=820)
sim = read.table("simulationPersistentAndTransientVar_copy2.out", header=TRUE)
sim = sim[sim$varSim<=3,]
sim = sim[sim$k<=0.2,]
sim$conv = sim$gr * 0
sim$conv[sim$gr<0.0001] = 1

simPar = sim[,c(4,6,8:12,21)]
simPar = simPar[simPar$a_mr<=30,]

simPar_long <- simPar %>% 
  gather(key   = map,
         value = varSim)
names(simPar_long) = c("mapi", "varSim", "var", "val")
simPar_long$mapi = as.factor(simPar_long$mapi)
simPar_long$varSim = as.factor(simPar_long$varSim)

sim$bestAIC = 0
for(svar in unique(sim$varSim))
  for(ss in unique(sim$sim)){
    tf = (sim$TMBAIC[sim$varSim==svar & sim$sim==ss] == min(sim$TMBAIC[sim$varSim==svar & sim$sim==ss]))
    sim$bestAIC[sim$varSim==svar & sim$sim==ss] = tf
  }


simPar_long$trueVal = 1 #dummy place column
simPar_long$trueVal[simPar_long$var=="X0"] = 5
simPar_long$trueVal[simPar_long$var=="beta1"] = 0.38
simPar_long$trueVal[simPar_long$var=="a"] = 8
simPar_long$trueVal[simPar_long$var=="a_mr"] = 15
simPar_long$trueVal[simPar_long$var=="k"] = 0.07
simPar_long$trueVal[simPar_long$var=="Linf"] = 95
simPar_long$var[simPar_long$var=="X0"] = 'X[0]'
simPar_long$var[simPar_long$var=="beta1"] = 'beta[bph]'
simPar_long$var[simPar_long$var=="a"] = 'mu[a]^hum'
simPar_long$var[simPar_long$var=="a_mr"] = 'mu[a]^mr'
simPar_long$var[simPar_long$var=="k"] = 'mu[k]'
simPar_long$var[simPar_long$var=="Linf"] = 'mu[L[infinity]]'

# levels(simPar_long$mapi) = c(expression(" "[1]~epsilon[k[i]]~","~epsilon[L[infinity,i]]~","~epsilon[k[i~","~j]]),
#                                expression(" "[2]~epsilon[L[infinity,i]]~","~epsilon[k[i~","~j]]),
#                                expression(" "[3]~epsilon[k[i]]~","~epsilon[k[i~","~j]]),
#                                expression(" "[4]~epsilon[k[i~","~j]]))

#
levels(simPar_long$mapi) = c(expression(""[1]~"both"),
                               expression(""[2]~epsilon[L[infinity,i]]),
                               expression(""[3]~epsilon[k[i]]),
                               expression(""[4]~"none"))

#Convergence and best AIC text at the bottom  
convTable = aggregate(sim$conv,list(varSim=sim$varSim, mapi=sim$map), sum)
convTable$var = rep('X[0]',nrow(convTable))
bestAIC = aggregate(sim$bestAIC, list(varSim=sim$varSim, mapi=sim$map), sum)
bestAIC$var = rep('X[0]',nrow(bestAIC))
convTable$x2 = bestAIC$x
convTable$mapi[convTable$mapi==1] = levels(simPar_long$mapi)[1]
convTable$mapi[convTable$mapi==2] = levels(simPar_long$mapi)[2]
convTable$mapi[convTable$mapi==3] = levels(simPar_long$mapi)[3]
convTable$mapi[convTable$mapi==4] = levels(simPar_long$mapi)[4]


#Labels for each box
myLabs = aggregate(list(val=simPar_long$val),
                   by=list(var=simPar_long$var, 
                           varSim=simPar_long$varSim, 
                           mapi=simPar_long$mapi), max)
myLabs = myLabs[myLabs$varSim==levels(simPar_long$varSim)[1],]
myLabs = myLabs[order(myLabs$var),]
myMax = aggregate(list(val=simPar_long$val),by=list(var=simPar_long$var), max)
myMin = aggregate(list(val=simPar_long$val),by=list(var=simPar_long$var), min)
for(i in unique(myLabs$var))
  myLabs$val[myLabs$var==i] = myMax$val[myMax$var==i] - 0.1*(myMax$val[myMax$var==i] - myMin$val[myMin$var==i])
myLabs$labels = LETTERS[1:nrow(myLabs)]


#Build the dumb ass plot
bp = ggplot(simPar_long, aes(x = varSim, y = val, fill=as.factor(varSim))) +
  geom_violin(size = .75) +
  facet_grid(var ~ mapi, scales="free_y", labeller=label_parsed) +
#  facet_grid(var ~ mapi, scales="free_y", labeller=bquote(c("1" = "x", "2" = "x", "3" = "x", "4" = "x"))) +
  xlab("\n Ratio of persistent to transient variability") +
  ylab("Estimates based on simulated data \n") +
  #scale_y_continuous(sec.axis = "test") + 
  labs(title="Random effects in estimation model") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_brewer(palette="Blues") +
  theme(legend.position="none") +
  geom_text(data=convTable, aes(label=x2, x = varSim, y=3.1),
            hjust=0.5, vjust=0,
            inherit.aes = FALSE)


# Draw with separate lines for each bar
 bp = bp + 
   #TRUE parameters
   geom_errorbar(aes(ymax=trueVal, ymin=trueVal), colour="#AA0000") +
   #Labels for the x-axis
   # scale_x_discrete(labels = c('1' = expression(paste(phi[k[i]]  ==  0.1, " , ", phi[k[ij]]  ==  0.7)),
   #                              '2' = expression(paste(phi[k[i]]  ==  0.4, " , ", phi[k[ij]]  ==  0.4)),
   #                              '3' = expression(paste(phi[k[i]]  ==  0.7, " , ", phi[k[ij]]  ==  0.1)))
   #  ) +
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
