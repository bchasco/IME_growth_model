#http://www.sthda.com/english/wiki/ggplot2-barplot-easy-bar-graphs-in-r-software-using-ggplot2
#install.packages("devtools")
library(devtools)
install_github("easyGgplot2", "kassambara")

require(easyGgplot2)
require(kassambara)

tiff("ggplot_modelSelectionExperiment.tiff", width=620, height=620)

sim = read.table("simulationModelSelection_copy.out", header=TRUE)
sim = sim[,c(1:7,23)]
sim$bestMod = 0
sim$trueMod[sim$trueMod==1] = "vB"
sim$trueMod[sim$trueMod==2] = "gompertz"
sim$trueMod[sim$trueMod==3] = "logistic"
sim$simMod[sim$simMod==1] = "vB"
sim$simMod[sim$simMod==2] = "gompertz"
sim$simMod[sim$simMod==3] = "logistic"


mods = unique(sim$trueMod)
ss = unique(sim$sim)
n = unique(sim$sampleSize)
for(m in mods)
  for(s in ss)
  for(ni in n){
    tf = sim$TMBAIC[sim$trueMod==m & sim$sim==s & sim$sampleSize==ni]==min(sim$TMBAIC[sim$trueMod==m & sim$sim==s & sim$sampleSize==ni])
    sim$bestMod[sim$trueMod==m & sim$sim==s & sim$sampleSize==ni] = tf+1-1
  }


convTable2 = aggregate(sim$bestMod,list(sampleSize=sim$sampleSize, trueMod=sim$trueMod, simMod = sim$simMod), sum)
convTable2$trueMod[convTable2$trueMod==1] = "vB"
convTable2$trueMod[convTable2$trueMod==2] = "gompertz"
convTable2$trueMod[convTable2$trueMod==3] = "logistic"

convTable = aggregate(sim$conv,list(sampleSize=sim$sampleSize, trueMod=sim$trueMod, simMod = sim$simMod), sum)
convTable$trueMod[convTable$trueMod==1] = "vB"
convTable$trueMod[convTable$trueMod==2] = "gompertz"
convTable$trueMod[convTable$trueMod==3] = "logistic"
convTable$bestMod = rep(0,nrow(convTable))
for(i in 1:nrow(convTable)){
  convTable$bestMod[convTable$trueMod==convTable2$trueMod[i] & convTable$simMod==convTable2$simMod[i] & convTable$sampleSize==convTable2$sampleSize[i]] = convTable2$x[i]
}

sim = sim[sim$bestMod==1,]

convTable$x2 = paste(convTable$bestMod,"/", convTable$x)
myLabels = convTable
bp = ggplot2.barplot(data=sim, xName="simMod", 
                faceting=TRUE,
                xtitle="Fitted model",
                ytitle="Count",
                face="plain",
                ylim=c(0,104),
                facetingVarNames=c("trueMod","sampleSize"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  theme(axis.text.x=element_text(size=16, face="plain"), 
        axis.text.y=element_text(size=12, face="plain"),
        strip.text.x = element_text(size = 14, face="plain"), 
        strip.text.y = element_text(size = 16, face="plain"), 
        axis.title=element_text(size=14, face="plain")) +
  geom_text(data=convTable, aes(label=x, x = simMod, y = bestMod+2),
          hjust=0.5, vjust=0, size=4.5,
          inherit.aes = FALSE)

print(bp)


dev.off()

