minObs = 1
d = source("readInRealData.r")

library(ggplot2)
library(viridis)
x = viridis(5)

myFiles = c("Rep_IntegratedHumerusMR.Rdata","Rep_NonIntegratedHumerus.Rdata", "Rep_NonIntegratedMR.Rdata")

# pdf("fig_ggplot_growthProcesses4_b.pdf", width=6.20, height=6.20)
png("figure_1.png", width=620, height=620)

source("func_simulation.r") #Get the growth functions to make the code easier to read.
#Labels for the different model runs
RandomEfs = c("Eps_kij, Eps_ki, Eps_Linfi", "Eps_kij, Eps_Linfi","Eps_kij, Eps_ki","Eps_kij")
GrowthProcess = c("von Bertalanffy", "gompertz", "logistic")

library(mvtnorm)
source("func_quantLenAtAge.r")
df = t(rep(NA,5))
pred = t(rep(NA,5))


myPoints = c(243)

dp = t(rep(NA,5))
for(myData in 1:3){
  load(myFiles[myData]) #load the results from the different data integration experiments
  for(map in 1:4){
    for(model in 1:3){
      tp = rep[[map]][[model]]
      if(length(sd[[map]][[model]])>1){
        if(myData<=2)
          myq = quantLenAtAge(sd[[map]][[model]]$value,
                            sd[[map]][[model]]$cov,
                            offset=0,
                            dataList$ni)  
        if(myData==3)
          myq = quantLenAtAge(sd[[map]][[model]]$value,
                              sd[[map]][[model]]$cov,
                              offset=0,
                              dataList$ni_mr)  
        df = rbind(df,cbind(rep(myData,ncol(myq)*2),
                            rep(map,ncol(myq)*2),
                            rep(model,ncol(myq)*2),
                            c(0:80,80:0),
                            c(myq[1,],-1*sort(-myq[3,]))))
      }
      
    }
  }
}
df = as.data.frame(df[2:nrow(df),])
names(df) = c("data","map","model","age","len")
levels(df$data) = c("Integrated", "Skeletochronology", "Tagged")
df$data[df$data=="1"] = "Integrated"
df$data[df$data=="2"] = "Skeletochronology"
df$data[df$data=="3"] = "Tagged"
df$map[df$map=="1"] = "both"
df$map[df$map=="2"] = "L_infty"
df$map[df$map=="3"] = "k"
df$map[df$map=="4"] = "none"

# levels(df$map) = c(expression(""[1]~"both"),
#                    expression(""[2]~epsilon[L[infinity,i]]),
#                    expression(""[3]~epsilon[k[i]]),
#                    expression(""[4]~"none"))

df$map = factor(df$map)
df$model = factor(df$model)

library(viridis)
x = viridis(5)

myLabels = data.frame(rep(levels(df$data),4),
                      rep(levels(df$map),each=3),
                      rep(0,12),
                      rep(100,12),LETTERS[1:12])
names(myLabels) = c("data", "map", "age","len","lab")

p <- ggplot(df, aes(x=age, y=len, group=model)) +
  xlab("\nAge (years)") +
  xlim(0,80) + 
  geom_polygon(data = df,
               # alpha = 0.3,
               aes(x=age,
                   y = len, 
                   group=model, 
                   fill=x[as.numeric(df$model)]),
               colour="grey") +
  geom_text(data=myLabels, aes(label=lab, y = len, x = age),
            hjust=0, vjust=0,
            inherit.aes = FALSE) +
  theme(panel.spacing.x=unit(0.0, "lines"), 
        panel.spacing.y=unit(0,"lines")) +  
  theme(panel.border = element_rect(color="black", fill=NA)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) + 
  scale_fill_manual( values = viridis(3,alpha = 0.5),
                     labels=c("logistic","gompertz","von Bertlanffy"),
                     name="model") +
  ylab("Straight carapace length (cm)\n") +
  facet_grid(map ~ data)

p <- p + theme(axis.text.x=element_text(size=12),
               axis.text.y=element_text(size=12),
               strip.text.x = element_text(size = 12),
               strip.text.y = element_text(size = 12),
               strip.background = element_rect(colour="black",fill="transparent"),
               axis.title=element_text(size=12)) 
# 
p <- p + labs(fill = "Growth process")

print(p)

dev.off()
