
minObs = 1 #This is only needed if you're not running the model from scratch.
d = source("readInRealData.r")
library(ggplot2)
png("fig_ggplot_age4.png", width=620, height=620)

myFiles = c("Rep_IntegratedHumerusMR.Rdata","Rep_NonIntegratedHumerus.Rdata", "Rep_NonIntegratedMR.Rdata")
source("func_simulation.r") #Get the growth functions to make the code easier to read.
#Labels for the different model runs
RandomEfs = c("Eps_kij, Eps_ki, Eps_Linfi", "Eps_kij, Eps_Linfi","Eps_kij, Eps_ki","Eps_kij")
GrowthProcess = c("von Bertalanffy", "gompertz", "logistic")
myModels = c("Integrated","Non-integrated","Non-integrated")
myAges = c("Age", "ASM")
myData = c("Skeletochronology", "Tagging")
library(mvtnorm)
source("func_quantLenAtAge.r")
source("func_drawASM.r")
source("func_drawASM_percent.r")


tm = aggregate(d$value$del_mrT, list(id = d$value$Lmr_iID), max)

df = data.frame(dataType=NA,modelType=NA,ageType=NA,age=NA)

for(mod in c(1,2,3)){ #aggregated, #stranding, #mark-recapture
  
  #change the data that was used for the different analayses
  load(myFiles[mod])
  loc = which(fres == min(fres), arr.ind = TRUE)
  #age
  age = c(rep[[loc[1]]][[loc[2]]]$a_i + d$value$nLAG - 1,
          rep[[loc[1]]][[loc[2]]]$a_i_mr + tm$x)
  
  print(paste(myModels[mod],"ai"))
  print(quantile(rep[[loc[1]]][[loc[2]]]$a_i + d$value$nLAG - 1,probs=c(0.025,0.5,0.975)))
  print(paste(myModels[mod],"ai_mr"))
  print(quantile(rep[[loc[1]]][[loc[2]]]$a_i_mr + tm$x, probs=c(0.025,0.5,0.975)))
  
  #age type
  dataType = factor(c(rep(myData[1],length(rep[[loc[1]]][[loc[2]]]$a_i)),
                     rep(myData[2],length(rep[[loc[1]]][[loc[2]]]$a_i_mr))))
  ageType = factor(rep(myAges[1],length(dataType))) #Age = 1, ASM = 2
  modelType = factor(rep(myModels[mod],length(dataType)))
  df_tmp = data.frame(dataType,modelType,ageType,age)
  df = rbind(df,df_tmp)
  
  
  for(ma in 2){#loc[1]){ #map
    for(mo in 1){#loc[2]){ #model
      model = mo
      if(mod==1){
        nn=dataList$ni
        ASMhum = drawASM_percent(sd[[ma]][[mo]]$value,
                         sd[[ma]][[mo]]$cov,
                         offset=0,
                         nn=nn,
                         mu_rSSM = 92,
                         sig_rSSM = 6.6)
        age = c(ASMhum$ASMi)
        dataType = factor(rep(myData[1],length(ASMhum$ASMi))) #hum vs MR
        ageType = factor(rep(myAges[2],length(dataType)))#age = 1, ASM = 2
        modelType = factor(rep(myModels[mod],length(dataType))) #Integrated, stranding or tagging
        df_tmp = data.frame(dataType,modelType,ageType,age)
        df = rbind(df,df_tmp)

        ASMhum = drawASM_percent(sd[[ma]][[mo]]$value,
                                 sd[[ma]][[mo]]$cov,
                                 offset=nn,
                                 nn=dataList$ni_mr,
                                 mu_rSSM = 92,
                                 sig_rSSM = 6.6)
        age = c(ASMhum$ASMi)
        dataType = factor(rep(myData[2],length(ASMhum$ASMi))) #hum vs MR
        ageType = factor(rep(myAges[2],length(dataType)))#age = 1, ASM = 2
        modelType = factor(rep(myModels[mod],length(dataType))) #Integrated, stranding or tagging
        df_tmp = data.frame(dataType,modelType,ageType,age)
        df = rbind(df,df_tmp)
        
      }
      if(mod==2){
        nn=dataList$ni
        ASMhum = drawASM_percent(sd[[ma]][[mo]]$value,
                                 sd[[ma]][[mo]]$cov,
                                 offset=0,
                                 nn=nn,
                                 mu_rSSM = 92,
                                 sig_rSSM = 6.6)
        age = c(ASMhum$ASMi)
        dataType = factor(rep(myData[1],length(ASMhum$ASMi))) #hum vs MR
        ageType = factor(rep(myAges[2],length(dataType)))#age = 1, ASM = 2
        modelType = factor(rep(myModels[mod],length(dataType))) #Integrated, stranding or tagging
        df_tmp = data.frame(dataType,modelType,ageType,age)
        df = rbind(df,df_tmp)
      }
      if(mod==3){
        ASMhum = drawASM_percent(sd[[ma]][[mo]]$value,
                                 sd[[ma]][[mo]]$cov,
                                 offset=0,
                                 nn=dataList$ni_mr,
                                 mu_rSSM = 92,
                                 sig_rSSM = 6.6)
        age = c(ASMhum$ASMi)
        dataType = factor(rep(myData[2],length(ASMhum$ASMi))) #hum vs MR
        ageType = factor(rep(myAges[2],length(dataType)))#age = 1, ASM = 2
        modelType = factor(rep(myModels[mod],length(dataType))) #Integrated, stranding or tagging
        df_tmp = data.frame(dataType,modelType,ageType,age)
        df = rbind(df,df_tmp)
      }
    }
  }
}

df = as.data.frame(df[is.na(df$age)==FALSE,])
names(df)

myLabels = data.frame(dataType = rep(unique(df$dataType), each= 2),
                      ageType = rep(unique(df$ageType), 2),
                      x = rep(60,4), 
                      y = c(0.08,0.08,0.14,0.14), 
                      labels=LETTERS[1:4])

library(viridis)
x = viridis(5)

p <- ggplot(df, aes(x=age, group = modelType)) +
  xlim(0,80) +
  xlab("Age (years)") +
  ylab("Smoothed density (proportion of simulated outcomes)") +
  # scale_fill_manual(values = x[c(1,2)])+
  scale_fill_manual(values = grey.colors(2,start = 0.0,end=0.9))+
  geom_density(aes(fill = modelType), alpha = 0.3) +
  facet_grid(dataType~ageType, 
             scales="free_y")+
  geom_text(data=myLabels, aes(label=labels, y = y, x = x),
          hjust=0, vjust=0,
          inherit.aes = FALSE)

p <- p + theme(axis.text.x=element_text(size=12),
               axis.text.y=element_text(size=12),
               strip.text.x = element_text(size = 12),
               strip.text.y = element_text(size = 12),
               axis.title=element_text(size=12))
p <- p + theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                         panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

p <- p + labs(fill = "Model type")

p <- p + theme(text = element_text(size = 18)) # this will change all text size 
# (except geom_text)
print(p)

# p <- p + geom_point(data=myLegend2, 
#                     x = rep(22,2),
#                     y = c(0.17,0.14)+0.02, 
#                     alpha = 0.6,
#                     fill=x[c(2,3)],
#                     shape = 22, size = 6,
#                     inherit.aes = FALSE)
# # # 
dev.off()
