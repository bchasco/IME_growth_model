
myFiles = c("Rep_IntegratedHumerusMR.Rdata")
source("func_simulation.r") #Get the growth functions to make the code easier to read.
#Labels for the different model runs
RandomEfs = c("Eps_kij, Eps_ki, Eps_Linfi", "Eps_kij, Eps_Linfi","Eps_kij, Eps_ki","Eps_kij")
GrowthProcess = c("von Bertalanffy", "gompertz", "logistic")

library(mvtnorm)
source("func_quantAgeAtLen.r")
icnt = 1

x1 =c(0.1,0.5,0.1,0.5) 
x2 = c(0.5,0.9,0.5,0.9)
y1 = c(0.5,0.5,0.1,0.1)
y2 = c(0.9,0.9,0.5,0.5)

ymax = 115
xmax = 70
load(myFiles[1]) #load the results from the different data integration experiments
bestModel = which(fres == min(fres), arr.ind = TRUE)
map = bestModel[1]
model = bestModel[2]
if(is.na(sd[[map]][[model]])==FALSE){
      myq = quantAgeAtLen(sd[[map]][[model]]$value,sd[[map]][[model]]$cov,offset=0,dataList$ni)  
}      

# pdf("fig_4.pdf", width=6.2, height=6.2)
png("figure_4.png", width=620, height=620)

matplot(t(myq), type="l", lty=c(2,1,2), las=1, ylab="Age (years)", xlab = "SCL (cm)",
        lwd=3, col="black",
        yaxs="i",
        xaxs="i")
myq[,40]

dev.off()
