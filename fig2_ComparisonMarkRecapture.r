myPNG=TRUE
if(myPNG==TRUE)
  pdf("fig_ComparisonMarkRecapture.pdf", width=4.00, height=4.00)
par(mfrow=c(1,1))

myFiles = c("Rep_IntegratedHumerusMR.Rdata",
            "Rep_NonIntegratedHumerus.Rdata",
            "Rep_NonIntegratedMR.Rdata")


load(myFiles[3]) #load the results from the different data integration experiments
bestModel = which(fres == min(fres), arr.ind = TRUE)
map = bestModel[1]
model = bestModel[2]
rep2 = rep[[map]][[model]]

nij= aggregate(dataList$Lmr_ijID,list(i=dataList$Lmr_iID),max)[,2]
aveL = c(0,(dataList$Lobs[2:length(dataList$Lobs)]+dataList$Lobs[1:(length(dataList$Lobs)-1)]))/2           
dL = c(0,rep2$Lpred[2:length(rep2$Lpred)] - rep2$Lpred[1:(length(rep2$Lpred)-1)])[dataList$Lmr_ijID>0]

obsdLdT = c(0,(dataList$Lobs[2:length(dataList$Lobs)] - dataList$Lobs[1:(length(dataList$Lobs)-1)]))/(dataList$del_mrT)
obsdLdT = obsdLdT[dataList$Lmr_ijID>0]
dLdT = dL/(dataList$del_mrT)[dataList$Lmr_ijID>0]

L = aveL[dataList$Lmr_ijID>0]
lRange = seq(0,100,10)
uRange = lRange + 10
comp = read.table("AvensCompData.dat", header=TRUE)

uprSD = NA
lwrSD = NA
mean = NA
obsMean = NA
upr = NA
lwr = NA

mySpace = 0.1
myN = 3
par(mai=c(1.5,1,0.1,0.1))
for(i in 1:length(lRange))
{
  if(length(dLdT[L>=lRange[i] & L<uRange[i]])>0){
    mean[i] =   mean(dLdT[L>=lRange[i] & L<uRange[i]])
    obsMean[i] = mean(obsdLdT[L>=lRange[i] & L<uRange[i]])
    mysd = sd(dLdT[L>=lRange[i] & L<uRange[i]])
    uprSD[i] = mean[i] + 2*mysd
    lwrSD[i] = mean[i] - 2*mysd
    upr[i] = max(dLdT[L>=lRange[i] & L<uRange[i]])
    lwr[i] = min(dLdT[L>=lRange[i] & L<uRange[i]])
  }
  
  if(length(dLdT[L>=lRange[i] & L<uRange[i]])==0){
    mean[i] = NA
    uprSD[i] = NA
    lwrSD[i] = NA
    upr[i] = NA
    lwr[i] = NA
  }
}

maxMean_mr = na.omit(mean[mean==max(na.omit(mean))])
maxlR_mr = na.omit(lRange[mean==maxMean_mr])
maxuR_mr = na.omit(uRange[mean==maxMean_mr])

plot(L,
     obsdLdT,
     type="p",
     xaxs="i",
     yaxs="i",
     las=1,
     axes=FALSE,
     ylab="",
     xlab="",
     col="lightgrey",
     lwd=1,
     cex=1.2,
     xlim=c(0,110),
     ylim=c(-2,25))
abline(h=0, lty=2)
par(new=TRUE)
plot(1,
     bg="transparent",
     type="n",
     xaxs="i",
     yaxs="i",
     las=1,
     axes=FALSE,
     ylab="",
     xlab="",
     xlim=c(0,length(lRange)*myN),
     ylim=c(-2,25))

mtext("Length bins (SCL; cm)", side=1, line=5.5)
text(30,20,"B")

for(i in 1:length(lRange)){
  axis(1, at = (i-1)*myN + 1, labels=paste0(lRange[i]," - <", uRange[i]),
       las=2)
}
axis(2, las=1)
for(i in 1:length(lRange)){
  rect(myN*(i-1)+mySpace,lwrSD[i],myN*(i-1)+1-mySpace,uprSD[i],
       # col=rgb(0.5,0.1,0.1,0.3),
       # border=rgb(0.5,0.1,0.1,0.3))
      col=grey(0.3,alpha=0.8),
      border='black', lwd=2)
segments(myN*(i-1)+mySpace,mean[i],myN*(i-1)+1-mySpace,mean[i],
           col="black",
           border=rgb(1,0,1),
           lwd=2)
}

for(i in 1:nrow(comp))
{
  if(comp$author[i]=="Bjorndal"){
    os = 1
    myCol = grey(0.1,alpha=0.3)

    rect(myN*(comp$inc[i]-1)+mySpace+os,comp$mean[i]-2*comp$SD[i],myN*(comp$inc[i]-1)+1-mySpace+os,comp$mean[i]+2*comp$SD[i],
         col=myCol,
         border='black', lwd=2)
    segments(myN*(comp$inc[i]-1)+mySpace+os,comp$mean[i],myN*(comp$inc[i]-1)+1-mySpace+os,comp$mean[i],
             col="black",
             lwd=2,
             border=rgb(0.1,1,0.5))
  }
  if(comp$author[i]=="McNeill"){
    os = 2
    myCol = grey(0,alpha=0)
    rect(myN*(comp$inc[i]-1)+mySpace+os,comp$lwr[i],myN*(comp$inc[i]-1)+1-mySpace+os,comp$upr[i],
         col=myCol,
         border='black', lwd=2)
    segments(myN*(comp$inc[i]-1)+mySpace+os,comp$mean[i],myN*(comp$inc[i]-1)+1-mySpace+os,comp$mean[i],
             col="black",
             lwd=2,
             border=rgb(1,0.5,0.5))
  } 
}

box()
legend(1,25,legend=c("Current study", "Bjorndal et al. (2013)", "Braun-McNeill et al. (2008)"),
       pch=15,
       bty="n",
       bg="transparent",
       pt.cex = 2.5,
       col=c(grey(0.3,alpha=0.8),grey(0.1,0.3),grey(0,alpha=0)))

legend(1,25,legend=c("Current study", "Bjorndal et al. (2013)", "Braun-McNeill et al. (2008)"),
       pch=0,
       bty="n",
       bg="transparent",
       pt.cex = 2.5,
       col=c(grey(0.3,alpha=0.7),rgb(0.5,0.5,0.1, 0.3),'black'))


par(new=TRUE)
plot(1,
     type="n",
     xaxs="i",
     yaxs="i",
     las=1,
     axes=FALSE,
     ylab="Growth increment (cm/year)",
     xlab="",
     xlim=c(0,110),
     ylim=c(-2,25))

if(myPNG==TRUE)
  dev.off()
