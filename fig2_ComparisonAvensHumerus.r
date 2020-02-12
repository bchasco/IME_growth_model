myPNG=TRUE

if(myPNG==TRUE)
  png("fig_ComparisonAvensHumerus.png", width=4.00, height=4.00)
par(mfrow=c(1,1))
load(myFiles[1]) #load the results from the different data integration experiments
bestModel = which(fres == min(fres), arr.ind = TRUE)
map = bestModel[1]
model = bestModel[2]
rep2 = rep[[map]][[model]]


transF = function(x, b1=rep2$beta1)
{
  return(((x)/b1))
}
obsdDdT = c(0,(transF(dataList$Xobs[2:length(dataList$Xobs)]) - transF(dataList$Xobs[1:(length(dataList$Xobs)-1)])))/dataList$del_xT

dLdT = c(0,rep2$L_ij[2:length(rep2$L_ij)] - rep2$L_ij[1:(length(rep2$L_ij)-1)])[dataList$x_ijID!=0] 
L = transF(dataList$Xobs)[dataList$x_ijID>0]

obsMean = NA
obsSd = NA
lRange = seq(0,100,10)
uRange = lRange + 10
comp = read.table("AvensCompData.dat", header=TRUE)

uprSD = NA
lwrSD = NA
mean = NA
upr = NA
lwr = NA

mySpace = 0.2
myN = 2
par(mai=c(1.5,1,0.1,0.1))
for(i in 1:length(lRange))
{
  if(length(dLdT[L>=lRange[i] & L<uRange[i]])>0){
    #i = 5
    obsMean[i] = mean((obsdDdT[dataList$x_ijID>0])[L>=lRange[i] & L<uRange[i]])
    obsSd[i] = sd((obsdDdT[dataList$x_ijID>0])[L>=lRange[i] & L<uRange[i]])
    mean[i] =   mean(na.omit(dLdT[L>=lRange[i] & L<uRange[i]]))
    sddLdT = sd(na.omit(dLdT[L>=lRange[i] & L<uRange[i]]))
    uprSD[i] = mean[i] + 2*sddLdT
    lwrSD[i] = mean[i] - 2*sddLdT
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
maxMean_hum = na.omit(mean[mean==max(na.omit(mean))])
maxlR_hum = na.omit(lRange[mean==maxMean_hum])
maxuR_hum = na.omit(uRange[mean==maxMean_hum])

plot(L,
     obsdDdT[dataList$x_ijID>0],
     xaxs="i",
     yaxs="i",
     axes=FALSE,
     ylab="",
     xlab="",
     ylim=c(-2,25),
     col="lightgrey",
     cex=1.2,
     xlim=c(0,110))
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
text(20,20,"A")

for(i in 1:length(lRange)){
  axis(1, at = (i-1)*myN + 1, labels=paste0(lRange[i]," - <", uRange[i]),
       las=2)
}
axis(2, las=1)

for(i in 1:length(lRange)){
  rect(myN*(i-1)+mySpace,lwrSD[i],myN*(i-1)+1-mySpace,uprSD[i],
       # myCol = rgb(0.5,0.1,0.1, 0.3)
       col=grey(0.1,alpha = 0.1),
       border="black",
       # border=rgb(0.1,0.5,0.1,0.5),
       lwd=2)
  segments(myN*(i-1)+mySpace,mean[i],myN*(i-1)+1-mySpace,mean[i],
           col="black",
           lwd=3,
           border=rgb(1,0,1,0.5))
}

for(i in 1:nrow(comp))
{
  if(comp$author[i]=="Avens"){
    os = 1
    # myCol = rgb(0.5,0.1,0.1, 0.3)
    myCol = grey(0.3, alpha=0.7)
    rect(myN*(comp$inc[i]-1)+mySpace+os,comp$mean[i]-2*comp$SD[i],myN*(comp$inc[i]-1)+1-mySpace+os,comp$mean[i]+2*comp$SD[i],
         col=myCol,
         border="black",
         lwd=2)
    segments(myN*(comp$inc[i]-1)+mySpace+os,comp$mean[i],myN*(comp$inc[i]-1)+1-mySpace+os,comp$mean[i],
             col="black",
             lwd=3,
             border=myCol)
  }
  if(comp$author[i]=="Bjorndal"){
    os = 3  
    myCol = rgb(0.5,0.5,1)
  } 
  if(comp$author[i]=="McNeill"){
    os = 4
    myCol = rgb(1,0.5,0.5)
  } 
  
  
}
box()
legend(1,25,legend=c("Current study", "Avens et al. (2013, 2015)"),
       pch=15,
       pt.cex = 2.5,
       bty="n",
       bg="transparent",
       col=c(grey(0.3,alpha=0.7),grey(0.1,alpha=0.1)))
abline(h=0, lty=2)
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
