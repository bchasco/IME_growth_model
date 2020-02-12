myFiles = c("Rep_IntegratedHumerusMR.Rdata","Rep_NonIntegratedHumerus.Rdata", "Rep_NonIntegratedMR.Rdata")

# pdf("fig_5.pdf", width=7.0, height=8.5)
png("figure_5.png", width=700, height=850)

load(myFiles[1])
loc = which(fres == min(fres), arr.ind = TRUE)
xx = rep[[loc[1]]][[loc[2]]]
par(fig=c(0.1,0.6,0.1,0.4), mai=c(0,0,0,0))
plot(xx$a_ij,
     xx$L_ij,
     ylim=c(10,110),
     xlim=c(0,80),
     xaxs="i",
     las=1,
     # col=rgb(0.5,0.1,0.5,0.3),
     col=grey(0.7,alpha=0.9),
     yaxs="i",
     axes=FALSE)
axis(1)
axis(2,las=1)
box()
text(70,100,"C")

mtext("Age (years)", side=1, line=3)
mtext("SCL (cm)", side=2, line=2.5)
points(xx$a_ij_mr,
       xx$Lpred,
       # col=rgb(0.1,0.5,0.5,0.3))
        col=grey(0.1,0.3))
par(fig=c(0.6,0.8,0.1,0.4), new=TRUE,mai=c(0,0,0,0))
hx = hist(xx$L_ij,breaks=0:120,plot=FALSE)
barplot(hx$counts/max(hx$counts),
        horiz = TRUE,
        yaxs="i",
        col=grey(0.7,alpha=0.9),
        border=grey(0.7,alpha=0.9),
        # col=rgb(0.5,0.1,0.5,0.3),
        # border=rgb(0.5,0.1,0.5,0.3),
        space=0,
        axes=FALSE,
        xaxs="i")
par(fig=c(0.6,0.8,0.1,0.4), mai=c(0,0,0,0),new=TRUE)
hx = hist(xx$Lpred,breaks=0:120,plot=FALSE)
barplot(hx$counts/max(hx$counts),
        horiz = TRUE,
        yaxs="i",
        col=grey(0.1,0.3),
        # col=rgb(0.1,0.5,0.5,0.3),
        border=grey(0.1,0.3),
        space=0,
        axes=FALSE,
        xaxs="i")


load(myFiles[2])
loc = which(fres == min(fres), arr.ind = TRUE)
xx = rep[[loc[1]]][[loc[2]]]
par(fig=c(0.1,0.6,0.4,0.7), mai=c(0,0,0,0), new=TRUE)
plot(xx$a_ij,
     xx$L_ij,
     ylim=c(10,110),
     xlim=c(0,80),
     xaxs="i",
     las=1,
     col=grey(0.7,0.9),
        yaxs="i",
     axes=FALSE)
axis(1,labels=FALSE)
axis(2,las=1)
box()
text(70,100,"B")

mtext("SCL (cm)", side=2, line=2.5)
par(fig=c(0.6,0.8,0.4,0.7), new=TRUE,mai=c(0,0,0,0))
hx = hist(xx$L_ij,breaks=0:120,plot=FALSE)
barplot(hx$counts/max(hx$counts),
        horiz = TRUE,
        yaxs="i",
        col=grey(0.7,0.9),
        border=grey(0.7,0.9),
        # col=rgb(0.5,0.1,0.5,0.3),
        # border=rgb(0.5,0.1,0.5,0.3),
        space=0,
        axes=FALSE,
        xaxs="i")

load(myFiles[3])
loc = which(fres == min(fres), arr.ind = TRUE)
xx = rep[[loc[1]]][[loc[2]]]
par(fig=c(0.1,0.6,0.7,1), mai=c(0,0,0,0), new=TRUE)
plot(xx$a_ij_mr,
     xx$Lpred,
     ylim=c(10,110),
     xlim=c(0,80),
     xaxs="i",
     las=1,
     col=grey(0.1,0.3),
     # col=rgb(0.1,0.5,0.5,0.3),
     yaxs="i",
     axes=FALSE)
axis(1,labels=FALSE)
axis(2,las=1)
box()
text(70,100,"A")

mtext("SCL (cm)", side=2, line=2.5)
par(fig=c(0.6,0.8,0.7,1), new=TRUE,mai=c(0,0,0,0))
hx = hist(xx$Lpred,breaks=10:110,plot=FALSE)
barplot(hx$counts/max(hx$counts),
        horiz = TRUE,
        yaxs="i",
        col=grey(0.1,0.3),
        border=grey(0.1,0.3),
        lwd=0.1,
        # col=rgb(0.1,0.5,0.5,0.3),
        # border=rgb(0.1,0.5,0.5,0.3),
        space=0,
        axes=FALSE,
        xaxs="i")

dev.off()
