#Read in the three data sets
#The is the humerus data
hx = read.table("lagDiameter2.dat", header=TRUE)
#Mark recapture data
tx = read.table("tagDat.dat", header=TRUE)
#Humerus data as it related to the body proportion hypothesis
sx = read.table("strandingData.dat", header=TRUE)

set.seed(200)

#all of the MR turtles with greater than minObs
t2 = sort(unique(tx$ChascoID[tx$RecapID>minObs]))
tij = rep(0,nrow(tx))
tij[tx$ChascoID%in%t2] = 1
tij[tx$RecapID==1] = 0
tx$tij = tij

#A lot of turtle show no growth, convergence is a lot easier if you jitter the data a little
#for those turtles with no growth.
for(i in 2:nrow(tx))
{
  if(tx$RecapID[i]==tx$RecapID[i-1] &  tx$SCL[i]==tx$SCL[i-1])
    tx$SCL[i] = tx$SCL[i]# + rnorm(1,0,0.01)
}

#Add a little jitter. This is necessary because a lot of the larger diameters have several
#repeated sizes the cause problems when estimating the transient effects.
hx = hx[,2:ncol(hx)]
for(i in 2:nrow(hx))
{
  if(hx$x_iID[i]==hx$x_iID[i-1] &  hx$D[i]==hx$D[i-1])
    hx$D[i] = hx$D[i]# + rnorm(1,0,0.0001)
}

#Create years associated with each LAG
sx[,7] = as.numeric(as.character(sx[,7]))
sx[,8] = as.numeric(as.character(sx[,8]))
yLAG = sx[1,7] - sx[1,8]:1 
for(i in 2:dim(sx)[1])
{
  yLAG = c(yLAG,sx[i,7] - sx[i,8]:1)
}
true_yLAG = yLAG
yLAG = yLAG - min(yLAG)

#Add the yLAG column and neritic or pelagic state based on diameter
hx = cbind(hx,yLAG)
hx$x_stateID = rep(0,dim(hx)[1])
hx$x_stateID[hx$D>20] = 1

#all of the stranded turtles with greater than min observations
h2 = sort(unique(hx$x_iID[hx$LAGID>(minObs-1)]))
hij = rep(0,nrow(hx))
hij[hx$x_iID%in%h2] = 1
hij[hx$LAGID==0] = 0
hx$hij = hij

#Put all of the information for the raw data into a data list
dataList = list(ni = max(hx[,1]), #Number of standed turtles
                nia = sum(sx$a0), #Number of stranded turtles without an annulus 
                a0 = sx$a0, #Does the turtle possess an annulus
                growthModel = 1, #1 = VB, 2 = gompertz, 3 = logistic
                rho_mod = 2, #DEPRECATED
                nLAG = sx$nLAG, #number of LAG for each turtle
                Xobs = hx$D, #Observed diameters of each LAG
                Xmax = sx$humerus, #Maximum diameter of humerus for BPH model
                SCLmax = sx$SCL, #maximum carapace length for stranded turtle for BPH model
                del_xT = hx$del_xT, #number of years between each LAG, almost always 1
                x_iID = hx$x_iID -1, #vector of the stranded turtle IDs,s starting with 0.
                x_ijID = hx$LAGID, #vector of LAG IDs, range from 0 to 33 depending on the number of LAG
                x_yrID = hx$yLAG, #The year each LAG formed, this is DEPRECATED because we no longer have a temporal vairable
                x_stateID = hx$x_stateID, #Neritic or pelagic, DEPRECATED because we no longer have this affect in the model
                hij = hx$hij, #Whether a transient random effect should be included in the 
                hum_obs_model = 1, #log-normal, normal, cv - DEPRECATED
                mr_obs_model = 1, #log-normal, normal, cv - DEPRECATED
                sex_ijID = rep(sx$sexID,sx$nLAG), #model no longer includes sex
                Lobs = tx$SCL, #Length at mark or recapture
                tij = tx$tij, #Whether a transient effect should be included in the k_ij parameter
                n_mr = length(tx$ChascoID), #this is DEPRECATED
                del_mrT = tx$delT/365.25, #Number of years between mark and recapture
                Lmr_iID = tx$ChascoID - 1, #Tagged turtle ID, starting at 0
                Lmr_ijID = tx$RecapID - 1, #Recapture ID, starting at 0
                Lvec = seq(75,95,1), #DEPRECATED  
                ni_mr = length(unique(tx$ChascoID))#Number of tagged turtles
)


