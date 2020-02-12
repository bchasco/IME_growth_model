# source("runIntegratedHumerusMR.r")
# source("runNonIntegratedMR.r")
# source("runNonIntegratedhumerus.r")

myFiles = c("Rep_IntegratedHumerusMR.Rdata",
            "Rep_NonIntegratedHumerus.Rdata",
            "Rep_NonIntegratedMR.Rdata")

#Labels for the different model runs
DataType = c("Integrated", "Humerus", "Mark-recapture")
RandomEfs = c("Eps_kij, Eps_ki, Eps_Linfi", "Eps_kij, Eps_Linfi","Eps_kij, Eps_ki","Eps_kij")
GrowthProcess = c("von Bertalanffy", "gompertz", "logistic")

#Output table
myTable = as.data.frame(matrix(NA,36,20))

icnt = 1
for(f in 1:length(myFiles)){
  load(myFiles[f]) #load the results from the different data integration experiments
  for(maps in 1:4){
    for(models in 1:3){
      if(is.na(rep[[maps]][[models]][1])) #Some models, particularly logistic were not estimable.
        myTable[icnt,] = c(DataType[f],
                           RandomEfs[maps],
                           GrowthProcess[models],rep("DNC",17))
      if(is.na(rep[[maps]][[models]][1])==FALSE){
        if(f==1){
          myTable[icnt,] = c(DataType[f],
                             RandomEfs[maps],
                             GrowthProcess[models],
                             round(fres[maps,models],1),
                             conv[maps,models],
                             maxgr[maps,models],
                             round(rep[[maps]][[models]]$Linf,1),
                             round(rep[[maps]][[models]]$k,3),
                             round(rep[[maps]][[models]]$L0,1),
                             round(rep[[maps]][[models]]$mu_a,1),
                             round(rep[[maps]][[models]]$mu_a_mr,1),
                             round(rep[[maps]][[models]]$beta1,2),
                             round(rep[[maps]][[models]]$sig_x,4),
                             round(rep[[maps]][[models]]$sig_L,4),
                             round(exp(rep[[maps]][[models]]$lnsig_bph),4),
                             round(rep[[maps]][[models]]$sig_A,2),
                             round(rep[[maps]][[models]]$sig_A_mr,2),
                             round(rep[[maps]][[models]]$sig_ki,2),
                             round(rep[[maps]][[models]]$sig_kij,2),
                             round(rep[[maps]][[models]]$sig_Linfi,2)
          )
        }
        if(f==2){
          myTable[icnt,] = c(DataType[f],RandomEfs[maps],GrowthProcess[models],
                             round(fres[maps,models],1),
                             conv[maps,models],
                             maxgr[maps,models],
                             round(rep[[maps]][[models]]$Linf,1),
                             round(rep[[maps]][[models]]$k,3),
                             round(rep[[maps]][[models]]$L0,1),
                             round(rep[[maps]][[models]]$mu_a,1),
                             "n.e.",
                             "n.e.",
                             round(rep[[maps]][[models]]$sig_x,4),
                             round(rep[[maps]][[models]]$sig_L,4),
                             "n.e.",
                             round(rep[[maps]][[models]]$sig_A,2),
                             "n.e.",
                             round(rep[[maps]][[models]]$sig_ki,2),
                             round(rep[[maps]][[models]]$sig_kij,2),
                             round(rep[[maps]][[models]]$sig_Linfi,2)
                             
          )
        }
        if(f==3){
          myTable[icnt,] = c(DataType[f],RandomEfs[maps],GrowthProcess[models],
                             round(fres[maps,models],1),
                             conv[maps,models],
                             maxgr[maps,models],
                             round(rep[[maps]][[models]]$Linf,1),
                             round(rep[[maps]][[models]]$k,3),
                             round(rep[[maps]][[models]]$L0,1),
                             "n.e.",
                             round(exp(rep[[maps]][[models]]$lnmu_a_mr),1),
                             "n.e.",
                             "n.e.",
                             round(rep[[maps]][[models]]$sig_L,4),
                             "n.e.",
                             "n.e.",
                             round(rep[[maps]][[models]]$sig_A_mr,2),
                             round(rep[[maps]][[models]]$sig_ki,2),
                             round(rep[[maps]][[models]]$sig_kij,2),
                             round(rep[[maps]][[models]]$sig_Linfi,2)
          )
        }
      }
      icnt = icnt + 1
    }
  }
}
write.table(myTable,
            "table_modelSummary.out",
            sep = "\t",
            row.names = FALSE, 
            col.names = FALSE, 
            quote = FALSE)
