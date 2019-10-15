# Script to calculate reference points with MCMC output
# Author: SDN Johnson
# Date: Wednesday, November 23, 2016

# source tools
source("mseRtools.r")
source("calcRefPointsMCMC.R")
source("mseRrefPoints.r")
library ( parallel ); library(coda)

recDevYrs <- 10:54
recDevPad <- rep(0,recDevYrs[1]-1)

# Read in MCMC output
mcmc <- as.mcmc(read.table ("mcout.dat", header = TRUE) )



# also, run MSY calculations on MCMC output
# Read in sim control file
.CTLFILE <- "simCtlFile.txt"
ctlPars <- .readParFile( .CTLFILE )
ctlList <- .createList( ctlPars )
cat( "\nMSG (runMSE) Parameter list created.\n" )
nCores <- detectCores()-1
cl <- makeCluster ( nCores )
mcmcMSY <- parApply (cl, X = mcmc, MARGIN = 1, FUN = .calcRefPointsFromMCMC, 
                      parList = ctlList$opMod )
# mcmcMSY <- apply ( X = mcmc, MARGIN = 1, FUN = .calcRefPointsFromMCMC, 
                      # parList = ctlList$opMod )
mcmcMSY <- do.call("rbind", mcmcMSY)
stopCluster(cl)
write.table(x =  mcmcMSY, file = "mcoutMSY.dat",row.names=FALSE,quote=FALSE)

refPtPost <- mcmcMSY[,c("Bmsy","Fmsy","Umsy","legUmsy","Dmsy","MSY")]

# Now make an output table for the reference points
means <- c("mean",apply(X = refPtPost, MARGIN = 2, FUN = mean))
meds  <- c("median",apply(X = refPtPost, MARGIN = 2, FUN = median))
sd    <- c("sd",apply(X = refPtPost, MARGIN = 2, FUN = sd))

RefPts <- as.data.frame(rbind(means,meds,sd))
names(RefPts)[1] <- "Statistic"

write.csv(RefPts, file = "refPtsSumm.csv")


