# Script for making rep files from the posterior
# distributions in mcoutMSY.dat, based on supplied
# variable names and quantiles.

# Author: SDN Johnson, November 30, 2016

# NOTE:
# Set var1, var2, var1probs, var2probs and recDevYrs before
# running.

library(coda);library(MASS);library(chemometrics);library(mvtnorm)
source("mseRtools.r")

recDevYrs <- 10:52
recDevPad <- rep(0,recDevYrs[1]-1)

# read in mcmc output
mcmc    <- as.mcmc(read.table ("mcoutMSY.dat", header = TRUE) )
mcmcDF  <- as.data.frame(mcmc)

# What two variables will we use, and what percentiles?
var1 <- "h"
var2 <- "spawnB54"
form <- paste(var2,"~",var1,sep = "")

var2Mat <- as.matrix(mcmc[,c(var1,var2)])


v1v2    <- lm ( as.formula(form), data = mcmcDF)
mean2d <- apply(X = var2Mat, FUN = mean, MARGIN = 2)
v1v2cov <- cov(var2Mat)

probs   <- c(0.1,0.9)
var1quants  <- quantile (mcmc[,var1], probs = probs)
v1m         <- mean (mcmc[,var1])
v2m         <- coef(v1v2)[1] + v1m*coef(v1v2)[2]

var2quants  <-  quantile(mcmc[,var2], probs=probs)

var2meanline  <-  coef(v1v2)[1] + var1quants*coef(v1v2)[2]

pairs = list (  mh_loB=c(v1m,var2quants[1]),
                mh_mB=c(v1m,v2m),
                mh_hiB=c(v1m,var2quants[2]),
                loh_mB=c(var1quants[1],var2meanline[1]),
                hih_mB=c(var1quants[2],var2meanline[2]))

densPairs <- lapply (  FUN = dmvnorm, X = pairs, mean = c(v1m,v2m),
                      sigma=v1v2cov)

densPairs <- as.numeric(densPairs) / sum(as.numeric(densPairs))
names(densPairs) <- names(pairs)

# Plot testing points on the joint posterior and join with an 
# ellipse of the appropriate Mahalanobis distance from the mean.
# 57% chosen visually by trial and error and

par(mfrow = c(1,1), mar = c(4,4,1,2), oma = c(0,0,0,0))

drawMahal(  var2Mat,
            center=mean2d,
            covariance=v1v2cov,
            quantile=c(0.57),las=1,
            lwd=c(2),col="grey75",xlab="",
            ylab="2016 Spawning Biomass (kt)", pch=16)
  abline (a=coef(v1v2)[1],b=coef(v1v2)[2],col="black",lty=2,lwd=3)
  abline (v = var1quants, lty =3, col = "black", lwd =2)
  for(p in pairs)
  {
    points(x=p[1],y=p[2], pch=16, col = "red", cex = 2)
  }
  # panLab(x=0.05,y=0.9,txt="(a)")

dev.new()

par(mfrow = c(2,1), mar = c(1,2,1,1) )

plot(x = mcmcDF[,var1], y = mcmcDF[,"MSY"],las=1,col="grey75",
        xlab="",ylab="Maximum Sustainable Yield (kt)", pch=16,
        ylim = range( mcmcDF[,"MSY"]))
  abline (v = var1quants, lty =3, col = "black", lwd =2)
  panLab(x=0.05,y=0.9,txt="(b)")

plot(x = mcmcDF[,var1], y = mcmcDF[,"M_f"],las=1,col="grey75",
        xlab="Steepness",ylab="Female Mortality (/yr)", pch=16,
        ylim = c(0.08,0.1))
  abline (v = var1quants, lty =3, col = "black", lwd =2)
  panLab(x=0.05,y=0.9,txt="(b)")


    
  
# # get rows for var1 and var2 quantiles
# rowIdx <- numeric(length(pairs))
# for (i in 1:length(pairs))
# {
#   p <- pairs[[i]]
#   v1q <- p[1]
#   v2q <- p[2]
#   diff1 <- abs(mcmc[,var1] - v1q)
#   diff2 <- abs(mcmc[,var2] - v2q)

#   # Take total difference, this is what we'll minimise to get the row
#   # number
#   totalDiff <- diff1+diff2
#   rowIdx[i] <- which (totalDiff == min(totalDiff))
# }

# # name the row idx vector
# names(rowIdx) <- names(pairs)

# # names to pull entries from rows
# recDevNames     <- paste("recDev",recDevYrs, sep = "")
# tauIndexNames   <- paste("tauIndex",1:3,sep ="")
# tauAges_mNames  <- paste("tauAges_m",1:3,sep="")
# tauAges_fNames  <- paste("tauAges_f",1:3,sep="")
# tauReleaseNames <- paste("tauReleases",1:5,sep="")
# qNames          <- paste("qg",1:3,sep = "")

# # Now make a list of rep files
# reps <- vector(mode="list",length=length(rowIdx))
# origRep   <- lisread ("sableopmod.rep")
# tauAges <- vector(mode="list",length=length(rowIdx))
# for (i in 1:length(rowIdx))
# {
#   reps[[i]] <- origRep
#   reps[[i]]$SSB0        <- mcmc[rowIdx[i],"SSB0"]
#   reps[[i]]$h           <- mcmc[rowIdx[i],"h"]
#   reps[[i]]$M_m         <- mcmc[rowIdx[i],c("M_m")]
#   reps[[i]]$M_f         <- mcmc[rowIdx[i],c("M_f")]
#   reps[[i]]$recDevs     <- c(recDevPad,mcmc[rowIdx[i],recDevNames])
#   reps[[i]]$tauIndex    <- mcmc[rowIdx[i],tauIndexNames]
#   reps[[i]]$tauAge_m    <- mcmc[rowIdx[i],tauAges_mNames]
#   reps[[i]]$tauAge_f    <- mcmc[rowIdx[i],tauAges_fNames]
#   reps[[i]]$tauRel      <- mcmc[rowIdx[i],tauReleaseNames]
#   reps[[i]]$qg          <- mcmc[rowIdx[i],qNames]

#   # compute average aging error (between males and females),
#   # needed for mseR, as there is only one error term
#   tauAges[[i]] <- (mcmc[rowIdx[i],tauAges_mNames] + mcmc[rowIdx[i],tauAges_fNames])/2
# }

# names(reps) <- names(rowIdx)
# names(tauAges) <- names(rowIdx)

# # Write element i of list object x to the active file. Note that any
# # elements must be numeric (as this is the only thing ADMB can read in),
# # also sub-list structure is removed and elements of lists are written as 
# # vectors, probably to be read in to ADMB as ragged arrays
# # usage:
# # cat ( "## writeAMDB output file, created ", format(Sys.time(), "%y-%m-%d %H:%M:%S"), "\n", sep = "", file = activeFile )
# # lapply ( X = seq_along(x), FUN = writeADMB, x, activeFile)
# writeADMB <- function( i, x, activeFile ) {
#   # Write element name
#   cat( paste( "# ", names(x)[[i]], "\n", sep = "" ), file=activeFile, append=TRUE )
#   # Write element value
#   if( is.matrix(x[[i]]) ) {
#     # Matrices
#     write.table( x[[i]], file=activeFile, row.names=FALSE, col.names=FALSE, append=TRUE )
#   } else if( is.list(x[[i]]) ) {
#     # Lists
#     lapply( x[[i]], "\n", FUN=cat, file=activeFile, append=TRUE )
#   } else {
#     # Vectors/Scalars
#     cat( x[[i]], "\n", file=activeFile, append=TRUE )  
#   }
# }  # end writeADMB

# # Use writeADMB to rewrite the repfile
# # Write out the rep files
# for (i in 1:length(reps))
# { 
#   fileName <- paste(names(reps)[i],".rep",sep="")
#   cat ( "## ", names(reps)[i]," Rep File, created ", 
#         format(Sys.time(), "%y-%m-%d %H:%M:%S"), "\n", 
#         sep = "", file = fileName, append = FALSE )
#   lapply (  X = seq_along(reps[[i]]), FUN = writeADMB, x = reps[[i]], 
#             activeFile=fileName )

#   # Now print the necessary values to screen
#   cat("Scenario: ", names(reps)[i],"\n")
  
#   cat("B0 = ", mcmc[rowIdx[i],"SSB0"],"\n")
#   cat("h = ",mcmc[rowIdx[i],"h"],"\n")
#   cat("M = ",mcmc[rowIdx[i],c("M_m","M_f")],"\n")
#   cat("tauIdx = ", mcmc[rowIdx[i],tauIndexNames],"\n")
#   cat("tauAges = ", tauAges[[i]], "\n", sep = "")
#   cat("tauRel = ", mcmc[rowIdx[i],tauReleaseNames],"\n")
#   cat("qg = ", mcmc[rowIdx[i],qNames],"\n")
# }


