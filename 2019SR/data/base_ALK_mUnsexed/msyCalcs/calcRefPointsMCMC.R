# ARK (03-Oct-16)
#
# This code extracted from mseRmod (2010).
# Inputs:
# parList - needed parameters for doing MSY or other ref point calcs such as
#           the fixed parameters (growth, maturity).
# mcmc    - the ADMB mcmc output file.

# Buried in mseRmod are these two lines, which calls .calcRefPointsFromMCMC.
#           parList <- fitTracker[[i]]$rep$refPoints
#           mcmcStats <- .calcRefPointsFromMCMC( parList, mcmc )
#
# Then, within .calcRefPointsFromMCMC was a call to calcRefPoints.
# Note the use of "ref" and "deref" from the "ref" package (pass by reference
# rather than pass by value) to avoid making copies.
#
# Probably this code could be made simpler when a GUI is not being used.
# 


.calcRefPointsFromMCMC <- function( obj, parList )
{
  library(ref)
  source("mseRrefPoints.r")
  source("mseRglobals.r")
  source("mseRtools.r")

  objNames <- names(obj)
  if ( is.null(dim(obj)) ) obj <- matrix(obj, nrow = 1, ncol = length(obj))
  colnames(obj) <- objNames

  # Change these... ***
  # resultNames <- c( "B0","h","M","MSY","Bmsy","Fmsy","Umsy","legUmsy","Dmsy",
  #                   "projSSB","projLegalB","legHR.T","slegHR.T",
  #                   "Uadj","legHarv","genTime" )
  resultNames <- c (  "Bmsy", "Fmsy", "Umsy", "legUmsy", "Dmsy" )

  result <- as.data.frame( matrix( NA, nrow=nrow(obj),ncol=length(resultNames) ) )
  names(result) <- resultNames
 
  cat( "\nMSG (.calcRefPointsFromMCMC) Calculating reference points for each draw.\n" )
          
  nDraws <- nrow(obj)
  nPrint <- 20

  t1 <- proc.time()

    
  for ( i in 1:nDraws )
  {
    # Extract the parameters from mcmc object and replace.
    
    # Unfished biomass.
    if ( any(objNames=="SSB0" ) )
      parList$B0 <- obj[i,"SSB0"]
    
    # Steepness.          
    if ( any(objNames=="h" ) )
      parList$rSteepness <- obj[i,"h" ]

    # Recruitment deviations.          
    idx <- grep( "recDev", objNames )
    if ( length( idx ) > 0 )
      parList$estRecDevs  <- obj[i,idx]
    
    # Natural mortality.  
    if ( any(objNames=="M_f") )
      parList$M[2] <- obj[ i, "M_f" ]

    if ( any(objNames=="M_m") )
      parList$M[1] <- obj[ i, "M_m" ]
      
            
    # Selectivities - this needs to grab the elements of each and assemble vectors
    # for calcRefPoints.
    idNames <- paste( "alpha_g",c(1:5),sep="" )
    alpha_g <- unlist( obj[ i,idNames ] )

    idNames <- paste( "beta_g",c(1:5),sep="" )
    beta_g <- unlist( obj[ i,idNames ] )


    selType           <- parList$repFile$selType
    parList$selType   <- selType
    parList$L50Cg1    <- alpha_g
    parList$L95Cg1    <- beta_g

    # Fix asymptotic pars
    parList$L95Cg1[which(selType == 1)] <- alpha_g[which(selType == 1)]
    parList$L50Cg1[which(selType == 1)] <- alpha_g[which(selType == 1)] - beta_g[which(selType == 1)]


    # Calculate reference points and stuff into result.
    refPoints <- calcRefPoints( parList )
            
    # result$B0[i]       <- obj$B0[i]
    # result$h[i]        <- obj$rSteepness[i]
    # result$M[i]        <- obj$M[i]
    result$MSY[i]      <- refPoints$yieldFmsy
    result$Bmsy[i]     <- refPoints$ssbFmsy
    result$Fmsy[i]     <- refPoints$Fmsy
    result$Umsy[i]     <- refPoints$Umsy
    result$legUmsy[i]  <- refPoints$legalHRFmsy
    result$Dmsy[i]     <- result$Bmsy[i] / obj[,"SSB0"]
    
    # Change these... ***
    # result$projSSB[i]    <- obj$projSSB[i]
    # result$projLegalB[i] <- obj$projLegalB[i]
    # #result$projDep[i]    <- obj$projDep[i]
    # result$legHR.T[i]    <- obj$legalHR.T[i]
    # result$slegHR.T[i]   <- obj$sublegalHR.T[i]
    # result$genTime[i]    <- refPoints$genTime
    
    # Calculate harvest control rule adjusted U for the i-th draw of the chain.
    # val <- list()
    # val$B0                <- obj$B0[i]
    
    # Change these... ***
    # val$SSB               <- obj$projSSB[i]
    # val$legalB            <- obj$projLegalB[i]
    # val$sublegalBt        <- obj$projSublegalB[i]
    # val$refPoints$ssbFmsy <- result$Bmsy[i]
    # #val$refPoints$Umsy    <- result$Umsy[i]
    # val$refPoints$Umsy    <- result$legUmsy[i]    
    
    # ruleList <- .calcDfoRule( val )  
   
    # result$Uadj[i] <- ruleList$precautionaryRemRate
    # result$legHarv[i] <- result$Uadj[i] * obj$projLegalB[i]   
    
    #cat( "\nDraw number = ", i, "\n" )
    
    # Print a message every X calls and for the first 5 to reassure progress.
    if ( i %% nPrint==0 | i <=5 )
    {
      t2 <- proc.time()
      cat( "\nMSG (.calcRefPointsFromMCMC) Completed ",i," of ",nDraws," draws.\n\n" )
      cat( "\nMSG (.calcRefPoints) took ",t2-t1," seconds for ",i,"calls.\n" )
      cat( "\nMSG (.calcRefPoints) Average time per call is ",(t2-t1)[1]/i," seconds.\n\n" )      
      print( result[i,] )
    }
  }

  retVal <- cbind( obj, result )
  retVal
}     # .calcRefPointsFromMCMC function
