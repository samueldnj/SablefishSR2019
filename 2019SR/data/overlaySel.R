source(here::here("mseRtools.R"))
source("base_ALK_mUnsexed/sableOpMod.r")

cond2016 <- lisread("2016CondRep/sableopmod.rep")
cond2019 <- lisread("base_ALK_mUnsexed/sableopmod.rep")


makeSelLen <- function( obj, gearNum )
{
  
  # Initial values.
  alpha_g1  <- obj$alpha_g1[gearNum]
  beta_g1   <- obj$beta_g1[gearNum] 

  # Get size limit and sel type
  sizeLim <- obj$sizeLim
  selType <- obj$selType

  # Calculate curves from parameters
  len <- seq(from = 32, to = 75, length = 100)


  Sel_l <- numeric(length = length(len))

  if(selType[gearNum] == 1)
  {
    Sel_l <- 1 / (1 + exp(-log(19) * (len - beta_g1) / (alpha_g1 - beta_g1)))
  }
  if(selType[gearNum] == 2)
  {
    Sel_l <- exp(-(0.5 * (len - alpha_g1)/beta_g1)^2 )
  }

  if(selType[gearNum] == 3)
  {
    Sel_l <- len ^(alpha_g1 - 1) * exp(-len/beta_g1)
    Sel_l <- Sel_l / max(Sel_l)
  }

  out <- cbind(len, Sel_l)

  return( out )
}


xRange  <- c(32, max( cond2019$lenAge_m, cond2019$lenAge_f ))
yRange  <- c(0,1)
sizeLim <- 55

# Make selectivities
trawl2016 <- makeSelLen( obj = cond2016, gearNum = 3)
trawl2019 <- makeSelLen( obj = cond2019, gearNum = 3)


plot( x = xRange, y = yRange, type = "n",
      xlab = "Length (cm)", ylab = "Selectivity at length" )
  lines( x = trawl2016[,1], y = trawl2016[,2], col = "grey40",
          lwd = 3, lty = 2 )
  lines( x = trawl2019[,1], y = trawl2019[,2], col = "black",
          lwd = 3, lty = 1 )
  abline( v = sizeLim, col = "red", lty = 2 )
  legend( x = "topright",
          legend = c( "2016 trawl selectivity (Normal)",
                      "2018 trawl Selectivity (Gamma)",
                      "Size limit"),
          lty = c(2,1,2),
          lwd = c(3,3,1),
          col = c("grey40","black","red"))


