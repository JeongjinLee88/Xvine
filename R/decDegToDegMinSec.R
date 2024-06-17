##  Internal function from graphicalExtremes package by Sebastian Engelke and Adrien S. Hitz and Nicola Gnecco and Manuel Hentschel.
decDegToDegMinSec <- function (decDeg, asString = FALSE) 
{
  deg <- floor(decDeg)
  decMin <- (decDeg - deg) * 60
  min <- floor(decMin)
  decSec <- (decMin - min) * 60
  if (asString) {
    return(paste(deg, min, decSec))
  }
  return(cbind(deg, min, decSec))
}