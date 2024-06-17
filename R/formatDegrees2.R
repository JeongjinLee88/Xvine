##  Internal function from graphicalExtremes package by Sebastian Engelke and Adrien S. Hitz and Nicola Gnecco and Manuel Hentschel.
formatDegrees2 <- function (decDeg, dirStrings, degString) 
{
  dms <- decDegToDegMinSec(decDeg)
  x <- paste0(dms[, 1], degString, dirStrings)
  if (!any(duplicated(x))) {
    return(x)
  }
  x <- paste0(dms[, 1], degString, " ", sprintf("%02d", dms[, 
                                                            2]), "'", dirStrings)
  if (!any(duplicated(x))) {
    return(x)
  }
  x <- paste0(dms[, 1], degString, " ", sprintf("%02d", dms[, 
                                                            2]), "' ", dms[, 3], "\"", dirStrings)
  return(x)
}