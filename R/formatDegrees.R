##  Internal function from graphicalExtremes package by Sebastian Engelke and Adrien S. Hitz and Nicola Gnecco and Manuel Hentschel.
formatDegrees <- function (decDeg, direction = "NS", latex = TRUE) 
{
  dir <- 1 + (decDeg < 0)
  dirStrings <- substring(direction, dir, dir)
  if (latex) {
    degString <- "^{\\circ}"
    delim <- "$"
    dirStrings <- paste0("\\mathrm{", dirStrings, "}")
  }
  else {
    degString <- "°"
    delim <- ""
  }
  decDeg <- abs(decDeg)
  isNa <- is.na(decDeg)
  degStrings <- rep("", length(decDeg))
  degStrings[!isNa] <- paste0(delim, formatDegrees2(decDeg[!isNa], 
                                                    dirStrings[!isNa], degString), delim)
  return(degStrings)
}