##  Internal function from graphicalExtremes package by Sebastian Engelke and Adrien S. Hitz and Nicola Gnecco and Manuel Hentschel.
fitInInterval <- function (x, xMin = -Inf, xMax = Inf) 
{
  if (any(xMax < xMin)) {
    stop("Make sure that xMax>=xMin!")
  }
  x <- pmax(x, xMin)
  x <- pmin(x, xMax)
  return(x)
}
