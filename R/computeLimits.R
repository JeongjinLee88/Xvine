##  Internal function from graphicalExtremes package by Sebastian Engelke and Adrien S. Hitz and Nicola Gnecco and Manuel Hentschel.
computeLimits <- function (xData, yData, xyRatio = 1, convertLatLong = TRUE, stretch = 1) 
{
  xRange <- range(xData)
  yRange <- range(yData)
  xMid <- mean(xRange)
  yMid <- mean(yRange)
  xRadius <- diff(xRange)/2
  yRadius <- diff(yRange)/2
  if (is.null(xyRatio)) {
    xScale <- 1
    yScale <- 1
    xyRatio0 <- NULL
  }
  else {
    xyRatio0 <- xyRatio
    if (convertLatLong) {
      xyRatio <- xyRatio/cos(pi/180 * yMid)
    }
    xScale <- xyRatio/(xRadius/yRadius)
    yScale <- 1/xScale
    xScale <- max(1, xScale)
    yScale <- max(1, yScale)
  }
  return(list(xlim = xMid + xRadius * xScale * c(-1, 1) * stretch, 
              ylim = yMid + yRadius * yScale * c(-1, 1) * stretch, 
              xyRatio = xyRatio0))
}