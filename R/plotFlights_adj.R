#' Create figures using flight delays data
#' 
#' @description
#' 
#' `plotFlights_adj()` creates useful figures such as the flight connections from the flight delay dataset.
#' This function is a slight adjusted version of the function [plotFlights()] in `graphicalExtremes` package.
#' 
#'
#' @param airportIndices Character string; the indices of the airports 
#' @param airports_sel the set of airports
#' @param connections_sel A three column data frame by \code{flightCountMatrixToConnectionList()}
#' @param graph An optional \code{igraph::graph} object
#' @param edgeColors Vector or symmetric matrix (character or numeric) used as colors for edges/connections
#' @param plotAirports Logical; whether to plot the specified airports
#' @param plotConnections Logical; whether to plot the specified connections
#' @param returnGGPlot a \code{ggplot2::ggplot} object is return if `TRUE`
#' @param useAirportNFlights Logical; whether to change the size of the circles indicating airports
#' @param useConnectionNFlights Logical; whether to change the size of the edges indicating connections 
#' @param minNFlights Numeric; plot connections with at least this specified number of flights
#' @param map String or data.frame or NULL as a background map
#' @param vertexColors Vector with IATA codes used as colors for the vertices/airports
#' @param vertexShapes Vector with IATA codes used as shapes for the vertices/airports
#' @param xyRatio Numeric; approximate X-Y ratio of the area
#' @param clipMap Logical or numeric scalar; whether to ignore the map image when specifying the axis limits of the figure
#' @param useLatex Logical; whether to format numbers as latex code
#' @param edgeAlpha 0 or 1; the alpha value used in plotting edges/connections
#'
#' @return A `ggp` object
#' @export
#' @import graphicalExtremes
#' @references 
#' Engelke S, Hitz A, Gnecco N, Hentschel M (2023). _graphicalExtremes:
#' Statistical Methodology for Graphical Extreme Value Models_. R package
#' version 0.3.0, <https://github.com/sebastian-engelke/graphicalExtremes>
#'
plotFlights_adj <- function (airportIndices = NULL, airports_sel = NULL, connections_sel = NULL, 
                             graph = NULL, edgeColors, plotAirports = TRUE, plotConnections = TRUE, 
                             returnGGPlot = FALSE, useAirportNFlights = FALSE, useConnectionNFlights = FALSE, 
                             minNFlights = 0, map = "state", vertexColors = NULL, vertexShapes = NULL, 
                             xyRatio = NULL, clipMap = FALSE, useLatex = FALSE, edgeAlpha = 0.4) 
{
  computeLimits<-fitInInterval<-formatDegrees<-NULL
  ggplotAvailable <- requireNamespace("ggplot2")
  if (!ggplotAvailable) {
    stop("ggplot2 needs to be installed")
  }
  flights <- flights
  if (is.null(airports_sel)) {
    airports_sel <- flights$airports
  }
  if (!is.null(graph)) {
    vNames <- NULL
    if (!is.null(vNames) && is.null(airportIndices)) {
      airportIndices <- vNames
    }
  }
  if (is.null(map) || is.na(map) || nchar(map) == 0) {
    map <- NULL
  }
  if (!is.null(airportIndices)) {
    airports_sel <- airports_sel[airportIndices, ]
  }
  IATAS <- airports_sel[, "IATA"]
  stretchMap <- 1
  if (is.numeric(clipMap)) {
    stretchMap <- fitInInterval(1 * clipMap, 0, Inf)
    clipMap <- (clipMap > 0)
  }
  computeConnections <- is.null(connections_sel) && plotConnections
  computeNFlights <- is.null(airports_sel$nFlights) && useAirportNFlights && 
    plotAirports
  if (computeConnections || computeNFlights) {
    nFlightMat <- apply(flights$flightCounts[IATAS, IATAS, 
    ], c(1, 2), sum)
  }
  if (computeConnections) {
    connections_sel <- flightCountMatrixToConnectionList(nFlightMat)
  }
  if (computeNFlights) {
    airports_sel$nFlights <- rowSums(nFlightMat) + colSums(nFlightMat)
  }
  if (!is.null(graph)) {
    nVertices <- nAirports <- nrow(airports_sel) # changed
    if (nVertices != nAirports) {
      stop(sprintf("Number of vertices (%d) and number of selected airports (%d) are different!", 
                   nVertices, nAirports))
    }
  }
  aesVertexColor <- NULL
  if (!is.null(vertexColors)) {
    aesVertexColor <- "vertexColor"
    airports_sel$vertexColor <- vertexColors[IATAS]
  }
  aesVertexShape <- NULL
  if (!is.null(vertexShapes)) {
    aesVertexShape <- "vertexShape"
    airports_sel$vertexShape <- as.character(vertexShapes[IATAS])
  }
  if (plotConnections) {
    if (is.null(graph)) {
      ind <- ((connections_sel$departureAirport %in% IATAS | 
                 connections_sel$arrivalAirport %in% IATAS) & 
                connections_sel$nFlights >= minNFlights)
      connections_sel <- connections_sel[ind, ]
    }
    else {
      m <- graph
      connections_graph <- data.frame(matrix(IATAS[m], 
                                             ncol = 2))
      airportColNames <- c("departureAirport", "arrivalAirport")
      colnames(connections_graph) <- airportColNames
      rownames(connections_sel) <- paste0(connections_sel$departureAirport, 
                                          "_", connections_sel$arrivalAirport)
      rownames(connections_graph) <- paste0(connections_graph$departureAirport, 
                                            "_", connections_graph$arrivalAirport)
      if (useConnectionNFlights) {
        connections_graph$nFlights <- 0
        connections_graph$nFlights <- connections_sel[rownames(connections_graph), 
                                                      "nFlights"]
        connections_graph$nFlights[is.na(connections_graph$nFlights)] <- 1
      }
      connections_sel <- connections_graph
    }
    connections_sel$x0 <- airports_sel[connections_sel$departureAirport, 
                                       "Longitude"]
    connections_sel$y0 <- airports_sel[connections_sel$departureAirport, 
                                       "Latitude"]
    connections_sel$x1 <- airports_sel[connections_sel$arrivalAirport, 
                                       "Longitude"]
    connections_sel$y1 <- airports_sel[connections_sel$arrivalAirport, 
                                       "Latitude"]
  }
  aesEdgeColor <- NULL
  if (plotConnections && !is.null(edgeColors)) {
    aesEdgeColor <- "edgeColors"
    if (is.matrix(edgeColors)) {
      if (is.null(dimnames(edgeColors))) {
        dimnames(edgeColors) <- list(IATAS, IATAS)
      }
      ind <- cbind(connections_sel$departureAirport, connections_sel$arrivalAirport)
      connections_sel$edgeColors <- edgeColors[ind]
    }
    else if (is.vector(edgeColors)) {
      connections_sel$edgeColors <- edgeColors
    }
    else {
      stop("Argument `edgeColors` must be a vector with one entry per edge, or a matrix.")
    }
  }
  aesSizeNodes <- NULL
  aesSizeEdges <- NULL
  if (useAirportNFlights) {
    aesSizeNodes <- "nFlights"
  }
  if (useConnectionNFlights) {
    aesSizeEdges <- "nFlights"
  }
  ggp <- (ggplot2::ggplot() + ggplot2::xlab(NULL) + ggplot2::ylab(NULL) + 
            ggplot2::scale_x_continuous(labels = function(x) formatDegrees(x, 
                                                                                               "EW", useLatex)) + ggplot2::scale_y_continuous(labels = function(x) formatDegrees(x, 
                                                                                                                                                                                                     "NS", useLatex)))
  if (!is.null(map)) {
    dmap <- ggplot2::map_data(map)
    ggp <- ggp + ggplot2::geom_polygon(data = dmap, ggplot2::aes_string(x = "long", 
                                                                        y = "lat", group = "group"), color = "grey65", fill = "#f9f9f9", 
                                       size = 0.2)
  }
  if (!is.null(xyRatio) || (clipMap && !is.null(map))) {
    xData <- airports_sel$Longitude
    yData <- airports_sel$Latitude
    if (!clipMap && !is.null(map)) {
      m <- ggplot2::map_data(map)
      xData <- c(xData, range(m$long))
      yData <- c(yData, range(m$lat))
    }
    limits <- computeLimits(xData, yData, xyRatio = xyRatio, 
                                                stretch = stretchMap)
    ggp <- ggp + ggplot2::coord_cartesian(xlim = limits$xlim, 
                                          ylim = limits$ylim)
    if (!is.null(limits$xyRatio)) {
      ggp <- ggp + ggplot2::theme(aspect.ratio = 1/limits$xyRatio)
    }
  }
  if (plotAirports) {
    ggp <- ggp + ggplot2::geom_point(data = airports_sel, 
                                     ggplot2::aes_string(x = "Longitude", y = "Latitude", 
                                                         size = aesSizeNodes, col = aesVertexColor, shape = aesVertexShape), 
                                     na.rm = TRUE, alpha = 1)
  }
  if (plotConnections) {
    ggp <- ggp + ggplot2::geom_segment(data = connections_sel, 
                                       ggplot2::aes_string(x = "x0", xend = "x1", y = "y0", 
                                                           yend = "y1", size = aesSizeEdges, col = aesEdgeColor), alpha = edgeAlpha,col=edgeColors)+
      theme(legend.position = "none")
    
  }
  if (returnGGPlot) {
    return(ggp)
  }
  #return(invisible(NULL))
  return(list("ggp"=ggp))
}

