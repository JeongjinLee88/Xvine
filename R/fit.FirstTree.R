#' Fit the first maximum spanning tree to multivariate inverted-Pareto data
#'
#' @description 
#' `fit.FirstTree()` fits the first maximum spanning tree to multivariate inverted-Pareto data set.
#'  The function takes a pair of variables for each edge and variables are passed to the function 'tcSelect'.
#'  After selecting the family of bivariate tail copula models, store the output of the fitted model in the list of MST
#'  , including parameter estimates, model classes, fitted models, and pseudo-observations.
#'  
#' @param MST A list of maximum spanning trees.
#' @param data.invPa A \eqn{n\times d} data matrix from multivariate inverse-Pareto distribution.
#' @param tcfamset Numeric vector; specifies the class of bivariate tail copula models that users consider.
#' @param selectioncrit Character string, indicating the selection criteria for the class of bivariate tail copula models.
#' @param weights Logical; whether weights should be assigned to observations when missing values exist.
#' @param si Numeric; a tuning parameter for mBIC (\eqn{\si_0=0.9}; default).
#'
#' @return a nested list object:
#' * the maximum spanning tree with a node set and edge set.
#' * the fitted tail copula model with selected classes, parameter estimates, log-likelihood values, AIC, and pseudo-observations.
#' @export
#'
fit.FirstTree <- function (MST, data.invPa, tcfamset, si, selectioncrit, weights = NA)
{
  d <- nrow(MST$E$nums)
  pc.data <- lapply(1:d, function(i) NULL)
  for (i in 1:d) { # i=1 2 3
    a <- MST$E$nums[i, ]
    pc.data[[i]]$zr1 <- data.invPa[, a[1]]
    pc.data[[i]]$zr2 <- data.invPa[, a[2]]
    if (is.null(MST$V$names[a[1]])) {
      MST$E$Copula.CondName.1[i] <- a[1]
    }
    else {
      MST$E$Copula.CondName.1[i] <- MST$V$names[a[1]]
    }
    if (is.null(MST$V$names[a[2]])) {
      MST$E$Copula.CondName.2[i] <- a[2]
    }
    else {
      MST$E$Copula.CondName.2[i] <- MST$V$names[a[2]]
    }
    if (is.null(MST$V$names[a[1]]) || is.null(MST$V$names[a[2]])) {
      MST$E$Copula.Name[i] <- paste(a[1], a[2], sep = " , ")
    }
    else {
      MST$E$Copula.Name[i] <- paste(MST$V$names[a[1]], 
                                    MST$V$names[a[2]], sep = " , ")
    }
  }
 
  pc.fits <- lapply(pc.data, tcSelect, familyset = tcfamset, si=si,
                    selectioncrit = selectioncrit)
  for (i in 1:d) {
    MST$E$Copula.param[[i]] <- c(pc.fits[[i]]$par, pc.fits[[i]]$par2)
    MST$E$Copula.type[i] <- pc.fits[[i]]$family
    MST$E$AIC[i] <- pc.fits[[i]]$AIC
    MST$E$BIC[i] <- pc.fits[[i]]$BIC
    MST$E$mBIC[i] <- pc.fits[[i]]$mBIC
    MST$E$EffectSize[i] <- pc.fits[[i]]$Eff_k
    MST$E$fits[[i]] <- pc.fits[[i]]
    MST$E$Copula.CondData.1[i] <- list(pc.fits[[i]]$CondOn.1)
    MST$E$Copula.CondData.2[i] <- list(pc.fits[[i]]$CondOn.2)
  }
  MST
}

tcSelect <- function (tcParams, ...) 
{
  return(tcFamSel(tcParams$zr1, tcParams$zr2, 
                   ...))
}
