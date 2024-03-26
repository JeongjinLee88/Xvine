#' Transform to multivariate Pareto samples via rank transformation
#' @description
#' `ParetoTransRank` transforms the original data to multivariate Pareto samples with suitable margins (the default: Pareto margins)
#' via rank transformation.
#' 
#' @details
#' Let \eqn{\boldsymbol(X)_{i}, i=1,\ldots,n} be independent samples from the unkown distribution \eqn{F}.
#' Transform to the uniform margin via rank transformation \eqn{\widehat{U}_{i,j}=1-(rank_{i,j}-0.5)/n}, where
#' \eqn{rank_{i,j}=\sum_{s=1}^{n}\mathbb{I}(X_{s,j}\le X_{i,j})} is the maximal rank of \eqn{X_{i,j}} among \eqn{X_{1,j},\ldots,X_{n,j}}.
#' The scaled points \eqn{(n/k)\widehat{U}_i} for \eqn{i=1,\ldots,n} such that \eqn{\min\widehat{U}_i<k/n} are pseudo-observations
#' from the inverted multivariate Pareto distribution. We set \eqn{k\in\{1,\ldots,n\}} 
#' such that \eqn{k} is large and \eqn{k/n} is small.
#' If we take the reciprocal of the rescaled samples, then we have Pareto scales.
#' If we take the negative log of the rescaled samples, then we have exponential scales.
#' 
#' @param data An \eqn{n\times d} data matrix
#' @param u_quan A numeric quantile, set a lower enough quantile \eqn{u_{quan}\in(0,1)}
#' @param scaleType A character specifying the type of scales for generated samples (the default is uniform scale: scaleType="U").
#' Other possible scales: 
#' * scaleType="P": Pareto scale
#' * scaleType="E": Exponential scale
#' 
#' @return An \eqn{N\times d} pseudo-observations for the inverted Pareto distribution (default)
#' @export
#'
#' @examples
#' StrMtx <- matrix(c(1, 1, 2, 2, 4,
#'              0, 2, 1, 3, 2,
#'              0, 0, 3, 1, 3,
#'              0, 0, 0, 4, 1,
#'              0, 0, 0, 0, 5),5,byrow = TRUE)
#' ParMtx <- matrix(c(0, 1.5, 2, 2.5, 2,
#'                   0, 0, 2, 2.5, 0.7,
#'                   0, 0, 0, 0.4, -0.3,
#'                   0, 0, 0, 0, 0.1,
#'                   0, 0, 0, 0, 0),5,byrow = TRUE)
#' FamMtx <- matrix(c(0, 1, 2, 3, 4,
#'                 0, 0, 3, 4, 1,
#'                 0, 0, 0, 3, 1,
#'                 0, 0, 0, 0, 1,
#'                 0, 0, 0, 0, 0),5,byrow = TRUE)
#' ## X-Vine specification (Vine structure, Bivariate parametric families, Parameters)
#' XVS=XVineSpec(M = StrMtx, Mmod = FamMtx, Mpar = ParMtx)
#' ##  Multivariate Pareto samples
#' Dat_P=ParetoSim(n = 5000, XVS = XVS) # Pareto scale
#' ##  Transform to pseudo-observations for inverted Pareto distribution
#' Dat_U=ParetoTransRank(data = Dat_P, u_quan = 0.2, scaleType = "U")
#' Dat_U
ParetoTransRank <- function(data, u_quan, scaleType="U")
{
  ##  Transform to uniform scale
  n <- nrow(data)
  k <- floor(n*u_quan)
  ##  Centered Rank trans to avoid boundary effects
  Uhat <- apply(data, 2, function(i) (n-rank(i)+0.5)/n)
  data_t <- (n/k)*Uhat[apply((n/k)*Uhat,1,min)<1,]
  if(scaleType=="U"){
    data_t=data_t
  }else if(scaleType=="P"){
    data_t=1/data_t
  }else if(scaleType=="E"){
    data_t=-log(data_t)
  }
  return(data_t)
}
