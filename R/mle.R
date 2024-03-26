#' Maximum likelihood estimates for bivariate (tail) copula densities
#' 
#' @description
#' `mleBiTC` returns ML estimates for bivariate tail copula densities and
#' `mleBiCop` returns ML estimates for bivariate copula densities.
#' 
#' 
#' @param ft A log-likelihood function (see: \code{BiLoglik()})
#' @param family An integer; indicates the type of bivariate (tail) copula densities.
#' Possible tail copula families include:
#' 1=Husler-Reiss
#' 2=Negative Logistic
#' 3=Logistic
#' 4=Dirichlet
#' Possible copula families include:
#' 0=Independence
#' 1=Gaussian
#' 3=Clayton
#' 4=Gumbel
#' 5=Frank
#' 6=Joe
#' 13=Survival Clayton
#' 14=Survival Gumbel
#' 16=Survival Joe
#' @param data An \eqn{N\times 2} data matrix from the bivariate inverted Pareto samples with the second column being less than 1
#' @param range A numeric vector; indicates the range of parameters
#'
#' @return A numeric; ML estimate(s)
#' @export
mleBiTC <- function(ft, family, data, range){
  ML.out=stats::optimise(f = ft,family=family,x=data,maximum = F,interval = range)
  ML=ML.out$minimum
  return(ML)  
}
mleBiCop <- function(ft, family, data, range){
  ML.out=stats::optimise(f = ft,family=family,u=data,maximum = F,interval = range)
  ML=ML.out$minimum
  return(ML)  
}


